module mod_io_meshdiag
    ! fesom.mesh.diag.zarr writer (M9 Stage 1) — a UGRID-1.0, xarray/ushow-readable Zarr analog of
    ! FESOM2's fesom.mesh.diag.nc (io_mesh_info.F90:write_mesh_info, the transcription oracle).
    !
    ! Transforms (cited vs io_mesh_info.F90):
    !   lon/lat          geo_coord_nod2D(1:2,:) / rad          (rad -> deg; :409 rbuffer/rad)
    !   nz/nz1           -zbar / -Z                            (CF positive-down; :336/:337)
    !   zbar_e/n_bottom  -zbar_e_bot / -zbar(nlevels_nod2D)    (:394/:400; n_bottom computed)
    !   face_nodes       myList_nod2D(elem2D_nodes)            (global node ids, :432)
    !   edge_nodes       myList_nod2D(edges)                   (:445)
    !   face_edges       myList_edge2D(elem_edges)             (FESOM2 writes RAW=local :456; we write
    !                                                           GLOBAL so the store is partition-INDEP;
    !                                                           identical at 1-rank where local==global)
    !   face_links       elem_neighbors/=self ? myList_elem2D : -999 (FESOM2 RAW :469; global here)
    !   edge_face_links  edge_tri(2)/=edge_tri(1) ? myList_elem2D : -999 (:488; slot 1 always)
    !
    ! Hybrid meshes (plan docs/plans/completed/2026-10-08-hybrid-meshes.md, Task 14): every per-face array
    ! has 4 columns (dim n4) in the slot order of the model (clockwise; edge k joins vertices k and
    ! k+1). A triangle has only 3: its 4th column holds the array's UGRID _FillValue (the Zarr
    ! fill_value) -- face_nodes / face_edges -1, face_links -999 (FESOM2's missing-neighbour value,
    ! used for both) -- while the model itself repeats the 1st vertex/edge/neighbour (user
    ! decision). gradient_sca_x/y: the 4 Green-Gauss coefficients (x = slots 1:4, y = 5:8 of
    ! gradient_sca), 0 in a triangle's 4th column. xarray decodes a declared _FillValue: with
    ! mask_and_scale (the default) these integer arrays come back as float with NaN.
    !   nod_in_elem2D    num>=i ? myList_elem2D : 0            (:508; N=max num; fill 0)
    !   nod_part/elem_part  mype stamp                         (:379/:387; partition-DEPENDENT by design)
    !
    ! Canonical global-id ordering + distributed-chunk-writers via mod_io_decomp => the store is
    ! partition-INDEPENDENT (dist_2 ≡ dist_8) with no rank-0 gather. Store-create ordering: rank 0
    ! creates the store + defines ALL arrays -> barrier -> writers write their chunks -> rank 0
    ! consolidates. Optional-`partit` (absent or npes==1 ⇒ 1-rank identity).
    use mpi
    use, intrinsic :: iso_fortran_env, only: real64
    use mod_precision,   only: WP
    use mod_constants,   only: rad
    use mod_mesh,        only: t_mesh
    use mod_partit,      only: t_partit
    use mod_part_bounds, only: is_multirank, local_dims
    use mod_io_zarr
    use mod_io_decomp
    implicit none
    private
    public :: meshdiag_write

    integer, parameter :: DEFAULT_CHUNK = 500000   ! horizontal chunk (clamped to entity N per decomp)

contains

    subroutine meshdiag_write(path, mesh, partit, chunk_horiz, n_writers)
        character(len=*), intent(in)           :: path
        type(t_mesh),     intent(in)           :: mesh
        type(t_partit),   intent(in), optional :: partit
        integer,          intent(in), optional :: chunk_horiz, n_writers

        type(t_zarr_store) :: store
        type(t_io_decomp)  :: Dn, De, Dg            ! node / elem / edge decomps
        type(t_zarr_attrs) :: gat
        logical :: mr
        integer :: mype, comm, ierr, C, nw, nl, N_max, i, c3, nv, nb
        integer :: nNodO, nNodL, nEdgeO, nEdgeL, nElemO, nElemL, nElemF
        ! zarr array handles (init on all ranks; defined on rank 0; used by writers)
        type(t_zarr_array) :: a_lon, a_lat, a_nz, a_nz1, a_elem_area, a_nlev_n, a_nlev_e, &
            a_nie_num, a_nod_part, a_elem_part, a_zbe_bot, a_zbn_bot, a_nod_area, &
            a_face_nodes, a_edge_nodes, a_edge_face_links, a_face_edges, a_face_links, &
            a_nod_in_elem, a_edge_cross, a_grad_x, a_grad_y, a_fmesh
        ! local field buffers
        real(WP), allocatable :: r1(:), r3(:,:)
        integer,  allocatable :: i1(:), i3(:,:)
        real(WP), allocatable :: depth(:)

        mr   = is_multirank(partit)
        mype = 0; comm = MPI_COMM_SELF
        if (present(partit)) then; mype = partit%mype; comm = partit%MPI_COMM_FESOM; end if
        C = DEFAULT_CHUNK; if (present(chunk_horiz)) C = chunk_horiz
        nw = 0;            if (present(n_writers))   nw = n_writers
        call read_env_int('FESOM3_CHUNK_HORIZ', C)
        call read_env_int('FESOM3_N_WRITERS',  nw)
        nl = mesh%nl

        call local_dims(mesh, partit, nNodO, nNodL, nEdgeO, nEdgeL, nElemO, nElemL, nElemF)

        ! decomps (chunk size clamped to N inside; read back D%C for the matching Zarr chunk)
        call decomp_init_entity(Dn, C, nw, DECOMP_NODE, mesh, partit)
        call decomp_init_entity(De, C, nw, DECOMP_ELEM, mesh, partit)
        call decomp_init_entity(Dg, C, nw, DECOMP_EDGE, mesh, partit)

        ! N_max for nod_in_elem2D (max adjacency over all ranks)
        N_max = 0
        do i = 1, nNodO; N_max = max(N_max, mesh%nod_in_elem2D_num(i)); end do
        if (mr) call MPI_Allreduce(MPI_IN_PLACE, N_max, 1, MPI_INTEGER, MPI_MAX, comm, ierr)

        ! ---------------- define phase (rank 0): store + all arrays ----------------
        store%path = path
        if (mype == 0) then
            call zarr_attrs_init(gat)
            call zattr_str(gat, 'Conventions', 'UGRID-1.0')
            call zattr_str(gat, 'description', 'FESOM3 mesh diagnostics (Zarr v2, M9)')
            call zarr_create_store(store, path, gat)
            call add_topology(store)
        end if
        call def2(store, a_lon,  'lon',  Dn, '<f8', 'nod2', 'longitude', 'degrees_east',  std='longitude')
        call def2(store, a_lat,  'lat',  Dn, '<f8', 'nod2', 'latitude',  'degrees_north', std='latitude')
        call def_global1(store, a_nz,  'nz',  nl,   '<f8', 'nz',  'depth of levels', 'meters', positive='down')
        call def_global1(store, a_nz1, 'nz1', nl-1, '<f8', 'nz1', 'depth of layers', '',       positive='down')
        call def2(store, a_elem_area, 'elem_area', De, '<f8', 'elem', 'element areas', '')
        call def2(store, a_nlev_n, 'nlevels_nod2D', Dn, '<i4', 'nod2', 'number of levels below nodes', '')
        call def2(store, a_nlev_e, 'nlevels',       De, '<i4', 'elem', 'number of levels below elements', '')
        call def2(store, a_nie_num,'nod_in_elem2D_num', Dn, '<i4', 'nod2', 'number of elements containing the node', '')
        call def2(store, a_nod_part, 'nod_part', Dn, '<i4', 'nod2', 'nodal partitioning at the cold start', '')
        call def2(store, a_elem_part,'elem_part', De, '<i4', 'elem', 'element partitioning at the cold start', '')
        call def2(store, a_zbe_bot, 'zbar_e_bottom', De, '<f8', 'elem', 'element bottom dep', '')
        call def2(store, a_zbn_bot, 'zbar_n_bottom', Dn, '<f8', 'nod2', 'nodal bottom depth', '')
        call def3(store, a_nod_area, 'nod_area', nl, Dn, '<f8', 'nz', 'nod2', 'nodal areas', '')
        call def3(store, a_face_nodes, 'face_nodes', 4, De, '<i4', 'n4', 'elem', &
            'Maps every face to its corner nodes (a triangle: 3, then _FillValue).', '', &
            cf='face_node_connectivity', sidx=.true., fillv=-1.0_real64, hasfill=.true.)
        call def3(store, a_edge_nodes, 'edge_nodes', 2, Dg, '<i4', 'n2', 'edg_n', &
            'Maps every edge to the two nodes that it connects', '', cf='edge_node_connectivity', sidx=.true.)
        call def3(store, a_face_edges, 'face_edges', 4, De, '<i4', 'n4', 'elem', &
            'Maps every face to its edges (edge k joins nodes k and k+1; a triangle: 3, then _FillValue).', &
            '', cf='face_edge_connectivity', sidx=.true., fillv=-1.0_real64, hasfill=.true.)
        call def3(store, a_face_links, 'face_links', 4, De, '<i4', 'n4', 'elem', &
            'neighbor faces for faces (across edge k; missing neighbor and a triangle 4th: _FillValue)', &
            '', cf='face_face_connectivity', sidx=.true., fillv=-999.0_real64, hasfill=.true.)
        call def3(store, a_edge_face_links, 'edge_face_links', 2, Dg, '<i4', 'n2', 'edg_n', &
            'neighbor faces for edges', '', cf='edge_face_connectivity')
        call def3(store, a_nod_in_elem, 'nod_in_elem2D', N_max, Dn, '<i4', 'N', 'nod2', &
            'elements containing the node', '')
        call def3(store, a_edge_cross, 'edge_cross_dxdy', 4, Dg, '<f8', 'n4', 'edg_n', 'edge cross distancess', '')
        call def3(store, a_grad_x, 'gradient_sca_x', 4, De, '<f8', 'n4', 'elem', &
            'x component of a gradient at nodes of an element', '')
        call def3(store, a_grad_y, 'gradient_sca_y', 4, De, '<f8', 'n4', 'elem', &
            'y component of a gradient at nodes of an element', '')

        if (mr) call MPI_Barrier(comm, ierr)

        ! ---------------- write phase ----------------
        ! global (replicated) 1-D coords + topology — rank 0 only
        if (mype == 0) then
            allocate(depth(nl))
            depth(1:nl) = real(-mesh%zbar(1:nl), WP);  call zarr_write_whole(store, a_nz,  depth(1:nl))
            depth(1:nl-1) = real(-mesh%Z(1:nl-1), WP); call zarr_write_whole(store, a_nz1, depth(1:nl-1))
            deallocate(depth)
            call zarr_write_whole(store, a_fmesh, [0])
        end if

        ! lon / lat (node, real)
        allocate(r1(max(1,nNodO)))
        do i = 1, nNodO; r1(i) = real(mesh%geo_coord_nod2D(1,i)/rad, WP); end do
        call wr2_r(store, a_lon, Dn, r1, 0.0_WP)
        do i = 1, nNodO; r1(i) = real(mesh%geo_coord_nod2D(2,i)/rad, WP); end do
        call wr2_r(store, a_lat, Dn, r1, 0.0_WP)
        ! zbar_n_bottom = -zbar(nlevels_nod2D(n))
        do i = 1, nNodO; r1(i) = real(-mesh%zbar(mesh%nlevels_nod2D(i)), WP); end do
        call wr2_r(store, a_zbn_bot, Dn, r1, 0.0_WP)
        deallocate(r1)

        ! node integer scalars
        allocate(i1(max(1,nNodO)))
        do i = 1, nNodO; i1(i) = mesh%nlevels_nod2D(i);    end do; call wr2_i(store, a_nlev_n,  Dn, i1, 0)
        do i = 1, nNodO; i1(i) = mesh%nod_in_elem2D_num(i);end do; call wr2_i(store, a_nie_num, Dn, i1, 0)
        do i = 1, nNodO; i1(i) = mype;                     end do; call wr2_i(store, a_nod_part, Dn, i1, 0)
        deallocate(i1)

        ! elem real / integer scalars
        allocate(r1(max(1,nElemO)))
        do i = 1, nElemO; r1(i) = real(mesh%elem_area(i), WP);  end do; call wr2_r(store, a_elem_area, De, r1, 0.0_WP)
        do i = 1, nElemO; r1(i) = real(-mesh%zbar_e_bot(i), WP);end do; call wr2_r(store, a_zbe_bot,   De, r1, 0.0_WP)
        deallocate(r1)
        allocate(i1(max(1,nElemO)))
        do i = 1, nElemO; i1(i) = mesh%nlevels(i); end do; call wr2_i(store, a_nlev_e,   De, i1, 0)
        do i = 1, nElemO; i1(i) = mype;            end do; call wr2_i(store, a_elem_part, De, i1, 0)
        deallocate(i1)

        ! nod_area (node, nl levels, real)
        allocate(r3(nl, max(1,nNodO)))
        ! FESOM3 bottom at vertices: mesh%area is 1-D (depth-independent). The diag file
        ! keeps its (nl, nod) schema, so broadcast the column value over the WET levels and
        ! leave the rest zero -- same shape, and still a faithful picture of the cell.
        r3 = 0.0_WP
        do i = 1, nNodO
            r3(mesh%ulevels_nod2D(i):mesh%nlevels_nod2D(i)-1, i) = real(mesh%area(i), WP)
        end do
        call wr3_r(store, a_nod_area, Dn, r3, nl, 0.0_WP)
        deallocate(r3)

        ! face_nodes / face_edges / face_links + gradient_sca_x/y (elem, 4 columns). Columns beyond
        ! the element's nv vertices get the array's fill value (a triangle's 4th).
        allocate(i3(4, max(1,nElemO)))
        do i = 1, nElemO
            nv = mesh%elem2D_nnodes(i)
            do c3 = 1, 4
                if (c3 <= nv) then; i3(c3,i) = gid_n(mesh%elem2D_nodes(c3,i))
                else;               i3(c3,i) = -1; end if
            end do
        end do
        call wr3_i(store, a_face_nodes, De, i3, 4, -1)
        do i = 1, nElemO
            nv = mesh%elem2D_nnodes(i)
            do c3 = 1, 4
                if (c3 <= nv) then; i3(c3,i) = gid_g(mesh%elem_edges(c3,i))
                else;               i3(c3,i) = -1; end if
            end do
        end do
        call wr3_i(store, a_face_edges, De, i3, 4, -1)
        do i = 1, nElemO
            nv = mesh%elem2D_nnodes(i)
            do c3 = 1, 4
                nb = mesh%elem_neighbors(c3,i)
                ! a boundary side stores the element itself (no zero index)
                if (c3 <= nv .and. nb /= i) then; i3(c3,i) = gid_e(nb)
                else;                             i3(c3,i) = -999; end if
            end do
        end do
        call wr3_i(store, a_face_links, De, i3, 4, -999)
        deallocate(i3)
        allocate(r3(4, max(1,nElemO)))
        do i = 1, nElemO; r3(1:4, i) = real(mesh%gradient_sca(1:4, i), WP); end do
        call wr3_r(store, a_grad_x, De, r3, 4, 0.0_WP)
        do i = 1, nElemO; r3(1:4, i) = real(mesh%gradient_sca(5:8, i), WP); end do
        call wr3_r(store, a_grad_y, De, r3, 4, 0.0_WP)
        deallocate(r3)

        ! edge_nodes / edge_face_links (edge, 2 levels) + edge_cross_dxdy (edge, 4 levels, real)
        allocate(i3(2, max(1,nEdgeO)))
        do i = 1, nEdgeO; do c3 = 1, 2; i3(c3,i) = gid_n(mesh%edges(c3,i)); end do; end do
        call wr3_i(store, a_edge_nodes, Dg, i3, 2, 0)
        do i = 1, nEdgeO
            do c3 = 1, 2
                ! a boundary edge stores its one element twice (edge_tri(2) == edge_tri(1))
                if (c3 == 1 .or. mesh%edge_tri(2,i) /= mesh%edge_tri(1,i)) then
                    i3(c3,i) = gid_e(mesh%edge_tri(c3,i))
                else
                    i3(c3,i) = -999
                end if
            end do
        end do
        call wr3_i(store, a_edge_face_links, Dg, i3, 2, -999)
        deallocate(i3)
        allocate(r3(4, max(1,nEdgeO)))
        do i = 1, nEdgeO; r3(1:4, i) = real(mesh%edge_cross_dxdy(1:4, i), WP); end do
        call wr3_r(store, a_edge_cross, Dg, r3, 4, 0.0_WP)
        deallocate(r3)

        ! nod_in_elem2D (node, N_max levels)
        allocate(i3(max(1,N_max), max(1,nNodO)))
        do i = 1, nNodO
            do c3 = 1, N_max
                if (mesh%nod_in_elem2D_num(i) >= c3) then; i3(c3,i) = gid_e(mesh%nod_in_elem2D(c3,i))
                else;                                      i3(c3,i) = 0; end if
            end do
        end do
        call wr3_i(store, a_nod_in_elem, Dn, i3, N_max, 0)
        deallocate(i3)

        if (mype == 0) call zarr_consolidate(store)
        if (mr) call MPI_Barrier(comm, ierr)

    contains

        ! ---- local->global id helpers (identity at 1-rank) ----
        integer function gid_n(l); integer, intent(in) :: l
            if (mr) then; gid_n = partit%myList_nod2D(l); else; gid_n = l; end if; end function
        integer function gid_e(l); integer, intent(in) :: l
            if (mr) then; gid_e = partit%myList_elem2D(l); else; gid_e = l; end if; end function
        integer function gid_g(l); integer, intent(in) :: l
            if (mr) then; gid_g = partit%myList_edge2D(l); else; gid_g = l; end if; end function

        ! ---- define helpers (init always; .zarray/.zattrs on rank 0) ----
        subroutine def2(st, arr, name, D, dtype, dim1, lname, units, std, fillv, hasfill)
            type(t_zarr_store), intent(inout) :: st
            type(t_zarr_array), intent(out)   :: arr
            character(len=*),   intent(in)    :: name, dtype, dim1, lname, units
            type(t_io_decomp),  intent(in)    :: D
            character(len=*),   intent(in), optional :: std
            real(real64),       intent(in), optional :: fillv
            logical,            intent(in), optional :: hasfill
            type(t_zarr_attrs) :: at
            real(real64) :: fv; logical :: hf
            ! fill_value: null by default so xarray never masks a valid 0 (e.g. nz surface, nod_part=0)
            fv = 0.0_real64; hf = .false.
            if (present(fillv)) fv = fillv
            if (present(hasfill)) hf = hasfill
            call zarr_array_init(arr, name, [D%N], [D%C], dtype, fill=fv, has_fill=hf)
            if (mype == 0) then
                call zarr_attrs_init(at)
                call zattr_str_arr(at, '_ARRAY_DIMENSIONS', [character(len=8) :: dim1])
                call zattr_str(at, 'long_name', lname)
                if (len_trim(units) > 0) call zattr_str(at, 'units', units)
                if (present(std)) call zattr_str(at, 'standard_name', std)
                call zarr_define_array(st, arr, at)
            end if
        end subroutine def2

        subroutine def3(st, arr, name, nlev, D, dtype, vdim, hdim, lname, units, cf, sidx, fillv, hasfill)
            type(t_zarr_store), intent(inout) :: st
            type(t_zarr_array), intent(out)   :: arr
            character(len=*),   intent(in)    :: name, dtype, vdim, hdim, lname, units
            integer,            intent(in)    :: nlev
            type(t_io_decomp),  intent(in)    :: D
            character(len=*),   intent(in), optional :: cf
            logical,            intent(in), optional :: sidx
            real(real64),       intent(in), optional :: fillv
            logical,            intent(in), optional :: hasfill
            type(t_zarr_attrs) :: at
            real(real64) :: fv; logical :: hf
            fv = 0.0_real64; hf = .false.
            if (present(fillv)) fv = fillv
            if (present(hasfill)) hf = hasfill
            call zarr_array_init(arr, name, [nlev, D%N], [nlev, D%C], dtype, fill=fv, has_fill=hf)
            if (mype == 0) then
                call zarr_attrs_init(at)
                call zattr_str_arr(at, '_ARRAY_DIMENSIONS', [character(len=8) :: vdim, hdim])
                call zattr_str(at, 'long_name', lname)
                if (len_trim(units) > 0) call zattr_str(at, 'units', units)
                if (present(cf)) call zattr_str(at, 'cf_role', cf)
                if (present(sidx)) then
                    if (sidx) call zattr_int(at, 'start_index', 1)
                end if
                call zarr_define_array(st, arr, at)
            end if
        end subroutine def3

        subroutine def_global1(st, arr, name, n, dtype, dim1, lname, units, positive)
            type(t_zarr_store), intent(inout) :: st
            type(t_zarr_array), intent(out)   :: arr
            character(len=*),   intent(in)    :: name, dtype, dim1, lname, units
            integer,            intent(in)    :: n
            character(len=*),   intent(in), optional :: positive
            type(t_zarr_attrs) :: at
            call zarr_array_init(arr, name, [n], [n], dtype, has_fill=.false.)
            if (mype == 0) then
                call zarr_attrs_init(at)
                call zattr_str_arr(at, '_ARRAY_DIMENSIONS', [character(len=8) :: dim1])
                call zattr_str(at, 'long_name', lname)
                if (len_trim(units) > 0) call zattr_str(at, 'units', units)
                if (present(positive)) call zattr_str(at, 'positive', positive)
                call zarr_define_array(st, arr, at)
            end if
        end subroutine def_global1

        subroutine add_topology(st)
            type(t_zarr_store), intent(inout) :: st
            type(t_zarr_attrs) :: at
            call zarr_array_init(a_fmesh, 'fesom_mesh', [1], [1], '<i4')
            call zarr_attrs_init(at)
            call zattr_str_arr(at, '_ARRAY_DIMENSIONS', [character(len=4) :: 'info'])
            call zattr_str(at, 'cf_role', 'mesh_topology')
            call zattr_str(at, 'long_name', 'Topology data of 2D unstructured mesh')
            call zattr_int(at, 'topology_dimension', 2)
            call zattr_str(at, 'node_coordinates', 'lon lat')
            call zattr_str(at, 'face_node_connectivity', 'face_nodes')
            call zattr_str(at, 'edge_node_connectivity', 'edge_nodes')
            call zattr_str(at, 'edge_face_connectivity', 'edge_face_links')
            call zattr_str(at, 'face_edge_connectivity', 'face_edges')
            call zattr_str(at, 'face_face_connectivity', 'face_links')
            call zarr_define_array(st, a_fmesh, at)
        end subroutine add_topology

        ! ---- write helpers (redistribute -> writer writes its chunks) ----
        subroutine wr2_r(st, arr, D, fld, fillv)
            type(t_zarr_store), intent(in) :: st
            type(t_zarr_array), intent(in) :: arr
            type(t_io_decomp),  intent(in) :: D
            real(WP),           intent(in) :: fld(:), fillv
            real(WP), allocatable :: buf(:)
            integer :: c, lo
            allocate(buf(max(1, D%w_nbuf)))
            call decomp_redistribute(D, fld, buf, fillv)
            do c = D%w_first_chunk, D%w_last_chunk
                lo = (c - D%w_first_chunk)*D%C + 1
                call zarr_write_chunk(st, arr, [c], buf(lo:lo+D%C-1))
            end do
        end subroutine wr2_r

        subroutine wr2_i(st, arr, D, fld, fillv)
            type(t_zarr_store), intent(in) :: st
            type(t_zarr_array), intent(in) :: arr
            type(t_io_decomp),  intent(in) :: D
            integer,            intent(in) :: fld(:), fillv
            integer, allocatable :: buf(:)
            integer :: c, lo
            allocate(buf(max(1, D%w_nbuf)))
            call decomp_redistribute(D, fld, buf, fillv)
            do c = D%w_first_chunk, D%w_last_chunk
                lo = (c - D%w_first_chunk)*D%C + 1
                call zarr_write_chunk(st, arr, [c], buf(lo:lo+D%C-1))
            end do
        end subroutine wr2_i

        subroutine wr3_r(st, arr, D, fld, nlev, fillv)
            type(t_zarr_store), intent(in) :: st
            type(t_zarr_array), intent(in) :: arr
            type(t_io_decomp),  intent(in) :: D
            integer,            intent(in) :: nlev
            real(WP),           intent(in) :: fld(:,:), fillv
            real(WP), allocatable :: buf(:,:)
            integer :: c, lo
            allocate(buf(nlev, max(1, D%w_nbuf)))
            call decomp_redistribute(D, fld, buf, fillv)
            do c = D%w_first_chunk, D%w_last_chunk
                lo = (c - D%w_first_chunk)*D%C + 1
                call zarr_write_chunk(st, arr, [0, c], buf(1:nlev, lo:lo+D%C-1))
            end do
        end subroutine wr3_r

        subroutine wr3_i(st, arr, D, fld, nlev, fillv)
            type(t_zarr_store), intent(in) :: st
            type(t_zarr_array), intent(in) :: arr
            type(t_io_decomp),  intent(in) :: D
            integer,            intent(in) :: nlev, fillv
            integer,            intent(in) :: fld(:,:)
            integer, allocatable :: buf(:,:)
            integer :: c, lo
            allocate(buf(nlev, max(1, D%w_nbuf)))
            call decomp_redistribute(D, fld, buf, fillv)
            do c = D%w_first_chunk, D%w_last_chunk
                lo = (c - D%w_first_chunk)*D%C + 1
                call zarr_write_chunk(st, arr, [0, c], buf(1:nlev, lo:lo+D%C-1))
            end do
        end subroutine wr3_i

    end subroutine meshdiag_write

    subroutine read_env_int(name, val)
        character(len=*), intent(in)    :: name
        integer,          intent(inout) :: val
        character(len=64) :: buf
        integer :: ios, tmp
        call get_environment_variable(name, buf, status=ios)
        if (ios == 0 .and. len_trim(buf) > 0) then
            read(buf, *, iostat=ios) tmp
            if (ios == 0) val = tmp
        end if
    end subroutine read_env_int

end module mod_io_meshdiag
