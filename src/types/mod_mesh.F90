module mod_mesh
    ! Static mesh type (decision D5). Structure from tracer_dwarf MOD_MESH
    ! (== FESOM2 v2.7.3 MOD_MESH), trimmed to v1 scope: cavity, iceberg (Z_3d_n_ib)
    ! and OASIS (lump2d_*/ind_*) fields are dropped.
    !
    ! Element-arity generalization (D2/D5): FESOM2 names are kept (elem2D_nodes,
    ! nod_in_elem2D, ...) so kernels transcribe verbatim, but elem2D_nodes /
    ! elem_edges / elem_neighbors are allocated (MAX_NV, elem2D) and a per-element
    ! vertex count elem2D_nnodes(:) is added; gradient_sca is (2*MAX_NV, elem2D).
    ! Hybrid meshes (mod_mesh_read): a triangle stores its 1st vertex again as the 4th
    ! (elem2D_nodes(4,e) == elem2D_nodes(1,e) <=> elem2D_nnodes(e) == 3), a quad 4
    ! distinct vertices, clockwise. Conventions without zero indices: elem_edges(k,e) is the
    ! edge between slots k and k+1 (1-2, 2-3, 3-4, 4-1; a triangle's 4th = its 1st),
    ! elem_neighbors(k,e) the element across it (missing = e itself; triangle 4th = 1st),
    ! edge_tri(2,ed) == edge_tri(1,ed) on a boundary edge (edge_cross_dxdy(3:4,ed) = 0),
    ! gradient_sca x in 1:4, y in 5:8 (a triangle's slots 4 and 8 are 0). In the
    ! multi-rank mesh coord_nod2D, geo_coord_nod2D, bc_index_nod2D, depth and
    ! nlevels_nod2D also hold the double-halo nodes
    ! (myDim+eDim+1 .. +eXDim_nod2D, when dist_N/format.out exists).
    !
    ! Geometry/coordinate arrays are MP (mesh precision, D3); == WP at the DP anchor.
    use mod_precision, only: WP, MP
    use, intrinsic :: iso_fortran_env, only: int32
#if !defined(USE_HALF_PRECISION)
    use mod_binary_arrays, only: write_bin_array, read_bin_array
#endif
    implicit none
    save

    integer, parameter :: MAX_NV       = 4    ! max element vertices (3=tri, 4=quad)
    integer, parameter :: MAX_ADJACENT = 32   ! max elements adjacent to a node

    type sparse_matrix
        integer :: nza = 0
        integer :: dim = 0
        real(kind=MP),  allocatable, dimension(:) :: values
        integer(int32), allocatable, dimension(:) :: colind
        integer(int32), allocatable, dimension(:) :: rowptr
        integer(int32), allocatable, dimension(:) :: colind_loc
        integer(int32), allocatable, dimension(:) :: rowptr_loc
        real(kind=MP),  allocatable, dimension(:) :: pr_values   ! inverse diagonal (Jacobi preconditioner), owned rows
    end type sparse_matrix

    type t_mesh
        ! ---- horizontal counts ----
        integer :: nod2D = 0
        integer :: elem2D = 0
        integer :: edge2D = 0, edge2D_in = 0
        real(kind=MP) :: ocean_area = 0.0_MP
        real(kind=MP) :: ocean_areawithcav = 0.0_MP

        ! ---- coordinates ----
        real(kind=MP), allocatable, dimension(:,:) :: coord_nod2D      ! (2,nod2D) rotated
        real(kind=MP), allocatable, dimension(:,:) :: geo_coord_nod2D  ! (2,nod2D) geographic

        ! ---- connectivity (arity-general: leading dim MAX_NV) ----
        integer, allocatable, dimension(:,:) :: elem2D_nodes    ! (MAX_NV,elem2D)
        integer, allocatable, dimension(:)   :: elem2D_nnodes   ! (elem2D) vertices per element
        integer, allocatable, dimension(:,:) :: edges           ! (2,edge2D)
        integer, allocatable, dimension(:,:) :: edge_tri        ! (2,edge2D) elems left/right of edge
        integer, allocatable, dimension(:,:) :: elem_edges      ! (MAX_NV,elem2D)
        integer, allocatable, dimension(:,:) :: elem_neighbors  ! (MAX_NV,elem2D)
        integer, allocatable, dimension(:,:) :: nod_in_elem2D   ! (MAX_ADJACENT,nod2D) dense
        integer, allocatable, dimension(:)   :: nod_in_elem2D_num ! (nod2D)

        ! ---- geometry ----
        real(kind=MP), allocatable, dimension(:)   :: elem_area
        ! R7 (FESOM3): edge_dxdy is in PHYSICAL measure, METRES -- FESOM2 stored radians
        ! and applied r_earth*mean(elem_cos) at each point of use. edge_len is the edge
        ! length in METRES, sqrt(edge_dxdy(1)^2 + edge_dxdy(2)^2). edge_cross_dxdy was
        ! already in metres. See compute_edge_geometry.
        real(kind=MP), allocatable, dimension(:,:) :: edge_dxdy         ! (2,edge2D) [m]
        real(kind=MP), allocatable, dimension(:)   :: edge_len          ! (edge2D)   [m]
        real(kind=MP), allocatable, dimension(:,:) :: edge_cross_dxdy   ! (4,edge2D) [m]
        real(kind=MP), allocatable, dimension(:)   :: elem_cos
        real(kind=MP), allocatable, dimension(:)   :: metric_factor
        real(kind=MP), allocatable, dimension(:,:) :: x_corners, y_corners
        real(kind=MP), allocatable, dimension(:)   :: depth
        real(kind=MP), allocatable, dimension(:,:) :: gradient_vec  ! velocity reconstruction
        real(kind=MP), allocatable, dimension(:,:) :: gradient_sca  ! (2*MAX_NV,elem2D) scalar gradient
        ! Kite (median-dual) weights: elem_w(k,e) = the share of the area of element e that
        ! belongs to its vertex k (sum 1; triangle 1/3,1/3,1/3,0), so an element mean is
        ! sum(elem_w(:,e)*f(elem2D_nodes(:,e))). nod_in_elem2D_w(k,n) = the weight of node n
        ! in its k-th element nod_in_elem2D(k,n), for node gathers
        ! sum_k nod_in_elem2D_w(k,n)*elem_area(e_k)*x(e_k) / area(n).
        real(kind=MP), allocatable, dimension(:,:) :: elem_w           ! (MAX_NV,owned elems)
        real(kind=MP), allocatable, dimension(:,:) :: nod_in_elem2D_w  ! (MAX_ADJACENT,nod2D)
        integer,       allocatable, dimension(:)   :: bc_index_nod2D

        ! ---- vertical structure ----
        !
        ! FESOM3 BOTTOM AT VERTICES. The VERTEX column is AUTHORITATIVE; the element
        ! vertical bounds are DERIVED from it and are never read from a mesh file.
        ! (FESOM2 was the other way round: elvls.out defined the element bottom and the
        ! node bottom was the MAX over adjacent elements. elvls.out is no longer read.)
        !
        ! Design note "Bottom implementation for FESOM3" (10 Sep 2026) name mapping --
        ! the note's identifiers are conceptual; these are the real ones:
        !     tlayer(v)      == ulevels_nod2D(v)        blayer(v)      == nlevels_nod2D(v)-1
        !     tlayer_elem(e) == ulevels(e)              blayer_elem(e) == nlevels(e)-1
        ! The note's `ulayer_edge` is min(blayer(ednodes)) -- a BOTTOM bound despite the
        ! `u`. It is never stored, so the naming inconsistency does not reach code.
        !
        ! LEVEL indexing (not layer indexing): nlevels* are level COUNTS, and a column's
        ! layers run from its upper level to its bottom level MINUS ONE:
        !     layers of vertex v :  nz = ulevels_nod2D(v) .. nlevels_nod2D(v)-1
        !     layers of element e:  nz = ulevels(e)       .. nlevels(e)-1
        ! A scalar cell (nz,v) exists iff nz is in the vertex range; hnode(nz,v) is its
        ! thickness and area(nz,v) its horizontal area (depth-independent, see below).
        !
        ! DERIVED in mod_mesh_read (setup_vertical / read_mesh_local):
        !     ulevels(e) = maxval(ulevels_nod2D(elem2D_nodes(1:nnodes,e)))
        !     nlevels(e) = minval(nlevels_nod2D(elem2D_nodes(1:nnodes,e)))
        ! so an element's layer range is exactly its FULLY WET prisms. Velocity DOF
        ! outside that range are never assembled or updated, which is what makes
        ! "velocities touching topography are zero" hold by construction rather than by
        ! scattered if-statements. Bottom drag therefore lands at nlevels(elem)-1 and the
        ! stiffness integration spans zbar(ulevels(e))..zbar_e_bot(e) with no extra code.
        !
        ! RETAINED, and NOT aliases of the vertex column:
        !     nlevels_nod2D_min(n) = min over e in adj(n) of nlevels(e)   ! 2-ring min
        !     ulevels_nod2D_max(n) = max over e in adj(n) of ulevels(e)
        ! These bound work that reaches ADJACENT ELEMENTS from a node, not the node's own
        ! cell, and they are several levels away from the vertex column over most of a
        ! real mesh (pi: mean -3.4, min -28). Do NOT substitute nlevels_nod2D for them.
        ! (The MUSCL read of tr_xy at the up/downwind elements, adv_tra_hor_muscl in
        ! oce_adv_tra_hor.F90, tests the up/downwind element's own wet range: tr_xy is
        ! uninitialized below an element's bottom.)
        !
        ! REQUIRED INVARIANT, asserted in setup_vertical / setup_vertical_local:
        !     maxval over e in adj(n) of nlevels(e) == nlevels_nod2D(n)
        ! i.e. every vertex's deepest scalar cell has at least one wet adjacent element.
        ! It holds because nlvls.out is exactly the max over adjacent elvls.out. This is
        ! not a quality metric: three UNGUARDED divides depend on it and a mesh that
        ! violates it yields NaN at oce_ale.F90:88 (tx/tvol), oce_ale.F90:377
        ! (Wvel/area) and oce_pressure_bv.F90:310 (1/(3*vol)).
        !
        integer :: nl = 0
        real(kind=MP), allocatable, dimension(:) :: zbar, Z, elem_depth
        ! ulevels/nlevels are DERIVED (see above), never read from file.
        integer, allocatable, dimension(:) :: ulevels, ulevels_nod2D, ulevels_nod2D_max
        integer, allocatable, dimension(:) :: nlevels, nlevels_nod2D, nlevels_nod2D_min

        ! ---- control-volume areas ----
        !
        ! FESOM3 BOTTOM AT VERTICES: these are 1-D (per vertex), NOT per level. The
        ! scalar cell is a straight prism -- one bottom level, and the full median-dual
        ! area at every wet layer -- so the horizontal area of a scalar column does not
        ! vary with depth. This is the design note's `area(1:myDim+eDim)`.
        !
        ! Keeping them 2-D would be worse than wasteful: `area(nz,n)` still READS as
        ! level-dependent at every call site, which is exactly how a wrong denominator
        ! hides. An adjacent element that is dry at some level is not missing data to be
        ! dropped from a node average -- its velocity/flux there is genuinely ZERO and it
        ! must contribute that zero carrying its FULL area weight. 1-D makes that the only
        ! expressible choice.
        real(kind=MP), allocatable, dimension(:) :: area, area_inv
        real(kind=MP), allocatable, dimension(:) :: areasvol, areasvol_inv
        real(kind=MP), allocatable, dimension(:)   :: mesh_resolution

        ! ---- elevation stiffness matrix ----
        type(sparse_matrix) :: ssh_stiff

        ! ---- node neighbourhood (CSR-like dense) ----
        integer :: nn_size = 0
        integer, allocatable, dimension(:)   :: nn_num
        integer, allocatable, dimension(:,:) :: nn_pos

        ! ---- ALE layer thickness / depths ----
        real(kind=MP), allocatable, dimension(:,:) :: hnode, hnode_new
        real(kind=MP), allocatable, dimension(:,:) :: zbar_3d_n, Z_3d_n
        real(kind=MP), allocatable, dimension(:,:) :: helem
        real(kind=MP), allocatable, dimension(:)   :: zbar_e_bot     ! (elem2D) depth of partial-cell bottom (full cells: zbar(nlevels(elem)))
        real(kind=MP), allocatable, dimension(:)   :: dhe, hbar, hbar_old

        ! ---- Coriolis ----
        real(kind=MP), allocatable, dimension(:) :: coriolis       ! at elements
        real(kind=MP), allocatable, dimension(:) :: coriolis_node  ! at nodes

#if !defined(USE_HALF_PRECISION)
    contains
        procedure :: write_unformatted => write_t_mesh
        procedure :: read_unformatted  => read_t_mesh
        generic   :: write(unformatted) => write_unformatted
        generic   :: read(unformatted)  => read_unformatted
#endif
    end type t_mesh

#if !defined(USE_HALF_PRECISION)
contains

    subroutine write_t_mesh(mesh, unit, iostat, iomsg)
        class(t_mesh), intent(in)    :: mesh
        integer,       intent(in)    :: unit
        integer,       intent(out)   :: iostat
        character(*),  intent(inout) :: iomsg
        write(unit, iostat=iostat, iomsg=iomsg) mesh%nod2D, mesh%elem2D, &
            mesh%edge2D, mesh%edge2D_in
        write(unit, iostat=iostat, iomsg=iomsg) mesh%ocean_area, mesh%ocean_areawithcav
        call write_bin_array(mesh%coord_nod2D,        unit, iostat, iomsg)
        call write_bin_array(mesh%geo_coord_nod2D,    unit, iostat, iomsg)
        call write_bin_array(mesh%elem2D_nodes,       unit, iostat, iomsg)
        call write_bin_array(mesh%elem2D_nnodes,      unit, iostat, iomsg)
        call write_bin_array(mesh%edges,              unit, iostat, iomsg)
        call write_bin_array(mesh%edge_tri,           unit, iostat, iomsg)
        call write_bin_array(mesh%elem_edges,         unit, iostat, iomsg)
        call write_bin_array(mesh%elem_neighbors,     unit, iostat, iomsg)
        call write_bin_array(mesh%nod_in_elem2D,      unit, iostat, iomsg)
        call write_bin_array(mesh%nod_in_elem2D_num,  unit, iostat, iomsg)
        call write_bin_array(mesh%elem_area,          unit, iostat, iomsg)
        call write_bin_array(mesh%edge_dxdy,          unit, iostat, iomsg)
        call write_bin_array(mesh%edge_len,           unit, iostat, iomsg)
        call write_bin_array(mesh%edge_cross_dxdy,    unit, iostat, iomsg)
        call write_bin_array(mesh%elem_cos,           unit, iostat, iomsg)
        call write_bin_array(mesh%metric_factor,      unit, iostat, iomsg)
        call write_bin_array(mesh%x_corners,          unit, iostat, iomsg)
        call write_bin_array(mesh%y_corners,          unit, iostat, iomsg)
        call write_bin_array(mesh%depth,              unit, iostat, iomsg)
        call write_bin_array(mesh%gradient_vec,       unit, iostat, iomsg)
        call write_bin_array(mesh%gradient_sca,       unit, iostat, iomsg)
        call write_bin_array(mesh%elem_w,             unit, iostat, iomsg)
        call write_bin_array(mesh%nod_in_elem2D_w,    unit, iostat, iomsg)
        call write_bin_array(mesh%bc_index_nod2D,     unit, iostat, iomsg)
        write(unit, iostat=iostat, iomsg=iomsg) mesh%nl
        call write_bin_array(mesh%zbar,               unit, iostat, iomsg)
        call write_bin_array(mesh%Z,                  unit, iostat, iomsg)
        call write_bin_array(mesh%elem_depth,         unit, iostat, iomsg)
        call write_bin_array(mesh%ulevels,            unit, iostat, iomsg)
        call write_bin_array(mesh%ulevels_nod2D,      unit, iostat, iomsg)
        call write_bin_array(mesh%ulevels_nod2D_max,  unit, iostat, iomsg)
        call write_bin_array(mesh%nlevels,            unit, iostat, iomsg)
        call write_bin_array(mesh%nlevels_nod2D,      unit, iostat, iomsg)
        call write_bin_array(mesh%nlevels_nod2D_min,  unit, iostat, iomsg)
        call write_bin_array(mesh%area,               unit, iostat, iomsg)
        call write_bin_array(mesh%area_inv,           unit, iostat, iomsg)
        call write_bin_array(mesh%areasvol,           unit, iostat, iomsg)
        call write_bin_array(mesh%areasvol_inv,       unit, iostat, iomsg)
        call write_bin_array(mesh%mesh_resolution,    unit, iostat, iomsg)
        write(unit, iostat=iostat, iomsg=iomsg) mesh%ssh_stiff%dim, mesh%ssh_stiff%nza
        call write_bin_array(mesh%ssh_stiff%rowptr,     unit, iostat, iomsg)
        call write_bin_array(mesh%ssh_stiff%colind,     unit, iostat, iomsg)
        call write_bin_array(mesh%ssh_stiff%values,     unit, iostat, iomsg)
        call write_bin_array(mesh%ssh_stiff%colind_loc, unit, iostat, iomsg)
        call write_bin_array(mesh%ssh_stiff%rowptr_loc, unit, iostat, iomsg)
        write(unit, iostat=iostat, iomsg=iomsg) mesh%nn_size
        call write_bin_array(mesh%nn_num,             unit, iostat, iomsg)
        call write_bin_array(mesh%nn_pos,             unit, iostat, iomsg)
        call write_bin_array(mesh%hnode,              unit, iostat, iomsg)
        call write_bin_array(mesh%hnode_new,          unit, iostat, iomsg)
        call write_bin_array(mesh%zbar_3d_n,          unit, iostat, iomsg)
        call write_bin_array(mesh%Z_3d_n,             unit, iostat, iomsg)
        call write_bin_array(mesh%helem,              unit, iostat, iomsg)
        call write_bin_array(mesh%zbar_e_bot,         unit, iostat, iomsg)
        call write_bin_array(mesh%dhe,                unit, iostat, iomsg)
        call write_bin_array(mesh%hbar,               unit, iostat, iomsg)
        call write_bin_array(mesh%hbar_old,           unit, iostat, iomsg)
        call write_bin_array(mesh%coriolis,           unit, iostat, iomsg)
        call write_bin_array(mesh%coriolis_node,      unit, iostat, iomsg)
    end subroutine write_t_mesh

    subroutine read_t_mesh(mesh, unit, iostat, iomsg)
        class(t_mesh), intent(inout) :: mesh
        integer,       intent(in)    :: unit
        integer,       intent(out)   :: iostat
        character(*),  intent(inout) :: iomsg
        read(unit, iostat=iostat, iomsg=iomsg) mesh%nod2D, mesh%elem2D, &
            mesh%edge2D, mesh%edge2D_in
        read(unit, iostat=iostat, iomsg=iomsg) mesh%ocean_area, mesh%ocean_areawithcav
        call read_bin_array(mesh%coord_nod2D,        unit, iostat, iomsg)
        call read_bin_array(mesh%geo_coord_nod2D,    unit, iostat, iomsg)
        call read_bin_array(mesh%elem2D_nodes,       unit, iostat, iomsg)
        call read_bin_array(mesh%elem2D_nnodes,      unit, iostat, iomsg)
        call read_bin_array(mesh%edges,              unit, iostat, iomsg)
        call read_bin_array(mesh%edge_tri,           unit, iostat, iomsg)
        call read_bin_array(mesh%elem_edges,         unit, iostat, iomsg)
        call read_bin_array(mesh%elem_neighbors,     unit, iostat, iomsg)
        call read_bin_array(mesh%nod_in_elem2D,      unit, iostat, iomsg)
        call read_bin_array(mesh%nod_in_elem2D_num,  unit, iostat, iomsg)
        call read_bin_array(mesh%elem_area,          unit, iostat, iomsg)
        call read_bin_array(mesh%edge_dxdy,          unit, iostat, iomsg)
        call read_bin_array(mesh%edge_len,           unit, iostat, iomsg)
        call read_bin_array(mesh%edge_cross_dxdy,    unit, iostat, iomsg)
        call read_bin_array(mesh%elem_cos,           unit, iostat, iomsg)
        call read_bin_array(mesh%metric_factor,      unit, iostat, iomsg)
        call read_bin_array(mesh%x_corners,          unit, iostat, iomsg)
        call read_bin_array(mesh%y_corners,          unit, iostat, iomsg)
        call read_bin_array(mesh%depth,              unit, iostat, iomsg)
        call read_bin_array(mesh%gradient_vec,       unit, iostat, iomsg)
        call read_bin_array(mesh%gradient_sca,       unit, iostat, iomsg)
        call read_bin_array(mesh%elem_w,             unit, iostat, iomsg)
        call read_bin_array(mesh%nod_in_elem2D_w,    unit, iostat, iomsg)
        call read_bin_array(mesh%bc_index_nod2D,     unit, iostat, iomsg)
        read(unit, iostat=iostat, iomsg=iomsg) mesh%nl
        call read_bin_array(mesh%zbar,               unit, iostat, iomsg)
        call read_bin_array(mesh%Z,                  unit, iostat, iomsg)
        call read_bin_array(mesh%elem_depth,         unit, iostat, iomsg)
        call read_bin_array(mesh%ulevels,            unit, iostat, iomsg)
        call read_bin_array(mesh%ulevels_nod2D,      unit, iostat, iomsg)
        call read_bin_array(mesh%ulevels_nod2D_max,  unit, iostat, iomsg)
        call read_bin_array(mesh%nlevels,            unit, iostat, iomsg)
        call read_bin_array(mesh%nlevels_nod2D,      unit, iostat, iomsg)
        call read_bin_array(mesh%nlevels_nod2D_min,  unit, iostat, iomsg)
        call read_bin_array(mesh%area,               unit, iostat, iomsg)
        call read_bin_array(mesh%area_inv,           unit, iostat, iomsg)
        call read_bin_array(mesh%areasvol,           unit, iostat, iomsg)
        call read_bin_array(mesh%areasvol_inv,       unit, iostat, iomsg)
        call read_bin_array(mesh%mesh_resolution,    unit, iostat, iomsg)
        read(unit, iostat=iostat, iomsg=iomsg) mesh%ssh_stiff%dim, mesh%ssh_stiff%nza
        call read_bin_array(mesh%ssh_stiff%rowptr,     unit, iostat, iomsg)
        call read_bin_array(mesh%ssh_stiff%colind,     unit, iostat, iomsg)
        call read_bin_array(mesh%ssh_stiff%values,     unit, iostat, iomsg)
        call read_bin_array(mesh%ssh_stiff%colind_loc, unit, iostat, iomsg)
        call read_bin_array(mesh%ssh_stiff%rowptr_loc, unit, iostat, iomsg)
        read(unit, iostat=iostat, iomsg=iomsg) mesh%nn_size
        call read_bin_array(mesh%nn_num,             unit, iostat, iomsg)
        call read_bin_array(mesh%nn_pos,             unit, iostat, iomsg)
        call read_bin_array(mesh%hnode,              unit, iostat, iomsg)
        call read_bin_array(mesh%hnode_new,          unit, iostat, iomsg)
        call read_bin_array(mesh%zbar_3d_n,          unit, iostat, iomsg)
        call read_bin_array(mesh%Z_3d_n,             unit, iostat, iomsg)
        call read_bin_array(mesh%helem,              unit, iostat, iomsg)
        call read_bin_array(mesh%zbar_e_bot,         unit, iostat, iomsg)
        call read_bin_array(mesh%dhe,                unit, iostat, iomsg)
        call read_bin_array(mesh%hbar,               unit, iostat, iomsg)
        call read_bin_array(mesh%hbar_old,           unit, iostat, iomsg)
        call read_bin_array(mesh%coriolis,           unit, iostat, iomsg)
        call read_bin_array(mesh%coriolis_node,      unit, iostat, iomsg)
    end subroutine read_t_mesh
#endif

end module mod_mesh
