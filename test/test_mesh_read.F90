program test_mesh_read
    ! Hybrid-mesh reading (plan 2026-10-08-hybrid-meshes, Task 6): read_mesh on a mesh
    ! with triangles and/or quads (FESOM3_MESH_DIR, default the NEMO-mixed cut-out) at
    ! 1 or N ranks, checked against the mesh files read independently here:
    !   - counts; every element (owned on some rank) has the vertices and the vertex
    !     count of its elem2d.out line (3 columns or 4, a triangle repeats its 1st
    !     vertex) -> the triangle/quad histogram; slot 4 == slot 1 <=> 3 vertices;
    !   - node levels == nlvls.out for owned, halo and double-halo nodes; element levels
    !     (owned + full halo) == min over the vertices; ulevels == 1;
    !   - zbar == -(levels of depth.out) if there is no aux3d.out;
    !   - dist_N/format.out: the double-halo nodes (eXDim_nod2D) are vertices of local
    !     elements and no owned/halo node, and carry their coordinates and levels;
    !   - no clockwise swap on a mesh with quads (fesom_ini oriented every element);
    !   - enforce_cw_orientation on anticlockwise and clockwise triangles and quads (an
    !     in-memory mesh): turned ones are clockwise with the cycle reversed at vertex 1.
    ! FESOM3_TEST_MALFORMED=k (1 rank): writes a small mesh with a defective elem2d.out
    ! and reads it: 1 = a line with a wrong number of node ids, 2 = a 1st line with
    ! neither 3 nor 4 ids, 3 = a non-integer id, 4 = an id outside 1..nod2D, 5 = a missing
    ! line. read_mesh must stop with a message (ctest PASS_REGULAR_EXPRESSION); reaching
    ! the end is the failure ('did NOT stop').
    use mpi
    use mod_precision,    only: WP, MP
    use mod_constants,    only: rad
    use mod_mesh,         only: t_mesh, MAX_NV
    use mod_partit,       only: t_partit
    use mod_partitioning, only: par_init, par_ex, set_partition, read_mesh_dims
    use mod_mesh_read,    only: read_mesh, enforce_cw_orientation
    implicit none

    character(len=512) :: mesh_dir, cbuf
    type(t_partit) :: partit
    type(t_mesh)   :: mesh
    integer :: ierr, nfail, nfail_g, nsw, nod2D, elem2D, edge2D, edge2D_in
    integer :: nNodO, nNodL, nNodX, nElemO, nElemF, i, g, ncol, ntri, nquad
    integer, allocatable :: enodes(:,:), ennod(:), nlvls(:), seen(:), seen_g(:)
    real(kind=WP), allocatable :: lon(:), lat(:), zlev(:)
    logical :: has_aux3d, any_quad

    call par_init(partit)
    nfail = 0
    cbuf = ''
    call get_environment_variable('FESOM3_TEST_MALFORMED', cbuf)
    if (len_trim(cbuf) > 0) call malformed_case(trim(cbuf))
    call orientation_case()

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) mesh_dir = '/work/ab0995/a270029/meshes/nemo025mix_cut'
    call set_partition(partit, trim(mesh_dir))
    call read_mesh(mesh, partit, trim(mesh_dir), 0.0_WP, 0.0_WP, 0.0_WP, 360.0_WP, &
                   force_rotation=.false., n_cw_swaps=nsw)
    nNodO  = partit%myDim_nod2D
    nNodL  = nNodO + partit%eDim_nod2D
    nNodX  = nNodL + partit%eXDim_nod2D
    nElemO = partit%myDim_elem2D
    nElemF = nElemO + partit%eDim_elem2D + partit%eXDim_elem2D

    ! ---- the mesh files, read independently ----
    call read_mesh_dims(trim(mesh_dir), nod2D, elem2D, edge2D, edge2D_in)
    call read_files()

    call check(mesh%nod2D == nod2D .and. mesh%elem2D == elem2D, 'global counts')
    ntri = count(ennod == 3); nquad = count(ennod == 4)
    any_quad = nquad > 0

    ! ---- elements: vertices and vertex count; histogram over the owners ----
    allocate(seen(elem2D), seen_g(elem2D))
    seen = 0
    do i = 1, nElemO
        g = partit%myList_elem2D(i)
        seen(g) = mesh%elem2D_nnodes(i)
        call check_elem(i, g)
    end do
    call MPI_Allreduce(seen, seen_g, elem2D, MPI_INTEGER, MPI_MAX, partit%MPI_COMM_FESOM, ierr)
    call check(all(seen_g == ennod), 'every element owned somewhere with its vertex count')
    if (partit%mype == 0) write(*,'(a,i0,a,i0,a,i0,a,i0)') 'test_mesh_read: npes=', partit%npes, &
        ' elements=', elem2D, ' triangles=', count(seen_g == 3), ' quads=', count(seen_g == 4)
    call check(nsw == 0 .or. .not. any_quad, 'no clockwise swap on a mesh with quads')

    ! ---- levels ----
    call check(all(mesh%nlevels_nod2D(1:nNodX) == nlvls(partit%myList_nod2D(1:nNodX))), &
               'nlevels_nod2D == nlvls.out (owned, halo, double halo)')
    call check(size(mesh%nlevels_nod2D) == nNodX .and. size(mesh%coord_nod2D, 2) == nNodX, &
               'node file arrays sized myDim+eDim+eXDim')
    do i = 1, nElemF
        g = partit%myList_elem2D(i)
        if (mesh%nlevels(i) /= minval(nlvls(enodes(1:ennod(g), g))) .or. mesh%ulevels(i) /= 1) then
            call check(.false., 'element levels == min over the vertices')
            exit
        end if
    end do
    if (.not. has_aux3d) call check(all(mesh%zbar == -zlev) .and. mesh%nl == size(zlev), &
                                    'zbar == -(depth.out levels)')

    ! ---- double-halo nodes ----
    call check(size(partit%myList_nod2D) == nNodX, 'myList_nod2D holds eXDim_nod2D')
    call check(xnodes_ok(), 'eXDim nodes: vertices of local elements, not owned/halo')
    call check(all(mesh%coord_nod2D(1, 1:nNodX) == real(lon(partit%myList_nod2D(1:nNodX)) * rad, MP)) &
         .and. all(mesh%coord_nod2D(2, 1:nNodX) == real(lat(partit%myList_nod2D(1:nNodX)) * rad, MP)), &
               'coordinates of owned, halo and double-halo nodes')
    if (partit%mype == 0) write(*,'(a,i0,a,i0,a,i0,a,l1)') 'test_mesh_read: rank 0 myDim_nod2D=', &
        nNodO, ' eDim=', partit%eDim_nod2D, ' eXDim=', partit%eXDim_nod2D, ' aux3d=', has_aux3d

    call MPI_Allreduce(nfail, nfail_g, 1, MPI_INTEGER, MPI_SUM, partit%MPI_COMM_FESOM, ierr)
    if (partit%mype == 0) then
        if (nfail_g == 0) then
            write(*,'(a,i0,a)') 'test_mesh_read (npes=', partit%npes, '): ALL PASS'
        else
            write(*,'(a,i0,a,i0,a)') 'test_mesh_read (npes=', partit%npes, '): ', nfail_g, ' FAILURE(S)'
        end if
    end if
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    if (nfail_g /= 0) error stop 1

contains

    subroutine check(cond, name)
        logical,          intent(in) :: cond
        character(len=*), intent(in) :: name
        if (.not. cond) then
            nfail = nfail + 1
            write(*,'(a,i0,a)') '  FAIL [rank ', partit%mype, ']: '//name
        end if
    end subroutine

    subroutine read_files()
        ! nod2d.out (lon/lat), elem2d.out (3/4 columns), nlvls.out, depth.out levels
        integer :: u, ios, id, bnd, j, ntok
        character(len=256) :: line
        allocate(lon(nod2D), lat(nod2D), enodes(MAX_NV, elem2D), ennod(elem2D), nlvls(nod2D))
        open(newunit=u, file=trim(mesh_dir)//'/nod2d.out', status='old', action='read')
        read(u,*)
        do j = 1, nod2D
            read(u,*) id, lon(j), lat(j), bnd
        end do
        close(u)
        open(newunit=u, file=trim(mesh_dir)//'/elem2d.out', status='old', action='read')
        read(u,*)
        read(u,'(a)') line
        ntok = 0
        do j = 1, len_trim(line)
            if (line(j:j) /= ' ' .and. (j == 1 .or. line(max(j-1,1):max(j-1,1)) == ' ')) ntok = ntok + 1
        end do
        ncol = ntok
        backspace(u)
        do j = 1, elem2D
            read(u,*) enodes(1:ncol, j)
        end do
        close(u)
        if (ncol == 3) enodes(4, :) = enodes(1, :)
        ennod = merge(3, 4, enodes(4, :) == enodes(1, :))
        open(newunit=u, file=trim(mesh_dir)//'/nlvls.out', status='old', action='read')
        read(u,*) nlvls
        close(u)
        inquire(file=trim(mesh_dir)//'/aux3d.out', exist=has_aux3d)
        if (.not. has_aux3d) then
            open(newunit=u, file=trim(mesh_dir)//'/depth.out', status='old', action='read', iostat=ios)
            read(u,*) j
            allocate(zlev(j))
            read(u,*) zlev
            close(u)
        end if
    end subroutine

    subroutine check_elem(i, g)
        ! owned element i (global g): the same vertices as the file (as a set: an element
        ! may have been turned), slot 4 == slot 1 exactly for the 3-vertex elements
        integer, intent(in) :: i, g
        integer :: v(MAX_NV), k
        logical :: same
        v = partit%myList_nod2D(mesh%elem2D_nodes(:, i))
        same = mesh%elem2D_nnodes(i) == ennod(g)
        do k = 1, MAX_NV
            same = same .and. any(v == enodes(k, g)) .and. any(enodes(:, g) == v(k))
        end do
        same = same .and. ((v(4) == v(1)) .eqv. (ennod(g) == 3))
        if (.not. same) then
            write(*,'(a,i0,a,4i8,a,4i8)') '  element ', g, ': read', v, ' file', enodes(:, g)
            call check(.false., 'element vertices and vertex count as in elem2d.out')
        end if
    end subroutine

    logical function xnodes_ok()
        ! every eXDim node is a vertex of some local element and not among 1..nNodL;
        ! every vertex of a local element is local (1..nNodX) when format.out was read
        logical, allocatable :: isloc(:), isvert(:)
        integer :: j, e, k
        allocate(isloc(nod2D), isvert(nod2D))
        isloc = .false.; isvert = .false.
        isloc(partit%myList_nod2D(1:nNodL)) = .true.
        do j = 1, nElemF
            e = partit%myList_elem2D(j)
            do k = 1, ennod(e)
                isvert(enodes(k, e)) = .true.
            end do
        end do
        xnodes_ok = .true.
        do j = nNodL+1, nNodX
            g = partit%myList_nod2D(j)
            if (isloc(g) .or. .not. isvert(g)) xnodes_ok = .false.
        end do
        if (partit%dist_format >= 1) then
            isloc(partit%myList_nod2D(1:nNodX)) = .true.
            if (any(isvert .and. .not. isloc)) xnodes_ok = .false.
        end if
    end function

    subroutine orientation_case()
        ! 7 nodes, 4 elements: an anticlockwise quad and triangle (turned) and their
        ! clockwise copies (kept); after enforce_cw_orientation every element has the
        ! signed area -1 or -1/2 (clockwise, a valid cycle) and slot 1 unchanged
        type(t_mesh) :: tm
        integer :: e, k, ns
        real(kind=WP) :: a, x(MAX_NV), y(MAX_NV)
        real(kind=WP), parameter :: area_ref(4) = [-1.0_WP, -0.5_WP, -1.0_WP, -0.5_WP]
        tm%nod2D = 7; tm%elem2D = 4
        allocate(tm%coord_nod2D(2, 7), tm%elem2D_nodes(MAX_NV, 4), tm%elem2D_nnodes(4))
        tm%coord_nod2D(1, :) = [0.0_WP, 1.0_WP, 1.0_WP, 0.0_WP, 2.0_WP, 2.0_WP, 1.5_WP]
        tm%coord_nod2D(2, :) = [0.0_WP, 0.0_WP, 1.0_WP, 1.0_WP, 0.0_WP, 1.0_WP, 2.0_WP]
        tm%elem2D_nodes(:, 1) = [1, 2, 3, 4]      ! quad, anticlockwise
        tm%elem2D_nodes(:, 2) = [2, 5, 6, 2]      ! triangle, anticlockwise (4th = 1st)
        tm%elem2D_nodes(:, 3) = [1, 4, 3, 2]      ! quad, clockwise
        tm%elem2D_nodes(:, 4) = [2, 6, 5, 2]      ! triangle, clockwise
        tm%elem2D_nnodes = [4, 3, 4, 3]
        call enforce_cw_orientation(tm, ns)
        call check(ns == 2, 'enforce_cw_orientation: 2 swaps (the anticlockwise quad and triangle)')
        call check(all(tm%elem2D_nodes(:, 1) == [1, 4, 3, 2]) .and. &
                   all(tm%elem2D_nodes(:, 2) == [2, 6, 5, 2]) .and. &
                   all(tm%elem2D_nodes(:, 3) == [1, 4, 3, 2]) .and. &
                   all(tm%elem2D_nodes(:, 4) == [2, 6, 5, 2]), 'enforce_cw_orientation: cycles')
        do e = 1, 4
            x = tm%coord_nod2D(1, tm%elem2D_nodes(:, e)); y = tm%coord_nod2D(2, tm%elem2D_nodes(:, e))
            a = 0.0_WP
            do k = 1, tm%elem2D_nnodes(e)
                a = a + 0.5_WP * (x(k) * y(mod(k, tm%elem2D_nnodes(e)) + 1) &
                                - x(mod(k, tm%elem2D_nnodes(e)) + 1) * y(k))
            end do
            if (a /= area_ref(e)) call check(.false., 'enforce_cw_orientation: clockwise signed area')
        end do
    end subroutine

    subroutine malformed_case(kind)
        ! 4 nodes, 3 elements with one defect (FESOM3_TEST_MALFORMED, see the header)
        character(len=*), intent(in) :: kind
        character(len=:), allocatable :: d
        integer :: u
        if (partit%npes /= 1) then
            write(*,'(a)') 'test_mesh_read: the malformed case needs 1 rank'
            error stop 1
        end if
        d = 'test_mesh_read_malformed'//kind
        call execute_command_line('mkdir -p '//d)
        open(newunit=u, file=d//'/nod2d.out', status='replace')
        write(u,'(a)') '4', '1 0.0 0.0 1', '2 1.0 0.0 1', '3 1.0 1.0 1', '4 0.0 1.0 1'
        close(u)
        open(newunit=u, file=d//'/elem2d.out', status='replace')
        select case (kind)
        case ('1'); write(u,'(a)') '3', '1 2 3', '1 3 4', '2 4'
        case ('2'); write(u,'(a)') '3', '1 2', '1 3 4', '2 3 4'
        case ('3'); write(u,'(a)') '3', '1 2 3', '1 3 x', '2 3 4'
        case ('4'); write(u,'(a)') '3', '1 2 3', '1 3 9', '2 3 4'
        case ('5'); write(u,'(a)') '3', '1 2 3', '1 3 4'
        case default
            write(*,'(a)') 'test_mesh_read: unknown FESOM3_TEST_MALFORMED='//kind
            error stop 1
        end select
        close(u)
        call read_mesh(mesh, partit, d, 0.0_WP, 0.0_WP, 0.0_WP, 360.0_WP, &
                       force_rotation=.false., n_cw_swaps=nsw)
        write(*,'(a)') 'test_mesh_read: read_mesh did NOT stop on a malformed elem2d.out'
        call par_ex(partit%MPI_COMM_FESOM, partit%mype)
        error stop 1
    end subroutine

end program test_mesh_read
