program test_halo_ext
    ! Hybrid meshes (Task 8): round trips of the edge exchange (com_edge2D) and the double
    ! node halo exchange (com_nod2D_full) read from dist_N/format.out.
    !   nodes: owned = global id, everything else -999 -> exchange_nod_full (2D, 3D) ->
    !          all myDim+eDim+eXDim_nod2D entries == myList_nod2D; owned_bounds nNodX.
    !   edges: an edge is set (= global id) only on its UNIQUE owner, the owner of
    !          edges(1,edge) (edges.out, global node owner gathered from the owned lists)
    !          -> exchange_edge (2D, 3D) -> all myDim+eDim_edge2D entries == myList_edge2D,
    !          incl. the edges in the owned range of two ranks (counted, must exist npes>1).
    !   lists: com_nod2D_full receives only halo slots, com_edge2D only foreign-owned edges.
    ! npes==1 (synthesis): the lists are empty, the exchanges are no-ops.
    ! A dist_N without format.out (the original pi) stops in exchange_edge with a message
    ! (ctest test_halo_ext_noformat_np2: PASS regex on it, FAIL regex on 'did NOT stop').
    use mpi
    use mod_precision,    only: WP
    use mod_mesh,         only: t_mesh
    use mod_partit,       only: t_partit
    use mod_partitioning, only: par_init, par_ex, set_partition, read_mesh_dims
    use mod_part_bounds,  only: owned_bounds
    use mod_halo,         only: exchange_edge, exchange_nod_full
    implicit none

    character(len=512) :: mesh_dir
    type(t_partit) :: partit
    type(t_mesh)   :: mesh
    integer :: ierr, nfail, nfail_g, u, k, e, ge, p
    integer :: nod2D, elem2D, edge2D, edge2D_in
    integer :: nNodO, nNodL, nNodX, nEdgeO, nElemO, nEdgeL, nX
    integer :: noverlap, noverlap_g
    integer, parameter :: NL = 5
    integer,       allocatable :: owner(:), edges(:,:), cnts(:), displs(:), ownlist(:)
    real(kind=WP), allocatable :: a(:), a3(:,:), ae(:), ae3(:,:)

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) mesh_dir = '/work/ab0995/a270029/meshes/nemo025q_cut'

    nfail = 0
    call par_init(partit)
    call set_partition(partit, trim(mesh_dir))
    call read_mesh_dims(trim(mesh_dir), nod2D, elem2D, edge2D, edge2D_in)
    nNodL  = partit%myDim_nod2D + partit%eDim_nod2D
    nX     = nNodL + partit%eXDim_nod2D
    nEdgeL = partit%myDim_edge2D + partit%eDim_edge2D

    ! ---- error case: FESOM2-layout dist_N (no format.out) -> the call stops ----
    if (partit%dist_format == 0) then
        allocate(ae(nEdgeL)); ae = 0.0_WP
        call exchange_edge(ae, partit)
        write(*,'(a)') 'test_halo_ext: exchange_edge without format.out did NOT stop'
        call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
    end if

    ! ---- owned_bounds: nNodX ----
    mesh%nod2D = nod2D; mesh%elem2D = elem2D; mesh%edge2D = edge2D
    call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit, nNodX=nNodX)
    call check(nNodX == nX, 'owned_bounds: nNodX == myDim+eDim+eXDim_nod2D')
    call check(nNodX == size(partit%myList_nod2D), 'owned_bounds: nNodX == size(myList_nod2D)')
    if (partit%npes == 1) call check(nNodX == nod2D, '1-rank: nNodX == nod2D')

    ! ---- receive lists: halo slots only ----
    call check(all(partit%com_nod2D_full%rlist > partit%myDim_nod2D .and. &
                   partit%com_nod2D_full%rlist <= nX), 'com_nod2D_full receives halo slots only')
    call check(all(partit%com_nod2D_full%slist >= 1 .and. &
                   partit%com_nod2D_full%slist <= partit%myDim_nod2D), 'com_nod2D_full sends owned only')
    call check(all(partit%com_edge2D%rlist >= 1 .and. partit%com_edge2D%rlist <= nEdgeL), &
               'com_edge2D receive slots in 1..myDim+eDim_edge2D')

    ! ---- double node halo: 2D ----
    allocate(a(nX))
    a = -999.0_WP
    a(1:partit%myDim_nod2D) = real(partit%myList_nod2D(1:partit%myDim_nod2D), WP)
    call exchange_nod_full(a, partit)
    call check(all(a == real(partit%myList_nod2D(1:nX), WP)), 'node full 2D round trip')

    ! ---- double node halo: 3D (level k = gid + k) ----
    allocate(a3(NL, nX))
    a3 = -999.0_WP
    do e = 1, partit%myDim_nod2D
        do k = 1, NL
            a3(k, e) = real(partit%myList_nod2D(e) + k, WP)
        end do
    end do
    call exchange_nod_full(a3, partit)
    call check(lev_ok(a3, partit%myList_nod2D(1:nX)), 'node full 3D round trip')

    ! ---- global node owner from the owned lists ----
    allocate(owner(nod2D), cnts(partit%npes), displs(partit%npes), ownlist(nod2D))
    call MPI_Allgather(partit%myDim_nod2D, 1, MPI_INTEGER, cnts, 1, MPI_INTEGER, &
                       partit%MPI_COMM_FESOM, ierr)
    displs(1) = 0
    do p = 2, partit%npes
        displs(p) = displs(p-1) + cnts(p-1)
    end do
    call MPI_Allgatherv(partit%myList_nod2D, partit%myDim_nod2D, MPI_INTEGER, ownlist, cnts, &
                        displs, MPI_INTEGER, partit%MPI_COMM_FESOM, ierr)
    owner = -1
    do p = 1, partit%npes
        owner(ownlist(displs(p)+1:displs(p)+cnts(p))) = p - 1
    end do
    call check(all(owner >= 0), 'every node owned by one rank')

    ! ---- global edges (edges.out, the numbering fesom_ini partitions) ----
    allocate(edges(2, edge2D))
    open(newunit=u, file=trim(mesh_dir)//'/edges.out', status='old', action='read')
    read(u,*) edges
    close(u)

    ! ---- edges: set on the unique owner only, exchange, compare ----
    allocate(ae(nEdgeL), ae3(NL, nEdgeL))
    ae = -999.0_WP; ae3 = -999.0_WP
    noverlap = 0
    do e = 1, nEdgeL
        ge = partit%myList_edge2D(e)
        if (owner(edges(1, ge)) == partit%mype) then
            ae(e) = real(ge, WP)
            do k = 1, NL
                ae3(k, e) = real(ge + k, WP)
            end do
        else if (e <= partit%myDim_edge2D) then
            noverlap = noverlap + 1      ! in my owned range, owned uniquely by another rank
        end if
    end do
    call check(all(owner(edges(1, partit%myList_edge2D(partit%com_edge2D%rlist))) /= partit%mype), &
               'com_edge2D receives foreign-owned edges only')
    call exchange_edge(ae, partit)
    call exchange_edge(ae3, partit)
    call check(all(ae == real(partit%myList_edge2D(1:nEdgeL), WP)), 'edge 2D round trip')
    call check(lev_ok(ae3, partit%myList_edge2D(1:nEdgeL)), 'edge 3D round trip')
    call MPI_Allreduce(noverlap, noverlap_g, 1, MPI_INTEGER, MPI_SUM, partit%MPI_COMM_FESOM, ierr)
    if (partit%npes > 1) call check(noverlap_g > 0, 'edges in the owned range of two ranks exist')

    call MPI_Allreduce(nfail, nfail_g, 1, MPI_INTEGER, MPI_SUM, partit%MPI_COMM_FESOM, ierr)
    if (partit%mype == 0) then
        write(*,'(a,i0,a,i0)') 'test_halo_ext: owned-range edges received from the unique owner: ', &
            noverlap_g, ', double-halo nodes on rank 0: ', partit%eXDim_nod2D
        if (nfail_g == 0) then
            write(*,'(a,i0,a)') 'test_halo_ext (npes=', partit%npes, '): ALL PASS'
        else
            write(*,'(a,i0,a,i0,a)') 'test_halo_ext (npes=', partit%npes, '): ', nfail_g, ' FAILURE(S)'
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

    logical function lev_ok(f, gid)
        ! f(k, j) == gid(j) + k on every level
        real(kind=WP), intent(in) :: f(:,:)
        integer,       intent(in) :: gid(:)
        integer :: j, k
        lev_ok = .true.
        do j = 1, size(f, 2)
            do k = 1, size(f, 1)
                if (f(k, j) /= real(gid(j) + k, WP)) lev_ok = .false.
            end do
        end do
    end function
end program test_halo_ext
