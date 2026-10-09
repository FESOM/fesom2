program test_partit
    ! M0.4 gate: 1-rank synthesis + dist_<NP>/ reading on the pi mesh.
    ! Runs at 1, 2, 8 ranks. Invariants: owned nodes/edges partition exactly
    ! (sum == global), owned+halo elements are redundant (sum >= global).
    ! Hybrid meshes: with dist_N/format.out the appended sections are read
    ! (eXDim_nod2D appended to myList_nod2D, com_nod2D_full, com_edge2D); the 1-rank
    ! synthesis provides them empty. The global counts come from the mesh headers, so
    ! FESOM3_MESH_DIR may point to any mesh (default pi). FESOM3_EXPECT_FORMAT=1: the
    ! appended sections are required (dist_N/format.out must be read).
    ! FESOM3_TEST_MANY_NEIGHBOURS=recv|send (np 2): a synthetic dist_2 with 33 neighbour PEs
    ! in com_nod2D_full must stop set_partition (the failure is 'did NOT stop').
    use mpi
    use mod_partit,        only: t_partit, com_struct
    use mod_partitioning,  only: par_init, par_ex, set_partition, read_mesh_dims
    use mod_io_zarr,       only: zarr_mkdir
    implicit none

    character(len=512) :: mesh_dir, cbuf
    type(t_partit) :: partit
    integer :: ierr, nfail, nfail_g
    integer :: sum_nod, sum_elem, sum_edge
    integer :: nod2D, elem2D, edge2D, edge2D_in, nNodL, nNodX, nEdgeL, ns, nr
    integer :: cnt(2), cnt_g(2)

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/home/a/a270088/port2/fesom2/tests/data/MESHES/pi'

    nfail = 0
    call par_init(partit)
    cbuf = ''
    call get_environment_variable('FESOM3_TEST_MANY_NEIGHBOURS', cbuf)
    if (len_trim(cbuf) > 0) call many_neighbours_case(trim(cbuf))
    call set_partition(partit, trim(mesh_dir))
    call read_mesh_dims(trim(mesh_dir), nod2D, elem2D, edge2D, edge2D_in)
    nNodL  = partit%myDim_nod2D + partit%eDim_nod2D
    nNodX  = nNodL + partit%eXDim_nod2D
    nEdgeL = partit%myDim_edge2D + partit%eDim_edge2D

    ! ---- global partition invariants (all rank counts) ----
    call MPI_Allreduce(partit%myDim_nod2D,  sum_nod,  1, MPI_INTEGER, MPI_SUM, partit%MPI_COMM_FESOM, ierr)
    call MPI_Allreduce(partit%myDim_elem2D, sum_elem, 1, MPI_INTEGER, MPI_SUM, partit%MPI_COMM_FESOM, ierr)
    call MPI_Allreduce(partit%myDim_edge2D, sum_edge, 1, MPI_INTEGER, MPI_SUM, partit%MPI_COMM_FESOM, ierr)

    ! Nodes are uniquely owned (sum == nod2D); elements and edges have
    ! boundary-redundant owned counts (sum >= global), verified on pi/dist_2.
    call check(sum_nod  == nod2D,  'sum(myDim_nod2D) == nod2D (unique)')
    call check(sum_elem >= elem2D, 'sum(myDim_elem2D) >= elem2D (boundary-redundant)')
    call check(sum_edge >= edge2D, 'sum(myDim_edge2D) >= edge2D (boundary-redundant)')

    ! ---- local sanity ----
    call check(partit%myDim_nod2D > 0, 'myDim_nod2D > 0')
    call check(size(partit%myList_nod2D) == nNodX, 'myList_nod2D size == myDim+eDim+eXDim')
    call check(all(partit%myList_nod2D >= 1 .and. partit%myList_nod2D <= nod2D), 'node gids in range')
    call check(all(partit%myList_elem2D >= 1 .and. partit%myList_elem2D <= elem2D), 'elem gids in range')

    cbuf = ''
    call get_environment_variable('FESOM3_EXPECT_FORMAT', cbuf)
    if (trim(cbuf) == '1') call check(partit%dist_format >= 1, 'FESOM3_EXPECT_FORMAT: appended sections read')

    ! ---- 1-rank synthesis specifics ----
    if (partit%npes == 1) then
        call check(partit%myDim_nod2D == nod2D, '1-rank: myDim_nod2D == nod2D')
        call check(partit%eDim_nod2D == 0 .and. partit%eXDim_nod2D == 0, '1-rank: eDim/eXDim_nod2D == 0')
        call check(partit%eDim_elem2D == 0 .and. partit%eXDim_elem2D == 0, '1-rank: elem halos == 0')
        call check(partit%myList_nod2D(1) == 1 .and. partit%myList_nod2D(nod2D) == nod2D, &
                   '1-rank: identity node map')
        call check(partit%com_nod2D%rPEnum == 0 .and. partit%com_nod2D%sPEnum == 0, '1-rank: no neighbours')
        call check(partit%dist_format == 1, '1-rank: appended lists provided')
        call check(empty_com(partit%com_nod2D_full) .and. empty_com(partit%com_edge2D), &
                   '1-rank: com_nod2D_full, com_edge2D empty')
    else
        ! multi-rank: every PE should have at least some halo and neighbours
        call check(partit%eDim_nod2D > 0, 'multi-rank: eDim_nod2D > 0')
        call check(partit%com_nod2D%rPEnum > 0, 'multi-rank: has receive neighbours')
        call check(size(partit%com_nod2D%rlist) == partit%eDim_nod2D, 'multi-rank: rlist size == eDim_nod2D')
        if (partit%dist_format >= 1) then
            ! com_nod2D_full receives each halo node (eDim + eXDim) once, sends owned ones
            call check(size(partit%com_nod2D_full%rlist) == partit%eDim_nod2D + partit%eXDim_nod2D, &
                       'format: com_nod2D_full rlist size == eDim+eXDim')
            call check(is_perm(partit%com_nod2D_full%rlist, partit%myDim_nod2D+1, nNodX), &
                       'format: com_nod2D_full receives every halo node once')
            call check(all(partit%com_nod2D_full%slist >= 1 .and. &
                           partit%com_nod2D_full%slist <= partit%myDim_nod2D), &
                       'format: com_nod2D_full sends owned nodes')
            call check(all(partit%com_edge2D%rlist >= 1 .and. partit%com_edge2D%rlist <= nEdgeL) .and. &
                       all(partit%com_edge2D%slist >= 1 .and. partit%com_edge2D%slist <= nEdgeL), &
                       'format: com_edge2D lists in 1..myDim+eDim_edge2D')
            call check(is_perm(partit%com_edge2D%rlist, 0, -1), 'format: com_edge2D receives an edge once')
            ! what is sent is received: global send count == global receive count
            ns = size(partit%com_nod2D_full%slist); nr = size(partit%com_nod2D_full%rlist)
            cnt = [ns - nr, size(partit%com_edge2D%slist) - size(partit%com_edge2D%rlist)]
            call MPI_Allreduce(cnt, cnt_g, 2, MPI_INTEGER, MPI_SUM, partit%MPI_COMM_FESOM, ierr)
            call check(all(cnt_g == 0), 'format: sum of sends == sum of receives (nod_full, edge)')
        else
            call check(partit%eXDim_nod2D == 0, 'no format.out: eXDim_nod2D == 0')
        end if
    end if

    call MPI_Allreduce(nfail, nfail_g, 1, MPI_INTEGER, MPI_SUM, partit%MPI_COMM_FESOM, ierr)
    if (partit%mype == 0) then
        write(*,'(a,i0,a,i0)') 'test_partit: '//trim(mesh_dir)//' npes=', partit%npes, &
            ' dist_format=', partit%dist_format
        if (nfail_g == 0) then
            write(*,'(a,i0,a)') 'test_partit (npes=', partit%npes, '): ALL PASS'
        else
            write(*,'(a,i0,a,i0,a)') 'test_partit (npes=', partit%npes, '): ', nfail_g, ' FAILURE(S)'
        end if
    end if

    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    if (nfail_g /= 0) error stop 1

contains
    subroutine many_neighbours_case(which)
        ! FESOM3_TEST_MANY_NEIGHBOURS=recv|send: a synthetic dist_2 (one node, element and
        ! edge per PE, written into test_partit_many_neighbours_<which>/ of the run directory:
        ! one directory per variant, so the two ctests can run concurrently) whose
        ! com_nod2D_full has 33 receive (or send) PEs; read_com must stop at the count, before
        ! reading past the com_struct arrays rPE/sPE(MAX_NEIGHBOR_PARTITIONS).
        character(len=*), intent(in) :: which
        character(len=:), allocatable :: d
        ! an empty com_struct: rPEnum, rPE, rptr, rlist, sPEnum, sPE, sptr, slist
        character(len=1), parameter :: empty(8) = ['0', ' ', '1', ' ', '0', ' ', '1', ' ']
        character(len=8) :: cpe
        integer :: u, p
        if (partit%npes /= 2) then
            write(*,'(a)') 'test_partit: the many-neighbours case needs 2 ranks'
            error stop 1
        end if
        if (which /= 'recv' .and. which /= 'send') then
            write(*,'(a)') 'test_partit: unknown FESOM3_TEST_MANY_NEIGHBOURS='//which
            error stop 1
        end if
        d = 'test_partit_many_neighbours_'//which
        if (partit%mype == 0) then
            call zarr_mkdir(d//'/dist_2')   ! POSIX mkdir: no fork after MPI_Init
            open(newunit=u, file=d//'/dist_2/rpart.out', status='replace')
            write(u,'(a)') '2', '1 1'
            close(u)
            open(newunit=u, file=d//'/dist_2/format.out', status='replace')
            write(u,'(a)') '1'
            close(u)
            do p = 0, 1
                write(cpe,'(i5.5)') p
                ! my_list: header, nodes (myDim, eDim, list), elements (myDim, eDim, eXDim,
                ! list), edges (myDim, eDim, list), eXDim_nod2D and its (empty) list
                open(newunit=u, file=d//'/dist_2/my_list'//trim(cpe)//'.out', status='replace')
                write(u,'(i0)') p, 1, 0, p+1, 1, 0, 0, p+1, 1, 0, p+1, 0
                write(u,'(a)') ' '
                close(u)
                ! com_info: header, com_nod2D, com_elem2D, com_elem2D_full (empty), then
                ! com_nod2D_full with 33 receive PEs (or none received and 33 send PEs)
                open(newunit=u, file=d//'/dist_2/com_info'//trim(cpe)//'.out', status='replace')
                write(u,'(i0)') p
                write(u,'(a)') empty, empty, empty
                if (which == 'send') write(u,'(a)') empty(1:4)
                write(u,'(a)') '33', repeat(' 1', 33)
                close(u)
            end do
        end if
        call MPI_Barrier(partit%MPI_COMM_FESOM, ierr)
        call set_partition(partit, d)
        write(*,'(a)') 'test_partit: set_partition did NOT stop on 33 neighbour PEs'
        call par_ex(partit%MPI_COMM_FESOM, partit%mype)
        error stop 1
    end subroutine

    subroutine check(cond, name)
        logical,          intent(in) :: cond
        character(len=*), intent(in) :: name
        if (.not. cond) then
            nfail = nfail + 1
            write(*,'(a,i0,a)') '  FAIL [rank ', partit%mype, ']: '//name
        end if
    end subroutine

    logical function empty_com(com)
        type(com_struct), intent(in) :: com
        empty_com = com%rPEnum == 0 .and. com%sPEnum == 0 .and. size(com%rlist) == 0 &
                    .and. size(com%slist) == 0
    end function

    logical function is_perm(list, lo, hi)
        ! the entries of list are distinct and, if lo <= hi, exactly lo..hi
        integer, intent(in) :: list(:), lo, hi
        logical, allocatable :: seen(:)
        integer :: i
        is_perm = .true.
        if (size(list) == 0) return
        allocate(seen(minval(list):maxval(list)))
        seen = .false.
        do i = 1, size(list)
            if (seen(list(i))) is_perm = .false.
            seen(list(i)) = .true.
        end do
        if (lo <= hi) is_perm = is_perm .and. size(list) == hi-lo+1 .and. &
                                minval(list) == lo .and. maxval(list) == hi
    end function
end program test_partit
