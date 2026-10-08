program test_halo
    ! M0.5 gate: halo-identity (set field=global id, exchange, assert halo==owner id)
    ! on pi/dist_2 and dist_8, for node 2D/3D/int and element exchanges; plus the
    ! stale-halo probe detects an injected corruption.
    use mpi
    use mod_precision,    only: WP
    use mod_partit,       only: t_partit
    use mod_partitioning, only: par_init, par_ex, set_partition
    use mod_halo,         only: exchange_nod, exchange_elem, stale_halo_max_nod
    implicit none

    character(len=512) :: mesh_dir
    type(t_partit) :: partit
    integer :: ierr, nfail, nfail_g, nn, ne, k
    integer, parameter :: NL = 5
    real(kind=WP), allocatable :: a(:), a3(:,:), ae(:)
    integer,       allocatable :: ai(:)
    real(kind=WP) :: d

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/home/a/a270088/port2/fesom2/tests/data/MESHES/pi'

    nfail = 0
    call par_init(partit)
    call set_partition(partit, trim(mesh_dir))
    nn = partit%myDim_nod2D + partit%eDim_nod2D   ! (myList_nod2D may hold eXDim_nod2D more)
    ne =partit%myDim_elem2D + partit%eDim_elem2D + partit%eXDim_elem2D

    ! ---- node 2D real halo-identity ----
    allocate(a(nn))
    a = -999.0_WP
    a(1:partit%myDim_nod2D) = real(partit%myList_nod2D(1:partit%myDim_nod2D), WP)
    call exchange_nod(a, partit)
    call check(maxval(abs(a - real(partit%myList_nod2D(1:nn), WP))) == 0.0_WP, 'node 2D halo-identity')

    ! ---- stale-halo probe ----
    d = stale_halo_max_nod(a, partit)
    call check(d == 0.0_WP, 'probe: fresh halo -> 0')
    if (partit%eDim_nod2D > 0) then
        a(partit%myDim_nod2D + 1) = a(partit%myDim_nod2D + 1) + 7.0_WP   ! inject stale halo
        d = stale_halo_max_nod(a, partit)
        call check(d > 0.0_WP, 'probe: stale halo detected')
        call exchange_nod(a, partit)                                     ! repair
        d = stale_halo_max_nod(a, partit)
        call check(d == 0.0_WP, 'probe: repaired halo -> 0')
    end if

    ! ---- node 2D integer halo-identity ----
    allocate(ai(nn))
    ai = -999
    ai(1:partit%myDim_nod2D) = partit%myList_nod2D(1:partit%myDim_nod2D)
    call exchange_nod(ai, partit)
    call check(all(ai == partit%myList_nod2D(1:nn)), 'node 2D int halo-identity')

    ! ---- node 3D real halo-identity (each level = gid) ----
    allocate(a3(NL, nn))
    a3 = -999.0_WP
    do k = 1, partit%myDim_nod2D
        a3(:, k) = real(partit%myList_nod2D(k), WP)
    end do
    call exchange_nod(a3, partit)
    call check(node3d_ok(a3, partit), 'node 3D halo-identity')

    ! ---- element 2D real halo-identity (com_elem2D fills myDim+eDim) ----
    allocate(ae(ne))
    ae = -999.0_WP
    ae(1:partit%myDim_elem2D) = real(partit%myList_elem2D(1:partit%myDim_elem2D), WP)
    call exchange_elem(ae, partit)
    k = partit%myDim_elem2D + partit%eDim_elem2D
    call check(maxval(abs(ae(1:k) - real(partit%myList_elem2D(1:k), WP))) == 0.0_WP, &
               'elem 2D halo-identity (small halo)')

    call MPI_Allreduce(nfail, nfail_g, 1, MPI_INTEGER, MPI_SUM, partit%MPI_COMM_FESOM, ierr)
    if (partit%mype == 0) then
        if (nfail_g == 0) then
            write(*,'(a,i0,a)') 'test_halo (npes=', partit%npes, '): ALL PASS'
        else
            write(*,'(a,i0,a,i0,a)') 'test_halo (npes=', partit%npes, '): ', nfail_g, ' FAILURE(S)'
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

    logical function node3d_ok(a3, partit)
        real(kind=WP),  intent(in) :: a3(:,:)
        type(t_partit), intent(in) :: partit
        integer :: j
        node3d_ok = .true.
        do j = 1, size(a3, 2)
            if (maxval(abs(a3(:, j) - real(partit%myList_nod2D(j), WP))) /= 0.0_WP) node3d_ok = .false.
        end do
    end function
end program test_halo
