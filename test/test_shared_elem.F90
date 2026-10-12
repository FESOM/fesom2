program test_shared_elem
    ! Shared-element sync (docs/plans/completed/20261012-shared-element-sync.md): the list
    ! partit%com_elem2D_shared built at the end of read_mesh_local (src/infra/mod_shared_elem.F90).
    !
    !   S1  the list against an independent oracle: a global owner map from an allgather of
    !       the owned node ids gives every owned element's owner set S and designated owner
    !       d = max(S). An element is in the receive block of d exactly when |S| > 1 and this
    !       rank is not d; it is in the send block of every other owner exactly when this
    !       rank is d. Blocks ascend by global id, partners ascend, every shared element has
    !       exactly one designated owner, and every rank in S owns the element (the partition
    !       rule the design relies on). At np 1 the list is empty.
    !   S2  exchange_elem_shared and shared_copy_max_elem. A (2,nl-1,elem) array and the
    !       (1,:,:,:) section of a (1,2,nl-1,elem) array hold blkval(mype, global id, j1, j2):
    !       before the exchange the block probe is positive (np > 1); after it every receive
    !       slot holds blkval(d, ...) of its designated owner d, the probe is 0, and the other
    !       owned and halo slots are bitwise unchanged. Rank-1 probe: a difference planted on
    !       one receive slot is returned exactly on that rank and 0 elsewhere. At np 1 the
    !       arrays stay unchanged and the probes return 0.
    !   S3  the mechanism, through the real kernels update_vel and compute_hbar_ale. uv and
    !       uv_rhs are hashed from the global id (equal on every copy), DELTA is added to
    !       uv_rhs on this rank's receive slots (the non-designated copies), d_eta = 0. With
    !       the sync on, the divergence-sum ratio r = |sum ssh_rhs_old| / sum|ssh_rhs_old|
    !       over owned nodes is < 1e-12 and the uv copy spread is 0; with it off (np > 1)
    !       r > 1e-8 and the spread is DELTA. At np 1, r < 1e-12 either way.
    !
    ! FESOM3_SHARED_SYNC (read by the builder): partit%sync_shared_elem must follow it; an
    ! invalid value must stop the builder (ctest regex pair with "did NOT stop").
    ! FESOM3_MESH_DIR (default pi), FESOM3_TEST_ROTATE=1 reads it rotated (setup_test_mesh).
    ! FESOM3_TEST_SHARED_CORRUPT=1 (np > 1): clear the com_nod2D receive blocks and build the
    ! list again; the builder must stop with "has a vertex with no known owner" (ctest regex
    ! pair with "did NOT stop").
    use mpi
    use mod_precision,    only: WP
    use mod_mesh,         only: t_mesh
    use mod_partit,       only: t_partit
    use mod_partitioning, only: par_init
    use mod_dyn,          only: t_dyn
    use mod_part_bounds,  only: owned_bounds
    use mod_shared_elem,  only: build_com_elem2D_shared
    use mod_halo,         only: exchange_elem_shared, shared_copy_max_elem
    use oce_ale,          only: update_vel, compute_hbar_ale
    use test_util,        only: setup_test_mesh, check, finish, sum_all, max_all, hval
    implicit none

    character(len=*), parameter :: PI_DIR = '/home/a/a270088/port2/fesom2/tests/data/MESHES/pi'
    type(t_partit) :: partit
    type(t_mesh)   :: mesh
    character(len=16) :: cbuf
    integer, allocatable :: gown(:)       ! oracle: owner rank of every global node

    call par_init(partit)
    call setup_test_mesh(mesh, partit, PI_DIR, cartesian=.false., partition_np1=.true.)

    cbuf = ''
    call get_environment_variable('FESOM3_SHARED_SYNC', cbuf)
    if (partit%npes > 1) then
        select case (trim(cbuf))
        case ('', '1')
            call check(partit%sync_shared_elem, 'FESOM3_SHARED_SYNC unset or 1: sync on', partit%mype)
        case ('0')
            call check(.not. partit%sync_shared_elem, 'FESOM3_SHARED_SYNC=0: sync off', partit%mype)
        case default
            write(*,'(a)') 'test_shared_elem: FESOM3_SHARED_SYNC='//trim(cbuf)// &
                           ': the builder did NOT stop'
        end select
    end if

    cbuf = ''
    call get_environment_variable('FESOM3_TEST_SHARED_CORRUPT', cbuf)
    if (trim(cbuf) == '1') then
        partit%com_nod2D%rPEnum = 0                     ! no eDim node has a known owner
        call build_com_elem2D_shared(mesh, partit)
        write(*,'(a)') 'test_shared_elem: corrupt com_nod2D: the builder did NOT stop'
        call finish(partit, 'test_shared_elem')
    end if

    call s1_list()
    call s2_exchange()
    call s3_mechanism()
    call finish(partit, 'test_shared_elem')

contains

    subroutine s1_list()
        integer :: nElemO, ngnod, e, k, p, q, r, nS, d, ierr, gid, nshared
        integer :: S(4)
        integer, allocatable :: cnts(:), displs(:), gids(:)
        integer, allocatable :: recv_from(:), recv_occ(:), send_cnt(:), send_bad(:)
        integer, allocatable :: ndesig(:), nowner(:)
        logical :: ok_order, ok_range

        nElemO = partit%myDim_elem2D
        if (partit%npes == 1) then
            call check(partit%com_elem2D_shared%rPEnum == 0 .and. &
                       partit%com_elem2D_shared%sPEnum == 0, 'S1 np 1: empty list', partit%mype)
            return
        end if

        ! ---- oracle: global owner map from the owned node ids of every rank ----
        allocate(cnts(0:partit%npes-1), displs(0:partit%npes-1))
        call MPI_Allgather(partit%myDim_nod2D, 1, MPI_INTEGER, cnts, 1, MPI_INTEGER, &
                           partit%MPI_COMM_FESOM, ierr)
        displs(0) = 0
        do p = 1, partit%npes - 1
            displs(p) = displs(p-1) + cnts(p-1)
        end do
        ngnod = sum(cnts)
        call check(ngnod == mesh%nod2D, 'S1 owned nodes cover the mesh once', partit%mype)
        allocate(gids(ngnod), gown(mesh%nod2D))
        call MPI_Allgatherv(partit%myList_nod2D, partit%myDim_nod2D, MPI_INTEGER, gids, cnts, &
                            displs, MPI_INTEGER, partit%MPI_COMM_FESOM, ierr)
        gown = -1
        do p = 0, partit%npes - 1
            gown(gids(displs(p)+1:displs(p)+cnts(p))) = p
        end do
        call check(all(gown >= 0), 'S1 every node has an owner', partit%mype)

        ! ---- the list as built: membership per owned element ----
        allocate(recv_from(nElemO), recv_occ(nElemO), send_cnt(nElemO), send_bad(nElemO))
        recv_from = -1; recv_occ = 0; send_cnt = 0; send_bad = 0
        ok_order = .true.; ok_range = .true.
        associate(c => partit%com_elem2D_shared)
            do p = 1, c%rPEnum
                if (p > 1) ok_order = ok_order .and. c%rPE(p) > c%rPE(p-1)
                do k = c%rptr(p), c%rptr(p+1) - 1
                    e = c%rlist(k)
                    if (e < 1 .or. e > nElemO) then
                        ok_range = .false.; cycle
                    end if
                    recv_from(e) = c%rPE(p); recv_occ(e) = recv_occ(e) + 1
                    if (k > c%rptr(p)) ok_order = ok_order .and. &
                        partit%myList_elem2D(e) > partit%myList_elem2D(c%rlist(k-1))
                end do
            end do
            do p = 1, c%sPEnum
                if (p > 1) ok_order = ok_order .and. c%sPE(p) > c%sPE(p-1)
                do k = c%sptr(p), c%sptr(p+1) - 1
                    e = c%slist(k)
                    if (e < 1 .or. e > nElemO) then
                        ok_range = .false.; cycle
                    end if
                    send_cnt(e) = send_cnt(e) + 1
                    call owners(e, S, nS)
                    if (all(S(1:nS) /= c%sPE(p)) .or. c%sPE(p) == partit%mype) send_bad(e) = 1
                    if (k > c%sptr(p)) ok_order = ok_order .and. &
                        partit%myList_elem2D(e) > partit%myList_elem2D(c%slist(k-1))
                end do
            end do
        end associate
        call check(ok_range, 'S1 list entries are owned elements', partit%mype)
        call check(ok_order, 'S1 blocks ascend by global id, partners ascend', partit%mype)

        ! ---- compare with the oracle ----
        allocate(ndesig(mesh%elem2D), nowner(mesh%elem2D))
        ndesig = 0; nowner = 0
        nshared = 0
        q = 0; r = 0
        do e = 1, nElemO
            gid = partit%myList_elem2D(e)
            call owners(e, S, nS)
            d = maxval(S(1:nS))
            nowner(gid) = 1
            if (nS > 1 .and. d == partit%mype) then
                ndesig(gid) = 1
                nshared = nshared + 1
            end if
            if (any(S(1:nS) == partit%mype)) then
                if (nS > 1 .and. d /= partit%mype) then
                    if (recv_occ(e) /= 1 .or. recv_from(e) /= d) q = q + 1
                    if (send_cnt(e) /= 0) q = q + 1
                else if (nS > 1) then
                    if (recv_occ(e) /= 0) q = q + 1
                    if (send_cnt(e) /= nS - 1 .or. send_bad(e) /= 0) q = q + 1
                else
                    if (recv_occ(e) /= 0 .or. send_cnt(e) /= 0) q = q + 1
                end if
            else
                r = r + 1                                  ! owned without an owned vertex
            end if
        end do
        call check(r == 0, 'S1 every owned element has an owned vertex', partit%mype)
        call check(q == 0, 'S1 receive/send membership == oracle', partit%mype)
        call MPI_Allreduce(MPI_IN_PLACE, ndesig, mesh%elem2D, MPI_INTEGER, MPI_SUM, &
                           partit%MPI_COMM_FESOM, ierr)
        call MPI_Allreduce(MPI_IN_PLACE, nowner, mesh%elem2D, MPI_INTEGER, MPI_SUM, &
                           partit%MPI_COMM_FESOM, ierr)
        q = 0; r = 0
        do e = 1, nElemO
            gid = partit%myList_elem2D(e)
            call owners(e, S, nS)
            if (nowner(gid) /= nS) r = r + 1             ! every vertex owner owns e
            if (nS > 1 .and. ndesig(gid) /= 1) q = q + 1
            if (nS == 1 .and. ndesig(gid) /= 0) q = q + 1
        end do
        call check(r == 0, 'S1 every vertex owner owns the element', partit%mype)
        call check(q == 0, 'S1 one designated owner per shared element', partit%mype)
        call check(all(nowner >= 1), 'S1 every element is owned', partit%mype)
        call sum_all(nshared, partit)
        call check(nshared > 0, 'S1 the partition has shared elements', partit%mype)
        if (partit%mype == 0) write(*,'(a,i0,a,i0)') '  S1: shared elements ', nshared, &
                                                     ' of ', mesh%elem2D
    end subroutine s1_list

    subroutine s2_exchange()
        integer :: nElemO, nElemF, nl, e, j1, j2, p, k, d, nbad, rplant, ierr
        real(kind=WP), allocatable :: a3(:,:,:), a3_0(:,:,:), ab(:,:,:,:), ab_0(:,:,:,:)
        real(kind=WP), allocatable :: a1(:)
        logical, allocatable :: is_recv(:)
        integer, allocatable :: dsg(:)
        real(kind=WP) :: pr, pr_ab
        real(kind=WP), parameter :: PLANT = 0.375_WP

        nElemO = partit%myDim_elem2D
        nElemF = nElem_full()
        nl = mesh%nl
        allocate(a3(2, nl-1, nElemF), ab(1, 2, nl-1, nElemF))
        do e = 1, nElemF
            do j2 = 1, nl - 1
                do j1 = 1, 2
                    a3(j1, j2, e)    = blkval(partit%mype, partit%myList_elem2D(e), j1, j2)
                    ab(1, j1, j2, e) = blkval(partit%mype, partit%myList_elem2D(e), j1, j2)
                end do
            end do
        end do
        a3_0 = a3; ab_0 = ab
        ! designated owner of every receive slot (the block's partner)
        allocate(is_recv(nElemF), dsg(nElemF))
        is_recv = .false.; dsg = -1
        do p = 1, partit%com_elem2D_shared%rPEnum
            do k = partit%com_elem2D_shared%rptr(p), partit%com_elem2D_shared%rptr(p+1) - 1
                e = partit%com_elem2D_shared%rlist(k)
                is_recv(e) = .true.; dsg(e) = partit%com_elem2D_shared%rPE(p)
            end do
        end do

        pr = shared_copy_max_elem(a3, partit); call max_all(pr, partit)
        pr_ab = shared_copy_max_elem(ab(1,:,:,:), partit); call max_all(pr_ab, partit)
        if (partit%npes > 1) then
            call check(pr > 0.0_WP .and. pr_ab > 0.0_WP, 'S2 probe > 0 before the exchange', partit%mype)
        else
            call check(pr == 0.0_WP .and. pr_ab == 0.0_WP, 'S2 np 1: probe 0', partit%mype)
        end if

        call exchange_elem_shared(a3, partit)
        call exchange_elem_shared(ab(1,:,:,:), partit)
        nbad = 0
        do e = 1, nElemF
            if (is_recv(e)) then
                d = dsg(e)
                do j2 = 1, nl - 1
                    do j1 = 1, 2
                        if (a3(j1, j2, e) /= blkval(d, partit%myList_elem2D(e), j1, j2)) nbad = nbad + 1
                        if (ab(1, j1, j2, e) /= blkval(d, partit%myList_elem2D(e), j1, j2)) nbad = nbad + 1
                    end do
                end do
            else
                if (any(a3(:,:,e) /= a3_0(:,:,e)) .or. any(ab(:,:,:,e) /= ab_0(:,:,:,e))) nbad = nbad + 1
            end if
        end do
        call check(nbad == 0, 'S2 receive slots hold the designated copy, the rest unchanged', partit%mype)
        pr = shared_copy_max_elem(a3, partit); call max_all(pr, partit)
        pr_ab = shared_copy_max_elem(ab(1,:,:,:), partit); call max_all(pr_ab, partit)
        call check(pr == 0.0_WP .and. pr_ab == 0.0_WP, 'S2 probe 0 after the exchange', partit%mype)
        if (partit%npes > 1) call check(count(is_recv) > 0 .or. partit%com_elem2D_shared%sPEnum > 0, &
                                        'S2 every rank takes part in the exchange', partit%mype)

        ! ---- rank-1 probe: plant PLANT on one receive slot of the lowest rank that has one ----
        allocate(a1(nElemF))
        a1 = real(partit%myList_elem2D(1:nElemF), WP)
        rplant = partit%npes
        if (count(is_recv) > 0) rplant = partit%mype
        if (partit%npes > 1) call MPI_Allreduce(MPI_IN_PLACE, rplant, 1, MPI_INTEGER, MPI_MIN, &
                                                partit%MPI_COMM_FESOM, ierr)
        if (partit%mype == rplant) a1(partit%com_elem2D_shared%rlist(1)) = &
                                   a1(partit%com_elem2D_shared%rlist(1)) + PLANT
        pr = shared_copy_max_elem(a1, partit)
        if (partit%mype == rplant) then
            call check(pr == PLANT, 'S2 rank-1 probe returns the planted difference', partit%mype)
        else
            call check(pr == 0.0_WP, 'S2 rank-1 probe 0 on the other ranks', partit%mype)
        end if
        if (partit%npes > 1) call check(rplant < partit%npes, 'S2 some rank has a receive slot', &
                                        partit%mype)
        if (partit%mype == 0) write(*,'(a)') '  S2: exchange and probes done'
    end subroutine s2_exchange

    subroutine s3_mechanism()
        type(t_dyn) :: dyn
        integer :: nNodO, nNodL, nEdgeO, nElemO, nElemF, nl, e, nz, gid, k, mode
        real(kind=WP), allocatable :: uv0(:,:,:)
        real(kind=WP) :: r(2), pr(2)
        logical :: sync_env
        real(kind=WP), parameter :: DELTA = 1.0e-3_WP, DT = 1800.0_WP

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        nElemF = nElem_full()
        nl = mesh%nl
        allocate(dyn%uv(2, nl-1, nElemF), dyn%uv_rhs(2, nl-1, nElemF))
        allocate(dyn%d_eta(nNodL), dyn%ssh_rhs_old(nNodL))
        allocate(mesh%helem(nl-1, nElemF), mesh%hbar(nNodL), mesh%hbar_old(nNodL), mesh%dhe(nElemF))
        do e = 1, nElemF
            gid = partit%myList_elem2D(e)
            do nz = 1, nl - 1
                mesh%helem(nz, e)  = mesh%zbar(nz) - mesh%zbar(nz+1)
                dyn%uv(1, nz, e)     = hval(gid, nz)
                dyn%uv(2, nz, e)     = hval(gid, nz + 97)
                dyn%uv_rhs(1, nz, e) = 0.01_WP*hval(gid, nz + 211)
                dyn%uv_rhs(2, nz, e) = 0.01_WP*hval(gid, nz + 307)
            end do
        end do
        do k = 1, size(partit%com_elem2D_shared%rlist)      ! the non-designated copies
            e = partit%com_elem2D_shared%rlist(k)
            dyn%uv_rhs(:, :, e) = dyn%uv_rhs(:, :, e) + DELTA
        end do
        dyn%d_eta = 0.0_WP
        uv0 = dyn%uv

        sync_env = partit%sync_shared_elem
        do mode = 1, 2                                      ! 1: sync on, 2: off
            partit%sync_shared_elem = (mode == 1)
            dyn%uv = uv0
            mesh%hbar = 0.0_WP; mesh%hbar_old = 0.0_WP; mesh%dhe = 0.0_WP
            call update_vel(dyn, mesh, DT, partit)
            call compute_hbar_ale(dyn, mesh, DT, partit)
            r(mode) = divsum_ratio(dyn%ssh_rhs_old(1:nNodO))
            pr(mode) = shared_copy_max_elem(dyn%uv, partit)
            call max_all(pr(mode), partit)
        end do
        partit%sync_shared_elem = sync_env

        call check(r(1) < 1.0e-12_WP, 'S3 sync on: divergence-sum ratio < 1e-12', partit%mype)
        call check(pr(1) == 0.0_WP, 'S3 sync on: uv copy spread 0', partit%mype)
        if (partit%npes > 1) then
            call check(r(2) > 1.0e-8_WP, 'S3 sync off: divergence-sum ratio > 1e-8', partit%mype)
            call check(abs(pr(2) - DELTA) < 1.0e-9_WP, 'S3 sync off: uv copy spread DELTA', partit%mype)
        else
            call check(r(2) < 1.0e-12_WP, 'S3 np 1: divergence-sum ratio < 1e-12', partit%mype)
        end if
        if (partit%mype == 0) write(*,'(a,2es10.2,a,2es10.2)') '  S3: ratio on/off', r, &
                                                               ', uv spread on/off', pr
    end subroutine s3_mechanism

    real(kind=WP) function divsum_ratio(rhs)
        ! |sum rhs| / sum|rhs| over the owned nodes of all ranks
        real(kind=WP), intent(in) :: rhs(:)
        real(kind=WP) :: s(2)
        s = [sum(rhs), sum(abs(rhs))]
        call sum_all(s, partit)
        divsum_ratio = abs(s(1)) / max(s(2), tiny(1.0_WP))
    end function divsum_ratio

    integer function nElem_full()
        ! local element count: owned + halo (+ double halo) at np > 1, all at np 1
        if (partit%npes == 1) then
            nElem_full = mesh%elem2D
        else
            nElem_full = partit%myDim_elem2D + partit%eDim_elem2D + partit%eXDim_elem2D
        end if
    end function nElem_full

    real(kind=WP) function blkval(rank, gid, j1, j2)
        ! rank-dependent value per element, component and level (exact in double precision)
        integer, intent(in) :: rank, gid, j1, j2
        blkval = 1.0e7_WP*real(rank, WP) + real(gid, WP) + 0.5_WP*real(j1, WP) + 0.03125_WP*real(j2, WP)
    end function blkval

    subroutine owners(e, S, nS)
        ! oracle owner set of owned element e from the global owner map gown
        integer, intent(in)  :: e
        integer, intent(out) :: S(4), nS
        integer :: k, o
        nS = 0
        do k = 1, mesh%elem2D_nnodes(e)
            o = gown(partit%myList_nod2D(mesh%elem2D_nodes(k, e)))
            if (all(S(1:nS) /= o)) then
                nS = nS + 1; S(nS) = o
            end if
        end do
    end subroutine owners

end program test_shared_elem
