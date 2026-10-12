module mod_shared_elem
    ! Shared elements (docs/plans/completed/20261012-shared-element-sync.md). FESOM2's partition rule
    ! gives an element to every rank that owns one of its vertices, so an element on a
    ! partition border has two or three OWNED copies. Each owner updates its own copy of a
    ! prognostic element field (uv); an order-dependent sum (the viscosity filter) makes the
    ! copies differ by roundoff, and the two ranks then compute different fluxes across
    ! their common edges -- the global divergence sum ssh_rhs_old is no longer zero and the
    ! linear free surface leaks (CORE2 results7/11/17/21, docs/HANDOFF.md).
    !
    ! build_com_elem2D_shared gives every shared element one DESIGNATED owner, the highest
    ! rank among its vertex owners (the copy the restart and output writers keep,
    ! mod_io_decomp decomp_redistribute_2d_r: last write wins), and builds
    ! partit%com_elem2D_shared: the designated owner sends the element to every other
    ! owner (send blocks), which overwrites its copy (receive blocks, local numbers
    ! 1..myDim_elem2D). Both sides build their blocks from the same owner sets, ordered by
    ! global element id with partners ascending, so they agree without communication; a
    ! one-time self-check proves it. The vertex owners come from com_nod2D, whose receive
    ! list covers every eDim halo node -- no MPI needed for the lists themselves.
    use mpi
    use mod_mesh,         only: t_mesh, MAX_NV
    use mod_partit,       only: t_partit, com_struct, MAX_NEIGHBOR_PARTITIONS
    use mod_partitioning, only: die
    use mod_halo,         only: exchange_elem_shared, allreduce_sum, allreduce_max
    implicit none
    private
    public :: build_com_elem2D_shared

contains

    subroutine build_com_elem2D_shared(mesh, partit)
        ! Build partit%com_elem2D_shared (re-entrant: the old lists are dropped first) and
        ! log the counts on rank 0. Every error stop goes through die (MPI_Abort).
        type(t_mesh),   intent(in)    :: mesh
        type(t_partit), intent(inout) :: partit
        integer :: npes, mype, nNodO, nNodL, nElemO
        integer :: i, j, k, p, e, d, nS, pass, ierr
        integer :: S(MAX_NV)
        integer :: nshared, ncopy, rmax, smax
        integer, allocatable :: own(:), order(:), nsend(:), nrecv(:), peer(:)
        integer, allocatable :: pos_s(:), pos_r(:), chk(:)
        character(len=256) :: msg
        character(len=16)  :: env
        integer :: env_len, ios

        npes   = partit%npes
        mype   = partit%mype
        nNodO  = partit%myDim_nod2D
        nNodL  = nNodO + partit%eDim_nod2D
        nElemO = partit%myDim_elem2D
        call reset_com(partit%com_elem2D_shared)

        ! ---- FESOM3_SHARED_SYNC: unset or 1 = update_vel syncs uv (default), 0 = off (the
        !      old unsynced copies, for A/B comparisons; the list is built anyway, so the
        !      copy-spread diagnostic still shows the drift). Value-based like
        !      FESOM3_MOMADV_OPT; any other value stops. No every-N option: every step with
        !      unequal copies leaks. ----
        partit%sync_shared_elem = .true.
        call get_environment_variable('FESOM3_SHARED_SYNC', env, length=env_len, status=ios)
        if (ios == 0 .and. env_len > 0) then
            select case (trim(adjustl(env)))
            case ('1')
                partit%sync_shared_elem = .true.
            case ('0')
                partit%sync_shared_elem = .false.
            case default
                call die('FESOM3_SHARED_SYNC must be 0 or 1', partit)
            end select
        end if

        ! ---- owner of every owned and eDim node (com_nod2D receive blocks) ----
        allocate(own(nNodL))
        own = -1
        own(1:nNodO) = mype
        do p = 1, partit%com_nod2D%rPEnum
            do k = partit%com_nod2D%rptr(p), partit%com_nod2D%rptr(p+1) - 1
                i = partit%com_nod2D%rlist(k)
                if (i < 1 .or. i > nNodL) &
                    call die('shared-element list: com_nod2D receives a node outside '// &
                             '1..myDim+eDim_nod2D', partit)
                own(i) = partit%com_nod2D%rPE(p)
            end do
        end do

        ! ---- owner sets: pass 1 counts per partner rank, pass 2 fills the blocks ----
        ! Owned elements in ascending global id, so every block comes out sorted.
        allocate(order(nElemO))
        call sort_by_key(partit%myList_elem2D(1:nElemO), order)
        allocate(nsend(0:npes-1), nrecv(0:npes-1), pos_s(0:npes-1), pos_r(0:npes-1))
        nsend = 0; nrecv = 0
        nshared = 0
        do pass = 1, 2
            if (pass == 2) call alloc_blocks()
            do i = 1, nElemO
                e = order(i)
                call owner_set(e)
                if (nS == 1) cycle
                d = maxval(S(1:nS))                      ! designated owner
                if (d == mype) then
                    if (pass == 1) nshared = nshared + 1
                    do j = 1, nS
                        p = S(j)
                        if (p == mype) cycle
                        if (pass == 1) then
                            nsend(p) = nsend(p) + 1
                        else
                            partit%com_elem2D_shared%slist(pos_s(p)) = e
                            pos_s(p) = pos_s(p) + 1
                        end if
                    end do
                else
                    if (pass == 1) then
                        nrecv(d) = nrecv(d) + 1
                    else
                        partit%com_elem2D_shared%rlist(pos_r(d)) = e
                        pos_r(d) = pos_r(d) + 1
                    end if
                end if
            end do
        end do

        ! ---- self-check, once: the counts agree pairwise (a mismatch would hang the
        !      exchange), then every received slot holds its own global id ----
        allocate(peer(0:npes-1))
        call MPI_Alltoall(nsend, 1, MPI_INTEGER, peer, 1, MPI_INTEGER, partit%MPI_COMM_FESOM, ierr)
        do p = 0, npes - 1
            if (peer(p) /= nrecv(p)) then
                write(msg,'(a,i0,a,i0,a,i0,a,i0)') 'shared-element list: rank ', p, ' sends ', &
                    peer(p), ' elements to rank ', mype, ', which expects ', nrecv(p)
                call die(trim(msg), partit)
            end if
        end do
        allocate(chk(nElemO))
        chk = partit%myList_elem2D(1:nElemO)
        chk(partit%com_elem2D_shared%rlist) = -1        ! receive slots are never sent
        call exchange_elem_shared(chk, partit)
        do p = 1, partit%com_elem2D_shared%rPEnum
            do k = partit%com_elem2D_shared%rptr(p), partit%com_elem2D_shared%rptr(p+1) - 1
                e = partit%com_elem2D_shared%rlist(k)
                if (chk(e) /= partit%myList_elem2D(e)) then
                    write(msg,'(a,i0,a,i0,a,i0,a,i0)') 'shared-element list: self-check failed '// &
                        'between rank ', partit%com_elem2D_shared%rPE(p), ' and rank ', mype, &
                        ': element ', partit%myList_elem2D(e), ' received ', chk(e)
                    call die(trim(msg), partit)
                end if
            end do
        end do

        ! ---- log: n counted once (at the designated owner), copies and maxima per rank ----
        ncopy = size(partit%com_elem2D_shared%rlist)
        rmax  = ncopy
        smax  = size(partit%com_elem2D_shared%slist)
        call allreduce_sum(nshared, partit)
        call allreduce_sum(ncopy, partit)
        call allreduce_max(rmax, partit)
        call allreduce_max(smax, partit)
        if (mype == 0) write(*,'(a,i0,a,i0,a,i0,a,i0)') 'shared-element sync: '// &
            trim(merge('ON ', 'OFF', partit%sync_shared_elem))//', ', nshared, &
            ' shared elements, ', ncopy, ' extra copies, max per rank: recv ', rmax, ', send ', smax

    contains

        subroutine owner_set(el)
            ! S(1:nS) = the distinct owners of el's vertices; el must have an owned vertex
            integer, intent(in) :: el
            integer :: kk, v
            nS = 0
            do kk = 1, mesh%elem2D_nnodes(el)
                v = mesh%elem2D_nodes(kk, el)
                if (v >= 1 .and. v <= nNodL) then
                    if (own(v) >= 0) then
                        if (all(S(1:nS) /= own(v))) then
                            nS = nS + 1
                            S(nS) = own(v)
                        end if
                        cycle
                    end if
                end if
                write(msg,'(a,i0,a)') 'shared-element list: element ', partit%myList_elem2D(el), &
                    ' has a vertex with no known owner'
                call die(trim(msg), partit)
            end do
            if (all(S(1:nS) /= mype)) then
                write(msg,'(a,i0,a,i0)') 'shared-element list: owned element ', &
                    partit%myList_elem2D(el), ' has no vertex owned by rank ', mype
                call die(trim(msg), partit)
            end if
        end subroutine owner_set

        subroutine alloc_blocks()
            ! partners ascending; rptr(1) = sptr(1) = 1 as read_com
            integer :: q, n
            character(len=16) :: cnum, cmax
            associate(c => partit%com_elem2D_shared)
                c%rPEnum = count(nrecv > 0)
                c%sPEnum = count(nsend > 0)
                if (max(c%rPEnum, c%sPEnum) > MAX_NEIGHBOR_PARTITIONS) then
                    write(cnum,'(i0)') max(c%rPEnum, c%sPEnum)
                    write(cmax,'(i0)') MAX_NEIGHBOR_PARTITIONS
                    call die('com_elem2D_shared: '//trim(cnum)//' partner PEs, outside '// &
                             '0..MAX_NEIGHBOR_PARTITIONS='//trim(cmax)//' (increase '// &
                             'MAX_NEIGHBOR_PARTITIONS in src/types/mod_partit.F90 and in '// &
                             'fesom_ini, then recompile)', partit)
                end if
                n = 0
                do q = 0, npes - 1
                    if (nrecv(q) == 0) cycle
                    n = n + 1
                    c%rPE(n) = q
                    pos_r(q) = c%rptr(n)
                    c%rptr(n+1) = c%rptr(n) + nrecv(q)
                end do
                n = 0
                do q = 0, npes - 1
                    if (nsend(q) == 0) cycle
                    n = n + 1
                    c%sPE(n) = q
                    pos_s(q) = c%sptr(n)
                    c%sptr(n+1) = c%sptr(n) + nsend(q)
                end do
                deallocate(c%rlist, c%slist)            ! the empty lists of reset_com
                allocate(c%rlist(c%rptr(c%rPEnum+1) - 1), c%slist(c%sptr(c%sPEnum+1) - 1))
            end associate
        end subroutine alloc_blocks

    end subroutine build_com_elem2D_shared

    subroutine reset_com(com)
        ! empty list (no partners), lists allocated with size 0
        type(com_struct), intent(inout) :: com
        com%rPEnum = 0; com%sPEnum = 0; com%nreq = 0
        com%rptr(1) = 1; com%sptr(1) = 1
        if (allocated(com%rlist)) deallocate(com%rlist)
        if (allocated(com%slist)) deallocate(com%slist)
        allocate(com%rlist(0), com%slist(0))
    end subroutine reset_com

    subroutine sort_by_key(key, idx)
        ! idx = the permutation of 1..size(key) with key(idx) ascending (heapsort; the keys,
        ! global element ids, are distinct)
        integer, intent(in)  :: key(:)
        integer, intent(out) :: idx(:)
        integer :: n, i, t
        n = size(key)
        idx = [(i, i = 1, n)]
        do i = n/2, 1, -1
            call sift(i, n)
        end do
        do i = n, 2, -1
            t = idx(1); idx(1) = idx(i); idx(i) = t
            call sift(1, i - 1)
        end do
    contains
        subroutine sift(first, last)
            integer, intent(in) :: first, last
            integer :: root, child, tt
            root = first
            do
                child = 2*root
                if (child > last) exit
                if (child < last) then
                    if (key(idx(child+1)) > key(idx(child))) child = child + 1
                end if
                if (key(idx(root)) >= key(idx(child))) exit
                tt = idx(root); idx(root) = idx(child); idx(child) = tt
                root = child
            end do
        end subroutine sift
    end subroutine sort_by_key

end module mod_shared_elem
