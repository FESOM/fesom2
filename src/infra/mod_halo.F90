module mod_halo
    ! Generic MPI halo exchange (decision D7). The generic exchange_nod/exchange_elem
    ! interfaces + com_struct ARE the language-agnostic portability contract; MPI is
    ! the single implementation (manual pack -> Isend/Irecv -> unpack). Halo is
    ! BROADCAST-ONLY: each halo entry is OVERWRITTEN with its owner's value (never
    ! additive). rlist/slist hold LOCAL indices; the field array is indexed locally
    ! (1..myDim+eDim). Structure follows the dwarf gen_halo_exchange.F90, written
    ! clean (GPU/g2g/OpenACC bloat omitted; luse_g2g kept as an inert foundation arg).
    !
    ! Deferred (optimizations, not needed for byte-identity): begin/end overlap split,
    ! multi-field batched exchange, com_elem2D_full default path.
    use mpi
    use mod_precision, only: WP, MPI_WP
    use mod_partit,    only: t_partit, com_struct
    implicit none
    private
    public :: exchange_nod, exchange_elem, exchange_elem_full, &
              exchange_edge, exchange_nod_full, exchange_elem_shared, &
              allreduce_sum, allreduce_max, stale_halo_max_nod, stale_halo_max_elem, &
              shared_copy_max_elem

    interface exchange_nod
        module procedure exchange_nod_2D_r, exchange_nod_3D_r, exchange_nod_2D_i, &
                         exchange_nod_blk_r
    end interface

    ! Cross-rank floating-point SUM reduction (the CG dot-products). Mirrors the
    ! FESOM2 oracle solver.F90's MPI_Allreduce(MPI_IN_PLACE, s, n, MPI_DOUBLE, MPI_SUM,
    ! MPI_COMM_FESOM): same library + comm size + op + 8-byte type => OpenMPI picks the
    ! same reduction tree => byte-identical given byte-identical per-rank partial sums
    ! (the same determinism FESOM2 relies on, LESSONS L6). Guard callers with
    ! is_multirank — at npes==1 the local sum already IS the global sum (identity).
    ! The integer variant sums COUNTS (mesh invariants, split statistics): exact.
    interface allreduce_sum
        module procedure allreduce_sum_r0, allreduce_sum_r1, allreduce_sum_i0
    end interface

    ! Cross-rank floating-point MAX reduction. Mirrors the FESOM2 extrap_nod3D /
    ! do_ic3d MPI_AllREDUCE(loc_max, glob_max, MPI_MAX, MPI_COMM_FESOM): MAX is a
    ! selection (no rounding) so it is exact + deterministic; glob_max only drives the
    ! extrapolation loop COUNT (same global max => same iteration count on both codes).
    ! Guard callers with is_multirank — at npes==1 the local max already IS the global.
    ! The integer variant takes maxima of COUNTS (setup statistics).
    interface allreduce_max
        module procedure allreduce_max_r0, allreduce_max_i0
    end interface

    interface exchange_elem
        module procedure exchange_elem_2D_r, exchange_elem_3D_r, exchange_elem_blk_r
    end interface

    ! Full element halo (com_elem2D_full = eDim+eXDim) — the FESOM2 convention picks
    ! com_elem2D_full when an element array is sized myDim+eDim+eXDim (the eXDim second
    ! halo layer that MUSCL/find_neighbors need: every halo node's full element list is
    ! then local). Same broadcast-only pack/unpack as exchange_elem, on com_elem2D_full.
    interface exchange_elem_full
        module procedure exchange_elem_full_2D_r, exchange_elem_full_2D_i, &
                         exchange_elem_full_3D_r
    end interface

    ! Hybrid meshes (dist_N/format.out): edge halo and double node halo. exchange_edge
    ! refreshes every held edge 1..myDim+eDim_edge2D that another rank owns uniquely (owner
    ! of edges(1,edge)) -- also edges in this rank's own 1..myDim_edge2D range;
    ! exchange_nod_full refreshes the halo nodes 1..eDim+eXDim_nod2D (arrays sized
    ! myDim+eDim+eXDim = nNodX, owned_bounds). Same broadcast-only cores as exchange_nod;
    ! a dist_N without format.out has no such lists and the call stops with a message.
    interface exchange_edge
        module procedure exchange_edge_2D_r, exchange_edge_3D_r
    end interface

    interface exchange_nod_full
        module procedure exchange_nod_full_2D_r, exchange_nod_full_3D_r
    end interface

    ! Shared elements (partit%com_elem2D_shared, src/infra/mod_shared_elem.F90): the
    ! designated owner's copy of an element is sent to the element's other owners, which
    ! overwrite their own copy (local numbers 1..myDim_elem2D). Same broadcast-only cores;
    ! a no-op on an empty list, and so at 1 rank.
    interface exchange_elem_shared
        module procedure exchange_elem_shared_2D_i, exchange_elem_shared_blk_r
    end interface

    ! Copy spread of the shared elements: this rank's max|arr - the designated owner's copy|
    ! over its receive slots of com_elem2D_shared (exchange-and-compare, as the stale-halo
    ! probe). The value is LOCAL (0 on a rank without receive slots); reduce it with
    ! allreduce_max. 0 <=> every copy equals the designated one.
    interface shared_copy_max_elem
        module procedure shared_copy_max_elem_2D_r, shared_copy_max_elem_blk_r
    end interface

    integer, parameter :: HALO_TAG = 1
    integer, parameter :: BLK_TAG  = 2   ! distinct tag for the rank-3 block exchange (core_blk_r)

contains

    ! ===================== node wrappers =====================
    subroutine exchange_nod_2D_r(arr, partit, luse_g2g)
        real(kind=WP),  intent(inout) :: arr(:)
        type(t_partit), intent(in)    :: partit
        logical, optional, intent(in) :: luse_g2g   ! GPU-direct hook (inert in v1)
        call core_2D_r(arr, partit%com_nod2D, partit%MPI_COMM_FESOM)
    end subroutine

    subroutine exchange_nod_3D_r(arr, partit, luse_g2g)
        real(kind=WP),  intent(inout) :: arr(:,:)   ! (nl, nod_size)
        type(t_partit), intent(in)    :: partit
        logical, optional, intent(in) :: luse_g2g
        call core_3D_r(arr, partit%com_nod2D, partit%MPI_COMM_FESOM)
    end subroutine

    subroutine exchange_nod_2D_i(arr, partit, luse_g2g)
        integer,        intent(inout) :: arr(:)
        type(t_partit), intent(in)    :: partit
        logical, optional, intent(in) :: luse_g2g
        call core_2D_i(arr, partit%com_nod2D, partit%MPI_COMM_FESOM)
    end subroutine

    subroutine exchange_nod_blk_r(arr, partit, luse_g2g)
        ! (d1, d2, nod_size) node field, e.g. UVnode/UVnode_rhs (2, nl-1, nod) — the
        ! leading two dims form a contiguous block transported per halo NODE (FESOM2's
        ! exchange_nod on a (2,nz,n) vector-at-nodes field). Broadcast-only on com_nod2D.
        real(kind=WP),  intent(inout) :: arr(:,:,:)
        type(t_partit), intent(in)    :: partit
        logical, optional, intent(in) :: luse_g2g
        call core_blk_r(arr, partit%com_nod2D, partit%MPI_COMM_FESOM)
    end subroutine

    ! ===================== element wrappers ==================
    subroutine exchange_elem_2D_r(arr, partit, luse_g2g)
        real(kind=WP),  intent(inout) :: arr(:)
        type(t_partit), intent(in)    :: partit
        logical, optional, intent(in) :: luse_g2g
        call core_2D_r(arr, partit%com_elem2D, partit%MPI_COMM_FESOM)
    end subroutine

    subroutine exchange_elem_3D_r(arr, partit, luse_g2g)
        real(kind=WP),  intent(inout) :: arr(:,:)
        type(t_partit), intent(in)    :: partit
        logical, optional, intent(in) :: luse_g2g
        call core_3D_r(arr, partit%com_elem2D, partit%MPI_COMM_FESOM)
    end subroutine

    ! Rank-3 (d1,d2,elem) element exchange over the STANDARD eDim halo (com_elem2D) —
    ! the FESOM2 exchange_elem(UV) for the (2,nl-1,elem) velocity. NOT com_elem2D_full:
    ! FESOM2's update_vel refreshes only the eDim halo, which is all the owned-edge
    ! kernels (viscosity etc.) read. core_blk_r is the same block core the node-block
    ! exchange_nod(UVnode) uses (proven), here driven by com_elem2D.
    subroutine exchange_elem_blk_r(arr, partit, luse_g2g)
        real(kind=WP),  intent(inout) :: arr(:,:,:)
        type(t_partit), intent(in)    :: partit
        logical, optional, intent(in) :: luse_g2g
        call core_blk_r(arr, partit%com_elem2D, partit%MPI_COMM_FESOM)
    end subroutine

    ! --------------- full element halo wrappers (com_elem2D_full) -----------
    subroutine exchange_elem_full_2D_r(arr, partit, luse_g2g)
        real(kind=WP),  intent(inout) :: arr(:)
        type(t_partit), intent(in)    :: partit
        logical, optional, intent(in) :: luse_g2g
        call core_2D_r(arr, partit%com_elem2D_full, partit%MPI_COMM_FESOM)
    end subroutine

    subroutine exchange_elem_full_2D_i(arr, partit, luse_g2g)
        integer,        intent(inout) :: arr(:)
        type(t_partit), intent(in)    :: partit
        logical, optional, intent(in) :: luse_g2g
        call core_2D_i(arr, partit%com_elem2D_full, partit%MPI_COMM_FESOM)
    end subroutine

    subroutine exchange_elem_full_3D_r(arr, partit, luse_g2g)
        real(kind=WP),  intent(inout) :: arr(:,:,:)   ! (d1, d2, entity) e.g. tr_xy(2,nl-1,*)
        type(t_partit), intent(in)    :: partit
        logical, optional, intent(in) :: luse_g2g
        call core_blk_r(arr, partit%com_elem2D_full, partit%MPI_COMM_FESOM)
    end subroutine

    ! --------------- shared-element wrappers (com_elem2D_shared) -----------
    subroutine exchange_elem_shared_2D_i(arr, partit)
        integer,        intent(inout) :: arr(:)       ! (elem), e.g. the builder's self-check
        type(t_partit), intent(in)    :: partit
        call core_2D_i(arr, partit%com_elem2D_shared, partit%MPI_COMM_FESOM)
    end subroutine

    subroutine exchange_elem_shared_blk_r(arr, partit)
        real(kind=WP),  intent(inout) :: arr(:,:,:)   ! (d1, d2, elem): uv (2,nl-1,elem) or
        type(t_partit), intent(in)    :: partit       ! the section uv_rhsAB(1,:,:,:)
        call core_blk_r(arr, partit%com_elem2D_shared, partit%MPI_COMM_FESOM)
    end subroutine

    ! --------------- hybrid-mesh wrappers (com_edge2D, com_nod2D_full) -----------
    subroutine exchange_edge_2D_r(arr, partit, luse_g2g)
        real(kind=WP),  intent(inout) :: arr(:)     ! (myDim+eDim_edge2D)
        type(t_partit), intent(in)    :: partit
        logical, optional, intent(in) :: luse_g2g
        call require_dist_format(partit, 'exchange_edge', 'com_edge2D')
        call core_2D_r(arr, partit%com_edge2D, partit%MPI_COMM_FESOM)
    end subroutine

    subroutine exchange_edge_3D_r(arr, partit, luse_g2g)
        real(kind=WP),  intent(inout) :: arr(:,:)   ! (nl, myDim+eDim_edge2D)
        type(t_partit), intent(in)    :: partit
        logical, optional, intent(in) :: luse_g2g
        call require_dist_format(partit, 'exchange_edge', 'com_edge2D')
        call core_3D_r(arr, partit%com_edge2D, partit%MPI_COMM_FESOM)
    end subroutine

    subroutine exchange_nod_full_2D_r(arr, partit, luse_g2g)
        real(kind=WP),  intent(inout) :: arr(:)     ! (nNodX)
        type(t_partit), intent(in)    :: partit
        logical, optional, intent(in) :: luse_g2g
        call require_dist_format(partit, 'exchange_nod_full', 'com_nod2D_full')
        call core_2D_r(arr, partit%com_nod2D_full, partit%MPI_COMM_FESOM)
    end subroutine

    subroutine exchange_nod_full_3D_r(arr, partit, luse_g2g)
        real(kind=WP),  intent(inout) :: arr(:,:)   ! (nl, nNodX)
        type(t_partit), intent(in)    :: partit
        logical, optional, intent(in) :: luse_g2g
        call require_dist_format(partit, 'exchange_nod_full', 'com_nod2D_full')
        call core_3D_r(arr, partit%com_nod2D_full, partit%MPI_COMM_FESOM)
    end subroutine

    subroutine require_dist_format(partit, caller, list)
        ! The appended lists exist only with dist_N/format.out (the 1-rank synthesis
        ! provides them empty); otherwise stop with a message instead of a silent no-op.
        type(t_partit),   intent(in) :: partit
        character(len=*), intent(in) :: caller, list
        integer :: ierr
        if (partit%dist_format >= 1) return
        write(*,'(a,i0,a)') caller//' [rank ', partit%mype, ']: '//list// &
            ' is not available: dist_N has no format.out (FESOM2 layout); '// &
            'regenerate dist_N with the fesom_ini of this repository'
        call MPI_Abort(partit%MPI_COMM_FESOM, 1, ierr)
    end subroutine require_dist_format

    ! ===================== allreduce (CG dot-products) =======
    subroutine allreduce_sum_r0(s, partit)
        real(kind=WP),  intent(inout) :: s
        type(t_partit), intent(in)    :: partit
        integer :: ierr
        call MPI_Allreduce(MPI_IN_PLACE, s, 1, MPI_WP, MPI_SUM, &
                           partit%MPI_COMM_FESOM, ierr)
    end subroutine allreduce_sum_r0

    subroutine allreduce_sum_r1(s, partit)
        real(kind=WP),  intent(inout) :: s(:)
        type(t_partit), intent(in)    :: partit
        integer :: ierr
        call MPI_Allreduce(MPI_IN_PLACE, s, size(s), MPI_WP, MPI_SUM, &
                           partit%MPI_COMM_FESOM, ierr)
    end subroutine allreduce_sum_r1

    subroutine allreduce_sum_i0(s, partit)
        integer,        intent(inout) :: s
        type(t_partit), intent(in)    :: partit
        integer :: ierr
        call MPI_Allreduce(MPI_IN_PLACE, s, 1, MPI_INTEGER, MPI_SUM, &
                           partit%MPI_COMM_FESOM, ierr)
    end subroutine allreduce_sum_i0

    subroutine allreduce_max_r0(s, partit)
        real(kind=WP),  intent(inout) :: s
        type(t_partit), intent(in)    :: partit
        integer :: ierr
        call MPI_Allreduce(MPI_IN_PLACE, s, 1, MPI_WP, MPI_MAX, &
                           partit%MPI_COMM_FESOM, ierr)
    end subroutine allreduce_max_r0

    subroutine allreduce_max_i0(s, partit)
        integer,        intent(inout) :: s
        type(t_partit), intent(in)    :: partit
        integer :: ierr
        call MPI_Allreduce(MPI_IN_PLACE, s, 1, MPI_INTEGER, MPI_MAX, &
                           partit%MPI_COMM_FESOM, ierr)
    end subroutine allreduce_max_i0

    ! ===================== cores =============================
    subroutine core_2D_r(arr, com, comm)
        real(kind=WP),    intent(inout) :: arr(:)
        type(com_struct), intent(in)    :: com
        integer,          intent(in)    :: comm
        real(kind=WP), allocatable :: sbuf(:), rbuf(:)
        integer, allocatable :: req(:)
        integer :: i, off, n, nreq, ierr
        if (com%rPEnum == 0 .and. com%sPEnum == 0) return
        allocate(req(com%rPEnum + com%sPEnum))
        allocate(rbuf(max(com%rptr(com%rPEnum+1)-1, 0)))
        allocate(sbuf(max(com%sptr(com%sPEnum+1)-1, 0)))
        nreq = 0
        do i = 1, com%rPEnum
            off = com%rptr(i); n = com%rptr(i+1) - com%rptr(i)
            nreq = nreq + 1
            call MPI_Irecv(rbuf(off), n, MPI_WP, com%rPE(i), HALO_TAG, comm, req(nreq), ierr)
        end do
        do i = 1, com%sPEnum
            off = com%sptr(i); n = com%sptr(i+1) - com%sptr(i)
            sbuf(off:off+n-1) = arr(com%slist(off:off+n-1))
            nreq = nreq + 1
            call MPI_Isend(sbuf(off), n, MPI_WP, com%sPE(i), HALO_TAG, comm, req(nreq), ierr)
        end do
        call MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE, ierr)
        do i = 1, com%rPEnum
            off = com%rptr(i); n = com%rptr(i+1) - com%rptr(i)
            arr(com%rlist(off:off+n-1)) = rbuf(off:off+n-1)   ! broadcast-only
        end do
    end subroutine core_2D_r

    subroutine core_2D_i(arr, com, comm)
        integer,          intent(inout) :: arr(:)
        type(com_struct), intent(in)    :: com
        integer,          intent(in)    :: comm
        integer, allocatable :: sbuf(:), rbuf(:), req(:)
        integer :: i, off, n, nreq, ierr
        if (com%rPEnum == 0 .and. com%sPEnum == 0) return
        allocate(req(com%rPEnum + com%sPEnum))
        allocate(rbuf(max(com%rptr(com%rPEnum+1)-1, 0)))
        allocate(sbuf(max(com%sptr(com%sPEnum+1)-1, 0)))
        nreq = 0
        do i = 1, com%rPEnum
            off = com%rptr(i); n = com%rptr(i+1) - com%rptr(i)
            nreq = nreq + 1
            call MPI_Irecv(rbuf(off), n, MPI_INTEGER, com%rPE(i), HALO_TAG, comm, req(nreq), ierr)
        end do
        do i = 1, com%sPEnum
            off = com%sptr(i); n = com%sptr(i+1) - com%sptr(i)
            sbuf(off:off+n-1) = arr(com%slist(off:off+n-1))
            nreq = nreq + 1
            call MPI_Isend(sbuf(off), n, MPI_INTEGER, com%sPE(i), HALO_TAG, comm, req(nreq), ierr)
        end do
        call MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE, ierr)
        do i = 1, com%rPEnum
            off = com%rptr(i); n = com%rptr(i+1) - com%rptr(i)
            arr(com%rlist(off:off+n-1)) = rbuf(off:off+n-1)
        end do
    end subroutine core_2D_i

    subroutine core_3D_r(arr, com, comm)
        real(kind=WP),    intent(inout) :: arr(:,:)   ! (nl, n)
        type(com_struct), intent(in)    :: com
        integer,          intent(in)    :: comm
        real(kind=WP), allocatable :: sbuf(:,:), rbuf(:,:)
        integer, allocatable :: req(:)
        integer :: i, off, n, nreq, ierr, nl, k
        if (com%rPEnum == 0 .and. com%sPEnum == 0) return
        nl = size(arr, 1)
        allocate(req(com%rPEnum + com%sPEnum))
        allocate(rbuf(nl, max(com%rptr(com%rPEnum+1)-1, 0)))
        allocate(sbuf(nl, max(com%sptr(com%sPEnum+1)-1, 0)))
        nreq = 0
        do i = 1, com%rPEnum
            off = com%rptr(i); n = com%rptr(i+1) - com%rptr(i)
            nreq = nreq + 1
            call MPI_Irecv(rbuf(1,off), nl*n, MPI_WP, com%rPE(i), HALO_TAG, comm, req(nreq), ierr)
        end do
        do i = 1, com%sPEnum
            off = com%sptr(i); n = com%sptr(i+1) - com%sptr(i)
            do k = 0, n-1
                sbuf(:, off+k) = arr(:, com%slist(off+k))
            end do
            nreq = nreq + 1
            call MPI_Isend(sbuf(1,off), nl*n, MPI_WP, com%sPE(i), HALO_TAG, comm, req(nreq), ierr)
        end do
        call MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE, ierr)
        do i = 1, com%rPEnum
            off = com%rptr(i); n = com%rptr(i+1) - com%rptr(i)
            do k = 0, n-1
                arr(:, com%rlist(off+k)) = rbuf(:, off+k)   ! broadcast-only
            end do
        end do
    end subroutine core_3D_r

    subroutine core_blk_r(arr, com, comm)
        ! Broadcast-only exchange of a (d1, d2, entity) array — the leading two dims
        ! form a contiguous block (blk=d1*d2) transported per halo entity (e.g.
        ! tr_xy(2,nl-1,*)). FLAT rank-1 MPI buffers + explicit copy (the proven core_2D_r
        ! pattern): pass a rank-1 element to MPI (clean sequence association), never a
        ! rank-3 array element.
        real(kind=WP),    intent(inout) :: arr(:,:,:)
        type(com_struct), intent(in)    :: com
        integer,          intent(in)    :: comm
        real(kind=WP), allocatable :: sbuf(:), rbuf(:)
        integer, allocatable :: req(:)
        integer :: i, off, n, nreq, ierr, d1, d2, blk, k, base, j1, j2, e
        if (com%rPEnum == 0 .and. com%sPEnum == 0) return
        d1 = size(arr, 1); d2 = size(arr, 2); blk = d1*d2
        allocate(req(com%rPEnum + com%sPEnum))
        allocate(rbuf(blk * max(com%rptr(com%rPEnum+1)-1, 0)))
        allocate(sbuf(blk * max(com%sptr(com%sPEnum+1)-1, 0)))
        nreq = 0
        do i = 1, com%rPEnum
            off = com%rptr(i); n = com%rptr(i+1) - com%rptr(i)
            nreq = nreq + 1
            call MPI_Irecv(rbuf((off-1)*blk+1), blk*n, MPI_WP, com%rPE(i), BLK_TAG, comm, req(nreq), ierr)
        end do
        do i = 1, com%sPEnum
            off = com%sptr(i); n = com%sptr(i+1) - com%sptr(i)
            do k = 0, n-1
                e = com%slist(off+k); base = (off-1+k)*blk
                do j2 = 1, d2
                    do j1 = 1, d1
                        sbuf(base + (j2-1)*d1 + j1) = arr(j1, j2, e)
                    end do
                end do
            end do
            nreq = nreq + 1
            call MPI_Isend(sbuf((off-1)*blk+1), blk*n, MPI_WP, com%sPE(i), BLK_TAG, comm, req(nreq), ierr)
        end do
        call MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE, ierr)
        do i = 1, com%rPEnum
            off = com%rptr(i); n = com%rptr(i+1) - com%rptr(i)
            do k = 0, n-1
                e = com%rlist(off+k); base = (off-1+k)*blk
                do j2 = 1, d2
                    do j1 = 1, d1
                        arr(j1, j2, e) = rbuf(base + (j2-1)*d1 + j1)
                    end do
                end do
            end do
        end do
    end subroutine core_blk_r

    ! ============ stale-halo probe (exchange-and-compare; used by test_halo) ============
    real(kind=WP) function stale_halo_max_nod(arr, partit) result(maxdiff)
        real(kind=WP),  intent(in) :: arr(:)
        type(t_partit), intent(in) :: partit
        maxdiff = probe(arr, partit%com_nod2D, partit%MPI_COMM_FESOM)
    end function

    real(kind=WP) function stale_halo_max_elem(arr, partit) result(maxdiff)
        real(kind=WP),  intent(in) :: arr(:)
        type(t_partit), intent(in) :: partit
        maxdiff = probe(arr, partit%com_elem2D, partit%MPI_COMM_FESOM)
    end function

    ! ============ shared-element copy spread (step diagnostic, test_shared_elem) ============
    real(kind=WP) function shared_copy_max_elem_2D_r(arr, partit) result(maxdiff)
        real(kind=WP),  intent(in) :: arr(:)
        type(t_partit), intent(in) :: partit
        maxdiff = probe(arr, partit%com_elem2D_shared, partit%MPI_COMM_FESOM)
    end function

    real(kind=WP) function shared_copy_max_elem_blk_r(arr, partit) result(maxdiff)
        real(kind=WP),  intent(in) :: arr(:,:,:)
        type(t_partit), intent(in) :: partit
        maxdiff = probe_blk(arr, partit%com_elem2D_shared, partit%MPI_COMM_FESOM)
    end function

    real(kind=WP) function probe_blk(arr, com, comm) result(maxdiff)
        ! probe on a (d1, d2, entity) array through core_blk_r
        real(kind=WP),    intent(in) :: arr(:,:,:)
        type(com_struct), intent(in) :: com
        integer,          intent(in) :: comm
        real(kind=WP), allocatable :: cpy(:,:,:)
        integer :: i, k, e
        maxdiff = 0.0_WP
        if (com%rPEnum == 0 .and. com%sPEnum == 0) return
        cpy = arr
        call core_blk_r(cpy, com, comm)
        do i = 1, com%rPEnum
            do k = com%rptr(i), com%rptr(i+1) - 1
                e = com%rlist(k)
                maxdiff = max(maxdiff, maxval(abs(arr(:,:,e) - cpy(:,:,e))))
            end do
        end do
    end function probe_blk

    real(kind=WP) function probe(arr, com, comm) result(maxdiff)
        ! Exchange a copy and report max|arr(halo) - exchanged(halo)|. >0 => stale halo.
        real(kind=WP),    intent(in) :: arr(:)
        type(com_struct), intent(in) :: com
        integer,          intent(in) :: comm
        real(kind=WP), allocatable :: cpy(:)
        integer :: i, off, n
        maxdiff = 0.0_WP
        cpy = arr
        call core_2D_r(cpy, com, comm)
        do i = 1, com%rPEnum
            off = com%rptr(i); n = com%rptr(i+1) - com%rptr(i)
            maxdiff = max(maxdiff, maxval(abs(arr(com%rlist(off:off+n-1)) - cpy(com%rlist(off:off+n-1)))))
        end do
    end function probe

end module mod_halo
