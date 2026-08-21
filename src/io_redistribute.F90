!> Redistribution schedule for parallel output and restart writing.
!>
!> FESOM's model decomposition gives each rank a sorted but non-contiguous set of
!> global indices: on NG5/dist_1024 a rank owns ~7229 nodes in ~130 runs of ~58,
!> spread over most of the global index range. A rank therefore cannot write its
!> data as one hyperslab, and because run counts differ between ranks, a per-run
!> nf90_put_var loop cannot be collective -- which is what parallel netCDF needs.
!>
!> So we redistribute first. Writer w owns the contiguous global range
!> [w*block+1, (w+1)*block], one MPI_Alltoallv moves the data there, and every
!> writer then issues exactly one put_var per level with identical call counts.
!> The block is UNIFORM (ceil(n_global/n_writers)) rather than balanced: netCDF
!> takes one chunk shape per variable, so only a uniform block lets the chunk
!> size along the node dimension equal the writer block. That equality is what
!> makes the later compressed-parallel case free of chunk redistribution inside
!> HDF5 -- each chunk lies entirely within one writer's range.
!>
!> This module is deliberately free of MOD_PARTIT: it takes plain arrays and an
!> MPI communicator, so it can be unit-tested in test/fortran_parallel, which
!> compiles no partit chain.
!>
!> Index lists passed in must be ASCENDING and must cover only owned entries --
!> for nodes that is myList_nod2D(1:myDim_nod2D), never the halo tail. For
!> elements it is myList_elem2D(myInd_elem2D_shrinked(:)), NOT myList_elem2D
!> itself: an element belongs to a PE if ANY of its nodes does (gen_comm.F90:265),
!> so sum(myDim_elem2D) > elem2D and there is no bijection to build a schedule on.
!> The shrinked set assigns each element to the rank owning its FIRST node
!> (oce_mesh.F90:888-898), which is the authoritative owner for this module.
!>
!> ONE MORE RULE FOR CALLERS THAT BUILD BOTH A NODAL AND AN ELEMENT SCHEDULE:
!> build the nodal one first and pass ITS n_writers as the request for the
!> element one. Do not hand the same namelist value to both. The requested count
!> is clamped to n_global/REDIST_MIN_BLOCK, and elem2D is about twice nod2D, so
!> one request can be halved for nodes and left alone for elements. The two
!> schedules then choose writers with different strides and no longer sit on the
!> same ranks -- while the caller has only one file, opened over one of the two
!> sets. Everything the other set holds is lost without an error: on core2 with
!> n_writers > 30 that is half of every element field (#969).

module io_redistribute
  use mpi
  implicit none
  private

  public :: redist_type
  public :: redist_block_size, redist_owner_of, redist_writer_rank
  public :: redist_effective_writers
  public :: redist_build, redist_free, redist_check
  public :: redist_exchange
  public :: redist_scatter_type
  public :: redist_build_scatter, redist_scatter, redist_scatter_free
  public :: REDIST_MIN_BLOCK

  !> Move a field from the model decomposition into block order.
  !>
  !> One call moves a WHOLE 3D field, not one level: FESOM stores 3D arrays as
  !> (nz, node), so a node's column is contiguous and the schedule's counts need
  !> only be scaled by nlev. That turns 69 level-wise exchanges into one, which
  !> is where most of the latency saving is. 2D fields are simply nlev = 1.
  interface redist_exchange
     module procedure redist_exchange_r4
     module procedure redist_exchange_r8
  end interface redist_exchange

  !> Smallest number of entries per writer block. Guards against pathological
  !> chunk sizes: core2 (nod2D = 126858) with 1024 writers would otherwise give
  !> 124 nodes, about 496 bytes per chunk, which explodes chunk-metadata volume
  !> and destroys compression -- and core2 is the correctness proxy, so that
  !> configuration will be exercised. Requesting more writers than this allows
  !> silently reduces the writer count; redist_build reports when it does.
  integer, parameter :: REDIST_MIN_BLOCK = 4096

  type :: redist_type
     integer :: n_global  = 0     !< total entries across all ranks
     integer :: n_writers = 0     !< effective writer count (after clamping)
     integer :: block     = 0     !< entries per writer, uniform
     integer :: npes      = 0
     integer :: mype      = -1
     logical :: is_writer = .false.
     integer :: widx      = -1    !< writer index, 0..n_writers-1, else -1
     integer :: first     = 0     !< first global index this writer owns (1-based)
     integer :: last      = -1    !< last global index this writer owns
     integer :: count     = 0     !< last - first + 1, may be 0
     integer :: nsend     = 0
     integer :: nrecv     = 0
     integer, allocatable :: scnt(:), sdsp(:)   !< 0:npes-1
     integer, allocatable :: rcnt(:), rdsp(:)   !< 0:npes-1
     integer, allocatable :: rloc(:)            !< 1:nrecv, offset within the block
  end type redist_type

  !> Schedule for the opposite direction: block order on the readers out to
  !> every rank that needs an entry.
  !>
  !> This is deliberately NOT redist_type reversed. A write schedule is a
  !> partition -- redist_build fails with ierr = 3 if a writer would not receive
  !> exactly its block -- because each global entry is written once. A read is
  !> not a partition: a rank needs every entry it holds locally, owned or halo,
  !> and neighbouring ranks need the same entry at the same time.
  !>
  !> For nodes the difference is only the halo, which exchange_nod3D could patch
  !> up afterwards. For ELEMENTS there is a second gap that no halo exchange
  !> closes, and it is the reason this type exists. myInd_elem2D_shrinked keeps
  !> only local elements whose FIRST node is owned (oce_mesh.F90:888-899), while
  !> an element enters the halo receive list only when NONE of its nodes belongs
  !> to this PE (gen_comm.F90:430). An element in 1:myDim_elem2D whose first node
  !> is a halo node is in neither set: inverting the write schedule would leave
  !> it holding whatever was in memory, and this rank would then pass that on,
  !> because the element does appear in its own send lists. Nothing catches this
  !> during normal running, where such elements are recomputed from nodal values
  !> every step; a restart read is the one place the value has to come from
  !> outside.
  !>
  !> So each rank simply asks for what it needs and readers serve the requests.
  !> Duplicate requests across ranks are expected, not an error.
  type :: redist_scatter_type
     integer :: n_global  = 0
     integer :: n_readers = 0     !< same blocks as the writers, so the same count
     integer :: block     = 0
     integer :: npes      = 0
     integer :: mype      = -1
     logical :: is_reader = .false.
     integer :: widx      = -1
     integer :: first     = 0     !< first global index this reader owns (1-based)
     integer :: last      = -1
     integer :: count     = 0     !< block length held by this reader, may be 0
     integer :: nwant     = 0     !< entries this rank needs (owned + halo)
     integer :: nserve    = 0     !< entries this rank must send out as a reader
     integer, allocatable :: qcnt(:), qdsp(:)   !< 0:npes-1, what we ask of each reader
     integer, allocatable :: scnt(:), sdsp(:)   !< 0:npes-1, what we serve each requester
     integer, allocatable :: sloc(:)            !< 1:nserve, offset within our block
     integer, allocatable :: qperm(:)           !< 1:nwant, request slot -> local position
  end type redist_scatter_type

contains

  !> Uniform block size. ceil(n_global / n_writers).
  pure function redist_block_size(n_global, n_writers) result(b)
    integer, intent(in) :: n_global, n_writers
    integer :: b
    if (n_writers <= 0) then
       b = n_global
    else
       b = (n_global + n_writers - 1) / n_writers
    end if
  end function redist_block_size

  !> Writer count actually used, after clamping to the rank count and to
  !> REDIST_MIN_BLOCK. Always >= 1.
  !>
  !> This clamp is why a node/element schedule pair must share one count: it
  !> depends on n_global, which differs between the two. See the module header.
  pure function redist_effective_writers(n_global, n_writers_req, npes) result(nw)
    integer, intent(in) :: n_global, n_writers_req, npes
    integer :: nw
    nw = n_writers_req
    if (nw <= 0)    nw = npes           ! 0 or negative means "all ranks"
    if (nw > npes)  nw = npes
    nw = min(nw, max(1, n_global / REDIST_MIN_BLOCK))
    if (nw < 1) nw = 1
  end function redist_effective_writers

  !> Writer index owning global index g (1-based g, 0-based result).
  pure function redist_owner_of(g, block, n_writers) result(w)
    integer, intent(in) :: g, block, n_writers
    integer :: w
    w = (g - 1) / block
    if (w > n_writers - 1) w = n_writers - 1
    if (w < 0) w = 0
  end function redist_owner_of

  !> MPI rank hosting writer w. Writers are spread by stride rather than packed
  !> onto ranks 0..n_writers-1, so they do not all land on the first few nodes.
  pure function redist_writer_rank(w, npes, n_writers) result(r)
    integer, intent(in) :: w, npes, n_writers
    integer :: r, stride
    stride = npes / max(1, n_writers)
    if (stride < 1) stride = 1
    r = w * stride
    if (r > npes - 1) r = npes - 1
  end function redist_writer_rank

  !> Build the schedule. glist(1:nloc) must be ascending global indices, owned
  !> entries only. Collective over comm.
  subroutine redist_build(this, glist, nloc, n_global, n_writers_req, comm, verbose, ierr)
    type(redist_type), intent(out) :: this
    integer, intent(in)  :: nloc, n_global, n_writers_req, comm
    integer, intent(in)  :: glist(nloc)
    logical, intent(in)  :: verbose
    integer, intent(out) :: ierr

    integer :: i, w, r, mpierr, prev
    integer, allocatable :: rglob(:)

    ierr = 0
    call MPI_Comm_size(comm, this%npes, mpierr)
    call MPI_Comm_rank(comm, this%mype, mpierr)

    this%n_global  = n_global
    this%n_writers = redist_effective_writers(n_global, n_writers_req, this%npes)
    this%block     = redist_block_size(n_global, this%n_writers)

    if (verbose .and. this%mype == 0 .and. this%n_writers /= min(max(n_writers_req,1), this%npes)) then
       write(*,'(a,i0,a,i0,a,i0,a)') &
          ' io_redistribute: writer count reduced from ', n_writers_req, ' to ', &
          this%n_writers, ' so each block holds at least ', REDIST_MIN_BLOCK, ' entries'
    end if

    ! --- is this rank a writer, and which range does it own -----------------
    this%is_writer = .false.
    this%widx      = -1
    do w = 0, this%n_writers - 1
       if (redist_writer_rank(w, this%npes, this%n_writers) == this%mype) then
          this%is_writer = .true.
          this%widx      = w
          exit
       end if
    end do
    if (this%is_writer) then
       this%first = this%widx * this%block + 1
       this%last  = min((this%widx + 1) * this%block, n_global)
       this%count = max(0, this%last - this%first + 1)
    else
       this%first = 0; this%last = -1; this%count = 0
    end if

    ! --- send plan -----------------------------------------------------------
    ! glist is ascending and redist_owner_of is monotone in g, so entries are
    ! already grouped by destination: only the run boundaries need counting and
    ! no sort or reorder buffer is required. Verified, not assumed.
    allocate(this%scnt(0:this%npes-1), this%sdsp(0:this%npes-1))
    allocate(this%rcnt(0:this%npes-1), this%rdsp(0:this%npes-1))
    this%scnt = 0

    prev = 0
    do i = 1, nloc
       if (glist(i) < 1 .or. glist(i) > n_global) then
          ierr = 1
          if (verbose) write(*,'(a,i0,a,i0)') ' io_redistribute: index out of range: ', glist(i), &
               ' on rank ', this%mype
          return
       end if
       if (glist(i) < prev) then
          ierr = 2                       ! not ascending -- the grouping assumption fails
          if (verbose) write(*,'(a,i0)') ' io_redistribute: index list not ascending on rank ', this%mype
          return
       end if
       prev = glist(i)
       w = redist_owner_of(glist(i), this%block, this%n_writers)
       r = redist_writer_rank(w, this%npes, this%n_writers)
       this%scnt(r) = this%scnt(r) + 1
    end do

    this%sdsp(0) = 0
    do r = 1, this%npes - 1
       this%sdsp(r) = this%sdsp(r-1) + this%scnt(r-1)
    end do
    this%nsend = sum(this%scnt)

    ! --- receive plan --------------------------------------------------------
    call MPI_Alltoall(this%scnt, 1, MPI_INTEGER, this%rcnt, 1, MPI_INTEGER, comm, mpierr)
    this%rdsp(0) = 0
    do r = 1, this%npes - 1
       this%rdsp(r) = this%rdsp(r-1) + this%rcnt(r-1)
    end do
    this%nrecv = sum(this%rcnt)

    if (this%nrecv /= this%count) then
       ierr = 3                          ! a writer must receive exactly its block
       if (verbose) write(*,'(a,i0,a,i0,a,i0)') ' io_redistribute: rank ', this%mype, &
            ' expects ', this%count, ' entries but will receive ', this%nrecv
       return
    end if

    ! --- where does each incoming entry land inside the block ---------------
    allocate(rglob(max(1, this%nrecv)))
    call MPI_Alltoallv(glist, this%scnt, this%sdsp, MPI_INTEGER, &
                       rglob, this%rcnt, this%rdsp, MPI_INTEGER, comm, mpierr)

    allocate(this%rloc(max(1, this%nrecv)))
    do i = 1, this%nrecv
       this%rloc(i) = rglob(i) - this%first + 1
       if (this%rloc(i) < 1 .or. this%rloc(i) > this%count) then
          ierr = 4
          if (verbose) write(*,'(a,i0,a,i0)') ' io_redistribute: received index ', rglob(i), &
               ' outside block on rank ', this%mype
          deallocate(rglob)
          return
       end if
    end do
    deallocate(rglob)
  end subroutine redist_build

  !> sbuf(nlev, nsend) in local ascending order -> out(nlev, count) in block
  !> order on writers. Non-writers pass a size-0 out and receive nothing.
  subroutine redist_exchange_r4(this, nlev, sbuf, out, comm, ierr)
    type(redist_type), intent(in) :: this
    integer, intent(in)  :: nlev, comm
    real(kind=4), intent(in)  :: sbuf(nlev, *)
    real(kind=4), intent(out) :: out(nlev, *)
    integer, intent(out) :: ierr
    real(kind=4), allocatable :: rbuf(:,:)
    integer, allocatable :: sc(:), sd(:), rc(:), rd(:)
    integer :: i, mpierr

    ierr = 0
    allocate(sc(0:this%npes-1), sd(0:this%npes-1), rc(0:this%npes-1), rd(0:this%npes-1))
    sc = this%scnt * nlev ; sd = this%sdsp * nlev
    rc = this%rcnt * nlev ; rd = this%rdsp * nlev
    allocate(rbuf(nlev, max(1, this%nrecv)))

    call MPI_Alltoallv(sbuf, sc, sd, MPI_REAL4, rbuf, rc, rd, MPI_REAL4, comm, mpierr)
    if (mpierr /= MPI_SUCCESS) ierr = mpierr

    ! received data is ordered by source rank; rloc puts it in block order
    do i = 1, this%nrecv
       out(1:nlev, this%rloc(i)) = rbuf(1:nlev, i)
    end do

    deallocate(rbuf, sc, sd, rc, rd)
  end subroutine redist_exchange_r4

  subroutine redist_exchange_r8(this, nlev, sbuf, out, comm, ierr)
    type(redist_type), intent(in) :: this
    integer, intent(in)  :: nlev, comm
    real(kind=8), intent(in)  :: sbuf(nlev, *)
    real(kind=8), intent(out) :: out(nlev, *)
    integer, intent(out) :: ierr
    real(kind=8), allocatable :: rbuf(:,:)
    integer, allocatable :: sc(:), sd(:), rc(:), rd(:)
    integer :: i, mpierr

    ierr = 0
    allocate(sc(0:this%npes-1), sd(0:this%npes-1), rc(0:this%npes-1), rd(0:this%npes-1))
    sc = this%scnt * nlev ; sd = this%sdsp * nlev
    rc = this%rcnt * nlev ; rd = this%rdsp * nlev
    allocate(rbuf(nlev, max(1, this%nrecv)))

    call MPI_Alltoallv(sbuf, sc, sd, MPI_REAL8, rbuf, rc, rd, MPI_REAL8, comm, mpierr)
    if (mpierr /= MPI_SUCCESS) ierr = mpierr

    do i = 1, this%nrecv
       out(1:nlev, this%rloc(i)) = rbuf(1:nlev, i)
    end do

    deallocate(rbuf, sc, sd, rc, rd)
  end subroutine redist_exchange_r8



  !> Build the read-side schedule. `want(1:nwant)` is every global index this
  !> rank needs a value for -- owned AND halo, in whatever order the model keeps
  !> them. Collective over comm.
  !>
  !> `n_readers_req` must be the value the WRITER used, so that first/last/count
  !> describe the same blocks the file was written in. Passing something else is
  !> not an error here, but then the reader blocks will not line up with the
  !> file's chunks and every read crosses chunk boundaries.
  !>
  !> Unlike redist_build this does NOT require `want` to be ascending. The write
  !> side gets away with the grouping shortcut because myList_nod2D(1:myDim) is
  !> sorted; append the halo tail and it is not, so the destination-sorted order
  !> is built explicitly and qperm remembers the way back.
  subroutine redist_build_scatter(this, want, nwant, n_global, n_readers_req, comm, verbose, ierr)
    type(redist_scatter_type), intent(out) :: this
    integer, intent(in)  :: nwant, n_global, n_readers_req, comm
    integer, intent(in)  :: want(nwant)
    logical, intent(in)  :: verbose
    integer, intent(out) :: ierr

    integer :: i, w, r, mpierr, slot
    integer, allocatable :: fill(:), qglob(:), sglob(:)

    ierr = 0
    call MPI_Comm_size(comm, this%npes, mpierr)
    call MPI_Comm_rank(comm, this%mype, mpierr)

    this%n_global  = n_global
    this%n_readers = redist_effective_writers(n_global, n_readers_req, this%npes)
    this%block     = redist_block_size(n_global, this%n_readers)
    this%nwant     = nwant

    ! Reader set and block range: identical rules to the writer side, so the
    ! two agree by construction rather than by coincidence.
    this%is_reader = .false.
    this%widx      = -1
    do w = 0, this%n_readers - 1
       if (redist_writer_rank(w, this%npes, this%n_readers) == this%mype) then
          this%is_reader = .true.
          this%widx      = w
          exit
       end if
    end do
    if (this%is_reader) then
       this%first = this%widx * this%block + 1
       this%last  = min((this%widx + 1) * this%block, n_global)
       this%count = max(0, this%last - this%first + 1)
    else
       this%first = 0; this%last = -1; this%count = 0
    end if

    ! --- what we ask of each reader -----------------------------------------
    allocate(this%qcnt(0:this%npes-1), this%qdsp(0:this%npes-1))
    allocate(this%scnt(0:this%npes-1), this%sdsp(0:this%npes-1))
    this%qcnt = 0

    do i = 1, nwant
       if (want(i) < 1 .or. want(i) > n_global) then
          ierr = 1
          if (verbose) write(*,'(a,i0,a,i0)') ' io_redistribute: requested index out of range: ', &
               want(i), ' on rank ', this%mype
          return
       end if
       w = redist_owner_of(want(i), this%block, this%n_readers)
       r = redist_writer_rank(w, this%npes, this%n_readers)
       this%qcnt(r) = this%qcnt(r) + 1
    end do

    this%qdsp(0) = 0
    do r = 1, this%npes - 1
       this%qdsp(r) = this%qdsp(r-1) + this%qcnt(r-1)
    end do

    ! Destination-sorted request list, plus the permutation back to local slots.
    allocate(qglob(max(1, nwant)), this%qperm(max(1, nwant)), fill(0:this%npes-1))
    fill = this%qdsp
    do i = 1, nwant
       w = redist_owner_of(want(i), this%block, this%n_readers)
       r = redist_writer_rank(w, this%npes, this%n_readers)
       fill(r) = fill(r) + 1
       slot = fill(r)
       qglob(slot)      = want(i)
       this%qperm(slot) = i
    end do
    deallocate(fill)

    ! --- what we must serve --------------------------------------------------
    call MPI_Alltoall(this%qcnt, 1, MPI_INTEGER, this%scnt, 1, MPI_INTEGER, comm, mpierr)
    this%sdsp(0) = 0
    do r = 1, this%npes - 1
       this%sdsp(r) = this%sdsp(r-1) + this%scnt(r-1)
    end do
    this%nserve = sum(this%scnt)

    if (this%nserve > 0 .and. .not. this%is_reader) then
       ierr = 2                       ! a non-reader was asked for data
       if (verbose) write(*,'(a,i0,a,i0)') ' io_redistribute: rank ', this%mype, &
            ' is not a reader but was asked for ', this%nserve
       deallocate(qglob)
       return
    end if

    allocate(sglob(max(1, this%nserve)))
    call MPI_Alltoallv(qglob, this%qcnt, this%qdsp, MPI_INTEGER, &
                       sglob, this%scnt, this%sdsp, MPI_INTEGER, comm, mpierr)
    deallocate(qglob)

    allocate(this%sloc(max(1, this%nserve)))
    do i = 1, this%nserve
       this%sloc(i) = sglob(i) - this%first + 1
       if (this%sloc(i) < 1 .or. this%sloc(i) > this%count) then
          ierr = 3
          if (verbose) write(*,'(a,i0,a,i0)') ' io_redistribute: asked for index ', sglob(i), &
               ' outside block on rank ', this%mype
          deallocate(sglob)
          return
       end if
    end do
    deallocate(sglob)
  end subroutine redist_build_scatter


  !> blk(nlev, count) in block order on the readers -> out(nlev, nwant) in the
  !> model's own local order on every rank. Non-readers pass a size-0 blk.
  !>
  !> Mirror image of redist_exchange: gather into the send buffer through sloc,
  !> one Alltoallv, then scatter through qperm.
  subroutine redist_scatter(this, nlev, blk, out, comm, ierr)
    type(redist_scatter_type), intent(in) :: this
    integer, intent(in)  :: nlev, comm
    real(kind=8), intent(in)  :: blk(nlev, *)
    real(kind=8), intent(out) :: out(nlev, *)
    integer, intent(out) :: ierr
    real(kind=8), allocatable :: sbuf(:,:), rbuf(:,:)
    integer, allocatable :: sc(:), sd(:), rc(:), rd(:)
    integer :: i, mpierr

    ierr = 0
    allocate(sc(0:this%npes-1), sd(0:this%npes-1), rc(0:this%npes-1), rd(0:this%npes-1))
    sc = this%scnt * nlev ; sd = this%sdsp * nlev
    rc = this%qcnt * nlev ; rd = this%qdsp * nlev

    allocate(sbuf(nlev, max(1, this%nserve)))
    do i = 1, this%nserve
       sbuf(1:nlev, i) = blk(1:nlev, this%sloc(i))
    end do

    allocate(rbuf(nlev, max(1, this%nwant)))
    call MPI_Alltoallv(sbuf, sc, sd, MPI_REAL8, rbuf, rc, rd, MPI_REAL8, comm, mpierr)
    if (mpierr /= MPI_SUCCESS) ierr = mpierr

    do i = 1, this%nwant
       out(1:nlev, this%qperm(i)) = rbuf(1:nlev, i)
    end do

    deallocate(sbuf, rbuf, sc, sd, rc, rd)
  end subroutine redist_scatter


  subroutine redist_scatter_free(this)
    type(redist_scatter_type), intent(inout) :: this
    if (allocated(this%qcnt))  deallocate(this%qcnt)
    if (allocated(this%qdsp))  deallocate(this%qdsp)
    if (allocated(this%scnt))  deallocate(this%scnt)
    if (allocated(this%sdsp))  deallocate(this%sdsp)
    if (allocated(this%sloc))  deallocate(this%sloc)
    if (allocated(this%qperm)) deallocate(this%qperm)
    this%is_reader = .false.
    this%widx   = -1
    this%count  = 0
    this%nwant  = 0
    this%nserve = 0
  end subroutine redist_scatter_free


  !> Confirm the schedule partitions the global index space exactly once.
  !> Collective. Allocates two n_global integer arrays, so this is a test and
  !> debug routine, not something to call every write at NG5 scale.
  function redist_check(this, comm) result(ok)
    type(redist_type), intent(in) :: this
    integer, intent(in) :: comm
    logical :: ok
    integer :: mpierr, total_recv, total_send
    integer, allocatable :: hit(:), hit_all(:)
    integer :: i

    call MPI_Allreduce(this%nrecv, total_recv, 1, MPI_INTEGER, MPI_SUM, comm, mpierr)
    call MPI_Allreduce(this%nsend, total_send, 1, MPI_INTEGER, MPI_SUM, comm, mpierr)
    ok = (total_recv == this%n_global) .and. (total_send == this%n_global)

    ! every global index claimed by exactly one writer
    allocate(hit(this%n_global), hit_all(this%n_global))
    hit = 0
    do i = 1, this%nrecv
       hit(this%first + this%rloc(i) - 1) = 1
    end do
    call MPI_Allreduce(hit, hit_all, this%n_global, MPI_INTEGER, MPI_SUM, comm, mpierr)
    ok = ok .and. all(hit_all == 1)
    deallocate(hit, hit_all)
  end function redist_check

  subroutine redist_free(this)
    type(redist_type), intent(inout) :: this
    if (allocated(this%scnt)) deallocate(this%scnt)
    if (allocated(this%sdsp)) deallocate(this%sdsp)
    if (allocated(this%rcnt)) deallocate(this%rcnt)
    if (allocated(this%rdsp)) deallocate(this%rdsp)
    if (allocated(this%rloc)) deallocate(this%rloc)
    this%is_writer = .false.
    this%widx = -1
    this%count = 0
  end subroutine redist_free

end module io_redistribute
