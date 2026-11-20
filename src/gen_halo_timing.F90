! Time spent in halo exchanges and in the global sums of the SSH solver, over the time loop
! (halo_timing in namelist.config). Each exchange is split into posting (from the call of an
! exchange_*_begin to its return: datatype setup, device packing, MPI_Isend/Irecv) and completing
! (exchange_nod_end / exchange_elem_end: MPI_Waitall and device unpacking), by kind: nodal or
! elemental, 2D (one value per point, or up to four values) or 3D (levels). Report: slowest rank
! and mean over ranks, written by rank 0 in fesom_finalize. Rank 0 also writes halo_ranks.txt: one
! line per rank with its work (from the caller), its post and wait time per category and phase (sea
! ice, ocean, other; set with ht_phase) and the mean, deviation and maximum over time steps of the
! step time outside exchanges and global sums (ht_step).
module halo_timing_mod
  use mpi
  use, intrinsic :: iso_fortran_env, only: real64, int64
  implicit none
  private
  public :: ht_on, ht_start, ht_report, ht_post, ht_done, ht_kind, ht_sum_mark, ht_sum_done, ht_phase, ht_step
  public :: HT_NOD2D, HT_NOD3D, HT_ELEM2D, HT_ELEM3D

  integer, parameter :: HT_NOD2D=1, HT_NOD3D=2, HT_ELEM2D=3, HT_ELEM3D=4, HT_SUM=5, NCAT=5
  character(len=26), parameter :: label(NCAT) = [character(len=26) :: &
       'nodal 2D exchanges', 'nodal 3D exchanges', 'element 2D exchanges', 'element 3D exchanges', &
       'SSH solver global sums']
  logical,        save :: ht_on = .false.
  logical,        save :: ht_sites = .false.                 ! per-call-site records (gen_halo_sites.c)
  integer,        save :: nl_ref = 0
  integer,        save :: pending(2) = 0                     ! category in flight: 1 nodal, 2 element
  integer,        save :: phase = 0                          ! 0 other, 1 sea ice, 2 ocean
  real(real64),   save :: t_post(NCAT) = 0, t_done(NCAT) = 0, t_mark = 0, t_loop0 = 0
  real(real64),   save :: p_post(NCAT,0:2) = 0, p_done(NCAT,0:2) = 0
  real(real64),   save :: s_t0 = 0, s_halo0 = 0, s_sum = 0, s_sq = 0, s_max = 0, s_ice = 0
  integer(int64), save :: s_n = 0
  integer(int64), save :: n_call(NCAT) = 0, n_field(NCAT) = 0, n_msg(NCAT) = 0

  interface
     subroutine fesom_halo_site(cat, iend, nfields, sec) bind(C, name='fesom_halo_site')
       import :: real64
       integer, value :: cat, iend, nfields
       real(real64), value :: sec
     end subroutine
     subroutine fesom_halo_site_dump(rank) bind(C, name='fesom_halo_site_dump')
       integer, value :: rank
     end subroutine
  end interface

contains

  ! Called at the start of every fesom_runloop; counts from the first call (a coupled model may
  ! call the run loop once per coupling step)
  subroutine ht_start(enable, nl, sites)
    logical, intent(in) :: enable
    integer, intent(in) :: nl
    logical, intent(in), optional :: sites
    if (ht_on .or. .not. enable) return
    ht_on = enable
    if (present(sites)) ht_sites = sites
    nl_ref = nl
    pending = 0
    t_post = 0; t_done = 0; n_call = 0; n_field = 0; n_msg = 0
    p_post = 0; p_done = 0; phase = 0
    s_sum = 0; s_sq = 0; s_max = 0; s_n = 0; s_halo0 = 0; s_ice = 0
    t_loop0 = MPI_Wtime()
    s_t0 = t_loop0
  end subroutine ht_start

  subroutine ht_phase(p)
    integer, intent(in) :: p
    phase = p
  end subroutine ht_phase

  ! End of a time step: time outside exchanges and global sums in this step; nice: owned nodes
  ! with sea ice in this step
  subroutine ht_step(nice)
    integer, intent(in) :: nice
    real(real64) :: now, h, c
    if (.not. ht_on) return
    s_ice = s_ice + nice
    now = MPI_Wtime()
    h = sum(t_post) + sum(t_done)
    c = (now - s_t0) - (h - s_halo0)
    s_sum = s_sum + c; s_sq = s_sq + c*c; s_max = max(s_max, c); s_n = s_n + 1
    s_t0 = now; s_halo0 = h
  end subroutine ht_step

  ! 3D if the level dimension has nl-1 or nl entries, otherwise 2D (one to four values per point)
  integer function ht_kind(side, nlev)
    integer, intent(in) :: side, nlev
    if (side == 1) then
       ht_kind = merge(HT_NOD3D, HT_NOD2D, nlev >= nl_ref-1)
    else
       ht_kind = merge(HT_ELEM3D, HT_ELEM2D, nlev >= nl_ref-1)
    end if
  end function ht_kind

  ! packed: all fields of the exchange travel in one message per neighbour
  subroutine ht_post(side, cat, nfields, nneigh, dt, packed)
    integer,      intent(in) :: side, cat, nfields, nneigh
    real(real64), intent(in) :: dt
    logical,      intent(in) :: packed
    pending(side) = cat
    if (ht_sites) call fesom_halo_site(cat, 0, nfields, dt)
    t_post(cat)  = t_post(cat) + dt
    p_post(cat, phase) = p_post(cat, phase) + dt
    n_call(cat)  = n_call(cat) + 1
    n_field(cat) = n_field(cat) + nfields
    n_msg(cat)   = n_msg(cat) + merge(1, nfields, packed)*int(nneigh, int64)
  end subroutine ht_post

  subroutine ht_done(side, dt)
    integer,      intent(in) :: side
    real(real64), intent(in) :: dt
    if (pending(side) > 0) then
       t_done(pending(side)) = t_done(pending(side)) + dt
       p_done(pending(side), phase) = p_done(pending(side), phase) + dt
    end if
    if (ht_sites .and. pending(side) > 0) call fesom_halo_site(pending(side), 1, 0, dt)
    pending(side) = 0
  end subroutine ht_done

  subroutine ht_sum_mark()
    t_mark = MPI_Wtime()
  end subroutine ht_sum_mark

  subroutine ht_sum_done()
    t_done(HT_SUM) = t_done(HT_SUM) + (MPI_Wtime() - t_mark)
    p_done(HT_SUM, phase) = p_done(HT_SUM, phase) + (MPI_Wtime() - t_mark)
    if (ht_sites) call fesom_halo_site(HT_SUM, 1, 0, MPI_Wtime() - t_mark)
    n_call(HT_SUM) = n_call(HT_SUM) + 1
  end subroutine ht_sum_done

  ! work, work_names: per-rank values from the caller (sizes, section times) for halo_ranks.txt
  subroutine ht_report(comm, mype, npes, work, work_names)
    integer,          intent(in)           :: comm, mype, npes
    real(real64),     intent(in), optional :: work(:)
    character(len=*), intent(in), optional :: work_names(:)
    real(real64)   :: loc(2*NCAT+1), mx(2*NCAT+1), sm(2*NCAT+1)
    integer(int64) :: cnt(3*NCAT), cmx(3*NCAT)
    integer        :: c, ierr
    if (.not. ht_on) return
    if (ht_sites) call fesom_halo_site_dump(mype)
    call ht_rank_file(comm, mype, npes, work, work_names)
    loc(1:NCAT) = t_post; loc(NCAT+1:2*NCAT) = t_done; loc(2*NCAT+1) = MPI_Wtime() - t_loop0
    cnt = [n_call, n_field, n_msg]
    call MPI_Reduce(loc, mx, size(loc), MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)
    call MPI_Reduce(loc, sm, size(loc), MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
    call MPI_Reduce(cnt, cmx, size(cnt), MPI_INTEGER8, MPI_MAX, 0, comm, ierr)
    if (mype /= 0) return
    write(*,'(a)') ' ___HALO TIMING (time loop; seconds; max = slowest rank, mean over ranks)___'
    write(*,'(a,f12.3,a,f12.3)') '   first time step to report: max ', mx(2*NCAT+1), '   mean ', sm(2*NCAT+1)/npes
    write(*,'(a26,3a12,4a12)') 'category', 'calls', 'fields', 'messages', 'post max', 'post mean', &
         'wait max', 'wait mean'
    do c = 1, NCAT
       write(*,'(a26,3i12,4f12.3)') label(c), cmx(c), cmx(NCAT+c), cmx(2*NCAT+c), mx(c), sm(c)/npes, &
            mx(NCAT+c), sm(NCAT+c)/npes
    end do
    write(*,'(a)') '   (messages = fields x receiving neighbours on the busiest rank; for the global sums'
    write(*,'(a)') '    "wait" is the whole MPI_Allreduce)'
  end subroutine ht_report

  subroutine ht_rank_file(comm, mype, npes, work, work_names)
    integer,          intent(in)           :: comm, mype, npes
    real(real64),     intent(in), optional :: work(:)
    character(len=*), intent(in), optional :: work_names(:)
    character(len=*), parameter :: cn(NCAT) = ['n2', 'n3', 'e2', 'e3', 'gs']
    character(len=*), parameter :: pn(0:2) = ['oth', 'ice', 'oce']
    real(real64), allocatable :: row(:), tab(:,:)
    real(real64) :: m
    integer :: nw, nr, c, p, k, r, ierr, u
    nw = 0
    if (present(work)) nw = size(work)
    nr = nw + 5 + 2*NCAT*3
    allocate(row(nr))
    if (nw > 0) row(1:nw) = work
    m = s_sum/max(1_int64, s_n)
    row(nw+1) = MPI_Wtime() - t_loop0
    row(nw+2) = m
    row(nw+3) = sqrt(max(0._real64, s_sq/max(1_int64, s_n) - m*m))
    row(nw+4) = s_max
    row(nw+5) = s_ice/max(1_int64, s_n)
    k = nw + 5
    do p = 0, 2
       do c = 1, NCAT
          row(k+1) = p_post(c, p); row(k+2) = p_done(c, p); k = k + 2
       end do
    end do
    if (mype == 0) then
       allocate(tab(nr, npes))
    else
       allocate(tab(1, 1))
    end if
    call MPI_Gather(row, nr, MPI_DOUBLE_PRECISION, tab, nr, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    if (mype /= 0) return
    open(newunit=u, file='halo_ranks.txt', status='replace', action='write')
    write(u, '(a)', advance='no') 'rank'
    do c = 1, nw
       write(u, '(1x,a)', advance='no') trim(work_names(c))
    end do
    write(u, '(a)', advance='no') ' loop step_mean step_sd step_max ice_nod_mean'
    do p = 0, 2
       do c = 1, NCAT
          write(u, '(5a)', advance='no') ' post_', cn(c), '_', pn(p), ''
          write(u, '(5a)', advance='no') ' wait_', cn(c), '_', pn(p), ''
       end do
    end do
    write(u, *)
    do r = 1, npes
       write(u, '(i6)', advance='no') r - 1
       do c = 1, nr
          write(u, '(1x,es13.6)', advance='no') tab(c, r)
       end do
       write(u, *)
    end do
    close(u)
  end subroutine ht_rank_file
end module halo_timing_mod
