! Time spent in halo exchanges and in the global sums of the SSH solver, over the time loop
! (halo_timing in namelist.config). Each exchange is split into posting (from the call of an
! exchange_*_begin to its return: datatype setup, device packing, MPI_Isend/Irecv) and completing
! (exchange_nod_end / exchange_elem_end: MPI_Waitall and device unpacking), by kind: nodal or
! elemental, 2D (one value per point, or up to four values) or 3D (levels). Report: slowest rank
! and mean over ranks, written by rank 0 in fesom_finalize.
module halo_timing_mod
  use mpi
  use, intrinsic :: iso_fortran_env, only: real64, int64
  implicit none
  private
  public :: ht_on, ht_start, ht_report, ht_post, ht_done, ht_kind, ht_sum_mark, ht_sum_done
  public :: HT_NOD2D, HT_NOD3D, HT_ELEM2D, HT_ELEM3D

  integer, parameter :: HT_NOD2D=1, HT_NOD3D=2, HT_ELEM2D=3, HT_ELEM3D=4, HT_SUM=5, NCAT=5
  character(len=26), parameter :: label(NCAT) = [character(len=26) :: &
       'nodal 2D exchanges', 'nodal 3D exchanges', 'element 2D exchanges', 'element 3D exchanges', &
       'SSH solver global sums']
  logical,        save :: ht_on = .false.
  integer,        save :: nl_ref = 0
  integer,        save :: pending(2) = 0                     ! category in flight: 1 nodal, 2 element
  real(real64),   save :: t_post(NCAT) = 0, t_done(NCAT) = 0, t_mark = 0, t_loop0 = 0
  integer(int64), save :: n_call(NCAT) = 0, n_field(NCAT) = 0, n_msg(NCAT) = 0

contains

  ! Called at the start of every fesom_runloop; counts from the first call (a coupled model may
  ! call the run loop once per coupling step)
  subroutine ht_start(enable, nl)
    logical, intent(in) :: enable
    integer, intent(in) :: nl
    if (ht_on .or. .not. enable) return
    ht_on = enable
    nl_ref = nl
    pending = 0
    t_post = 0; t_done = 0; n_call = 0; n_field = 0; n_msg = 0
    t_loop0 = MPI_Wtime()
  end subroutine ht_start

  ! 3D if the level dimension has nl-1 or nl entries, otherwise 2D (one to four values per point)
  integer function ht_kind(side, nlev)
    integer, intent(in) :: side, nlev
    if (side == 1) then
       ht_kind = merge(HT_NOD3D, HT_NOD2D, nlev >= nl_ref-1)
    else
       ht_kind = merge(HT_ELEM3D, HT_ELEM2D, nlev >= nl_ref-1)
    end if
  end function ht_kind

  subroutine ht_post(side, cat, nfields, nneigh, dt)
    integer,      intent(in) :: side, cat, nfields, nneigh
    real(real64), intent(in) :: dt
    pending(side) = cat
    t_post(cat)  = t_post(cat) + dt
    n_call(cat)  = n_call(cat) + 1
    n_field(cat) = n_field(cat) + nfields
    n_msg(cat)   = n_msg(cat) + int(nfields, int64)*nneigh
  end subroutine ht_post

  subroutine ht_done(side, dt)
    integer,      intent(in) :: side
    real(real64), intent(in) :: dt
    if (pending(side) > 0) t_done(pending(side)) = t_done(pending(side)) + dt
    pending(side) = 0
  end subroutine ht_done

  subroutine ht_sum_mark()
    t_mark = MPI_Wtime()
  end subroutine ht_sum_mark

  subroutine ht_sum_done()
    t_done(HT_SUM) = t_done(HT_SUM) + (MPI_Wtime() - t_mark)
    n_call(HT_SUM) = n_call(HT_SUM) + 1
  end subroutine ht_sum_done

  subroutine ht_report(comm, mype, npes)
    integer, intent(in) :: comm, mype, npes
    real(real64)   :: loc(2*NCAT+1), mx(2*NCAT+1), sm(2*NCAT+1)
    integer(int64) :: cnt(3*NCAT), cmx(3*NCAT)
    integer        :: c, ierr
    if (.not. ht_on) return
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
end module halo_timing_mod
