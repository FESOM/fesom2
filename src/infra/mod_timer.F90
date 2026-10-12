module mod_timer
    !===========================================================================
    ! Permanent per-component wall-clock timing for FESOM3.
    !
    ! A small fixed registry of named timers; timer_start/timer_stop accumulate
    ! wall time (MPI_Wtime) per component; timer_report reduces mean/min/max across
    ! ranks and prints a hierarchical FESOM2-style breakdown (top-level components
    ! plus the ocean sub-steps), with per-step ms and % of the loop. timer_reset
    ! zeroes the accumulators so a caller can also emit periodic (windowed) monitoring.
    !
    ! Overhead is negligible: each timer_start/stop is one MPI_Wtime (~tens of ns) and
    ! NO collective — the single reduction happens only in timer_report (end of run, or
    ! once per monitoring window). Enabled by default; FESOM3_TIMING=0 silences output.
    !
    ! Mirrors FESOM2's `rtime_*` breakdown (fesom_module.F90 + oce_ale.F90 t0..t10):
    !   mix/pres, dyn u,v,w, dyn ssh, solve ssh (subset of dyn ssh), GM/Redi, tracer.
    !===========================================================================
    use, intrinsic :: iso_fortran_env, only: real64
    use mpi
    implicit none
    private

    ! ---- timer ids (pre-order tree: each parent is immediately followed by its children, so the
    !      report prints with correct indentation by simple iteration. Use the NAMES at call sites;
    !      the integer values may be reordered freely.) ------------------------
    integer, parameter, public :: TMR_OCEAN2ICE   =  1   ! top-level driver components
    integer, parameter, public :: TMR_FORCING     =  2
    integer, parameter, public :: TMR_FRC_SBC     =  3   ! forcing sub-steps (children of FORCING)
    integer, parameter, public :: TMR_FRC_OPEN    =  4   ! crossing read split (children of FRC_SBC)
    integer, parameter, public :: TMR_FRC_READ    =  5   !   nc slice read + DEFLATE inflate + halo
    integer, parameter, public :: TMR_FRC_INTRP2  =  6   !   spatial bilinear coeff build (do ii)
    integer, parameter, public :: TMR_FRC_INTERP  =  7
    integer, parameter, public :: TMR_FRC_BULK    =  8
    integer, parameter, public :: TMR_FRC_STRESS  =  9
    integer, parameter, public :: TMR_ICE         = 10
    integer, parameter, public :: TMR_FLUXES      = 11
    integer, parameter, public :: TMR_STEP_OCE    = 12   ! umbrella: whole ocean step
    integer, parameter, public :: TMR_OCE_MIXPRES = 13   ! ocean sub-steps (children of STEP_OCE)
    integer, parameter, public :: TMR_OCE_DYN     = 14
    integer, parameter, public :: TMR_OCE_SSH     = 15
    integer, parameter, public :: TMR_OCE_SOLVE   = 16   ! SUBSET of SSH (the CG solve) — not summed
    integer, parameter, public :: TMR_OCE_SHSYNC  = 17   ! SUBSET of SSH (update_vel's shared-element uv sync)
    integer, parameter, public :: TMR_OCE_GMREDI  = 18
    integer, parameter, public :: TMR_OCE_TRACER  = 19
    integer, parameter, public :: NTIMER          = 19

    character(len=24), parameter :: TNAME(NTIMER) = [character(len=24) :: &
        'ocean2ice', 'forcing', 'sbc / crossing', 'open/reopen', 'read+inflate', 'interp loop', &
        'time-interp', 'bulk NCAR', 'wind/ice stress', &
        'ice', 'oce_fluxes', 'step_oce (ocean)', &
        'mix, pres, EOS', 'dynamics u,v,w', 'dynamics ssh', '(of which) ssh solve', &
        '(of which) shared sync', 'GM / Redi', 'tracer' ]
    ! parent id for the indented report (0 = top-level component)
    integer, parameter :: TPARENT(NTIMER) = [ 0, 0, TMR_FORCING, TMR_FRC_SBC, TMR_FRC_SBC, TMR_FRC_SBC, &
        TMR_FORCING, TMR_FORCING, TMR_FORCING, &
        0, 0, 0, &
        TMR_STEP_OCE, TMR_STEP_OCE, TMR_STEP_OCE, TMR_OCE_SSH, TMR_OCE_SSH, TMR_STEP_OCE, TMR_STEP_OCE ]
    ! .true. => this timer is a subset of its parent (don't add to the parent's child-sum check)
    logical, parameter :: TSUBSET(NTIMER) = [ .false.,.false.,.false., .true.,.true.,.true., &
        .false.,.false.,.false., &
        .false.,.false.,.false., &
        .false.,.false.,.false., .true., .true., .false.,.false. ]

    real(real64) :: tacc(NTIMER) = 0.0_real64   ! accumulated seconds
    real(real64) :: tbeg(NTIMER) = 0.0_real64   ! start stamp of an open interval
    integer(8)   :: tcnt(NTIMER) = 0_8          ! number of intervals (for sanity)
    logical      :: enabled = .true.            ! FESOM3_TIMING=0 silences the report

    public :: timer_init, timer_start, timer_stop, timer_report, timer_reset

contains

    subroutine timer_init()
        ! Read the on/off switch once (default ON). Accumulation is always cheap; this only
        ! controls whether timer_report prints.
        character(len=8) :: env
        integer :: ln, ios
        call get_environment_variable('FESOM3_TIMING', env, length=ln, status=ios)
        if (ios == 0 .and. ln > 0) enabled = (trim(env) /= '0')
        call timer_reset()
    end subroutine timer_init

    subroutine timer_reset()
        tacc = 0.0_real64
        tcnt = 0_8
    end subroutine timer_reset

    subroutine timer_start(id)
        integer, intent(in) :: id
        tbeg(id) = MPI_Wtime()
    end subroutine timer_start

    subroutine timer_stop(id)
        integer, intent(in) :: id
        tacc(id) = tacc(id) + (MPI_Wtime() - tbeg(id))
        tcnt(id) = tcnt(id) + 1_8
    end subroutine timer_stop

    subroutine timer_report(comm, mype, npes, nsteps, label)
        !-----------------------------------------------------------------------
        ! Collective: ALL ranks must call. Reduces mean(=sum/npes)/min/max of each
        ! timer across `comm`, then rank 0 prints the hierarchical breakdown with
        ! per-step ms and % of the loop. `nsteps` normalizes to per-step; `label`
        ! tags the block (e.g. 'FINAL' or 'STEP 4800').
        !-----------------------------------------------------------------------
        integer,          intent(in) :: comm, mype, npes, nsteps
        character(len=*), intent(in) :: label
        real(real64) :: tmean(NTIMER), tmin(NTIMER), tmax(NTIMER), loop_s, denom
        real(real64) :: perstep_ms, pct
        integer :: i, ierr, ns

        if (.not. enabled) return
        ns = max(nsteps, 1)

        tmean = tacc; call MPI_Allreduce(MPI_IN_PLACE, tmean, NTIMER, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
        tmean = tmean / real(max(npes,1), real64)
        tmin  = tacc; call MPI_Allreduce(MPI_IN_PLACE, tmin,  NTIMER, MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)
        tmax  = tacc; call MPI_Allreduce(MPI_IN_PLACE, tmax,  NTIMER, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)

        if (mype /= 0) return

        ! loop denominator = sum of the top-level (non-subset) components' mean time
        loop_s = 0.0_real64
        do i = 1, NTIMER
            if (TPARENT(i) == 0) loop_s = loop_s + tmean(i)
        end do
        denom = max(loop_s, tiny(1.0_real64))

        write(*,'(a)') ''
        write(*,'(a)') '  ==================================================================================='
        write(*,'(3a,i0,a,i0,a)') '  TIMING [', trim(label), '] -- ', ns, ' step(s), ', npes, ' rank(s)'
        write(*,'(a)') '  component                       mean_ms/step    min_ms/step    max_ms/step   %loop'
        write(*,'(a)') '  -----------------------------------------------------------------------------------'
        do i = 1, NTIMER
            perstep_ms = 1.0e3_real64 * tmean(i) / ns
            pct        = 100.0_real64 * tmean(i) / denom
            if (TPARENT(i) == 0) then
                write(*,'(a,a20,3f15.4,f8.1)') '  ', TNAME(i), &
                    perstep_ms, 1.0e3_real64*tmin(i)/ns, 1.0e3_real64*tmax(i)/ns, pct
            else if (TSUBSET(i)) then
                write(*,'(a,a24,3f15.4)') '        > ', TNAME(i), &
                    perstep_ms, 1.0e3_real64*tmin(i)/ns, 1.0e3_real64*tmax(i)/ns
            else
                write(*,'(a,a22,3f15.4,f8.1)') '      > ', TNAME(i), &
                    perstep_ms, 1.0e3_real64*tmin(i)/ns, 1.0e3_real64*tmax(i)/ns, pct
            end if
        end do
        write(*,'(a)') '  -----------------------------------------------------------------------------------'
        write(*,'(a,a20,f15.4)') '  ', 'LOOP TOTAL', 1.0e3_real64*loop_s/ns
        write(*,'(a)') '  ==================================================================================='
    end subroutine timer_report

end module mod_timer
