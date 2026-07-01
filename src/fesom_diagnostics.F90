!=============================================================================!
!                   FESOM2 Scalar-Diagnostics Module
!=============================================================================!
! Lightweight, reusable facility to record arbitrary per-call SCALAR quantities
! (solver iterations, residuals, CFL numbers, event counts, ...) and report
! them, MPI-reduced, into fesom.stats next to the timing table.
!
! It mirrors the fesom_profiler patterns (name->index registry, per-entry
! accumulation, reduce-then-print) but accumulates VALUES rather than wall time.
! Storage is double precision (PROF_WP) regardless of the model working
! precision WP, so SP and DP diagnostics are directly comparable.
!
! Usage (guard call sites with #if defined(FESOM_PROFILING) for zero cost when
! profiling is compiled out, matching the profiler's call sites):
!     use fesom_diagnostics
!     call diag_count("ssh_cg_iters", iter)          ! integer
!     call diag_count("ssh_cg_resid", real(res,8))   ! real
!
! Report is emitted by diag_report(), called from fesom_profiler_report().
!=============================================================================!

module fesom_diagnostics
    use mpi
    use, intrinsic :: iso_fortran_env, only: real32, real64
    implicit none
    private

    integer, parameter :: PROF_WP = selected_real_kind(15, 307)  ! double precision
    integer, parameter :: MAX_DIAGS = 64

    type :: diag_stat
        character(len=100) :: name = ""
        real(kind=PROF_WP) :: sum_val   = 0.0_PROF_WP
        real(kind=PROF_WP) :: sumsq_val = 0.0_PROF_WP
        real(kind=PROF_WP) :: min_val   = HUGE(1.0_PROF_WP)
        real(kind=PROF_WP) :: max_val   = -HUGE(1.0_PROF_WP)
        integer            :: count     = 0
    end type diag_stat

    type(diag_stat), save :: diags(MAX_DIAGS)
    integer, save :: num_diags = 0
    logical, save :: diag_on   = .false.

    public :: diag_init, diag_enabled, diag_reset, diag_count, diag_report

    ! Generic: pass an integer, a real(real32), or a real(real64) with no cast
    interface diag_count
        module procedure diag_count_i
        module procedure diag_count_r4
        module procedure diag_count_r8
    end interface diag_count

contains

    subroutine diag_init(enabled)
        logical, intent(in) :: enabled
        diag_on = enabled
        call diag_reset()
    end subroutine diag_init

    logical function diag_enabled()
        diag_enabled = diag_on
    end function diag_enabled

    subroutine diag_reset()
        integer :: i
        do i = 1, MAX_DIAGS
            diags(i)%name      = ""
            diags(i)%sum_val   = 0.0_PROF_WP
            diags(i)%sumsq_val = 0.0_PROF_WP
            diags(i)%min_val   =  HUGE(1.0_PROF_WP)
            diags(i)%max_val   = -HUGE(1.0_PROF_WP)
            diags(i)%count     = 0
        end do
        num_diags = 0
    end subroutine diag_reset

    !------------------------------------------------------------------ generics
    subroutine diag_count_i(name, value)
        character(len=*), intent(in) :: name
        integer,          intent(in) :: value
        call accumulate(name, real(value, PROF_WP))
    end subroutine diag_count_i

    subroutine diag_count_r4(name, value)
        character(len=*),   intent(in) :: name
        real(kind=real32),  intent(in) :: value
        call accumulate(name, real(value, PROF_WP))
    end subroutine diag_count_r4

    subroutine diag_count_r8(name, value)
        character(len=*),   intent(in) :: name
        real(kind=real64),  intent(in) :: value
        call accumulate(name, real(value, PROF_WP))
    end subroutine diag_count_r8

    !--------------------------------------------------------------- accumulate
    subroutine accumulate(name, val)
        character(len=*),   intent(in) :: name
        real(kind=PROF_WP), intent(in) :: val
        integer :: i
        if (.not. diag_on) return
        i = find_or_create_diag(trim(name))
        if (i < 1) return
        diags(i)%sum_val   = diags(i)%sum_val   + val
        diags(i)%sumsq_val = diags(i)%sumsq_val + val*val
        diags(i)%min_val   = min(diags(i)%min_val, val)
        diags(i)%max_val   = max(diags(i)%max_val, val)
        diags(i)%count     = diags(i)%count + 1
    end subroutine accumulate

    integer function find_or_create_diag(name) result(index)
        character(len=*), intent(in) :: name
        integer :: i
        index = -1
        do i = 1, num_diags
            if (trim(diags(i)%name) == name) then
                index = i
                return
            endif
        end do
        if (num_diags < MAX_DIAGS) then
            num_diags = num_diags + 1
            index = num_diags
            diags(index)%name = trim(name)
        endif
    end function find_or_create_diag

    !------------------------------------------------------------------- report
    ! Reduce across ranks and (on rank 0) append a PER-CALL DIAGNOSTICS block.
    ! Assumes all ranks created the same diagnostics in the same order (the same
    ! assumption fesom_profiler_report makes for timing sections); true for
    ! diagnostics recorded on a globally-collective code path (e.g. the SSH CG).
    subroutine diag_report(unit, mpi_comm, mpi_rank)
        integer, intent(in) :: unit, mpi_comm, mpi_rank
        integer :: ierr, i, gnmax
        real(kind=PROF_WP), allocatable :: lsum(:), lmin(:), lmax(:), lsq(:)
        real(kind=PROF_WP), allocatable :: gsum(:), gmin(:), gmax(:), gsq(:)
        integer, allocatable :: lcnt(:), gcnt(:)
        real(kind=PROF_WP) :: mean, var, std

        if (.not. diag_on) return

        ! Guard against ranks disagreeing on the diagnostic count.
        call MPI_Allreduce(num_diags, gnmax, 1, MPI_INTEGER, MPI_MAX, mpi_comm, ierr)
        if (gnmax == 0) return
        if (num_diags /= gnmax) return  ! inconsistent registries -> skip safely

        allocate(lsum(num_diags), lmin(num_diags), lmax(num_diags), lsq(num_diags), lcnt(num_diags))
        allocate(gsum(num_diags), gmin(num_diags), gmax(num_diags), gsq(num_diags), gcnt(num_diags))
        do i = 1, num_diags
            lsum(i) = diags(i)%sum_val
            lsq(i)  = diags(i)%sumsq_val
            lmin(i) = diags(i)%min_val
            lmax(i) = diags(i)%max_val
            lcnt(i) = diags(i)%count
        end do

        call MPI_Allreduce(lsum, gsum, num_diags, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm, ierr)
        call MPI_Allreduce(lsq,  gsq,  num_diags, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm, ierr)
        call MPI_Allreduce(lcnt, gcnt, num_diags, MPI_INTEGER,          MPI_SUM, mpi_comm, ierr)
        call MPI_Allreduce(lmin, gmin, num_diags, MPI_DOUBLE_PRECISION, MPI_MIN, mpi_comm, ierr)
        call MPI_Allreduce(lmax, gmax, num_diags, MPI_DOUBLE_PRECISION, MPI_MAX, mpi_comm, ierr)

        if (mpi_rank == 0) then
            write(unit, '(A)') ''
            write(unit, '(A)') repeat('=', 90)
            write(unit, '(A)') '================ PER-CALL DIAGNOSTICS (scalar counters) ================'
            write(unit, '(A)') 'Mean/call = sum over all calls & ranks / total call count.'
            write(unit, '(A)') 'Min/Max across ranks & calls; for a globally-identical value Min==Max.'
            write(unit, '(A)') repeat('=', 90)
            write(unit, '(A28,A14,A12,A12,A12,A12)') &
                'Name', 'Mean/call', 'Min', 'Max', 'Std', 'Count'
            write(unit, '(A)') repeat('-', 90)
            do i = 1, num_diags
                if (gcnt(i) <= 0) cycle
                mean = gsum(i) / real(gcnt(i), PROF_WP)
                var  = gsq(i) / real(gcnt(i), PROF_WP) - mean*mean
                if (var < 0.0_PROF_WP) var = 0.0_PROF_WP
                std  = sqrt(var)
                write(unit, '(A28,F14.4,F12.4,F12.4,F12.4,I12)') &
                    adjustl(diags(i)%name), mean, gmin(i), gmax(i), std, gcnt(i)
            end do
            write(unit, '(A)') repeat('-', 90)
            write(unit, '(A)') ''
        endif

        deallocate(lsum, lmin, lmax, lsq, lcnt, gsum, gmin, gmax, gsq, gcnt)
    end subroutine diag_report

end module fesom_diagnostics
