!=============================================================================!
!                   FESOM2 Enhanced Profiler Module
!=============================================================================!
! Enhanced profiling system with decorator pattern support, comprehensive
! statistics, and load balance analysis for parallel computations.
!
! Key Features:
! - Decorator-style profiling with FESOM_PROFILE_START/END macros
! - Advanced statistics: mean, std dev, load balance metrics
! - Thread-safe design for OpenMP compatibility  
! - Zero overhead when disabled via preprocessor
! - Comprehensive load balance analysis
!=============================================================================!

module fesom_profiler
    use g_config, only: runid, ResultPath
    implicit none
    private
    
    ! Define profiler-specific precision for time measurements
    ! Use double precision to ensure compatibility with MPI_Wtime and MPI operations
    integer, parameter :: PROF_WP = selected_real_kind(15, 307)  ! Double precision
    
    ! MPI data type for PROF_WP - will be set to MPI_DOUBLE_PRECISION at runtime
    ! This ensures compatibility with our double precision timing variables
    
    ! Public interface
    public :: fesom_profiler_init, fesom_profiler_finalize
    public :: fesom_profiler_start, fesom_profiler_end
    public :: fesom_profiler_report, fesom_profiler_reset
    public :: fesom_profiler_enabled, fesom_profiler_set_timesteps
    public :: fesom_profiler_set_timestep_size
    
    ! Maximum number of profiling sections
    integer, parameter :: MAX_PROFILE_SECTIONS = 256
    
    ! Profiling statistics type
    type :: profile_stats
        character(len=100) :: name = ""
        real(kind=PROF_WP) :: total_time = 0.0_PROF_WP
        real(kind=PROF_WP) :: min_time = HUGE(1.0_PROF_WP)
        real(kind=PROF_WP) :: max_time = 0.0_PROF_WP
        real(kind=PROF_WP) :: sum_squares = 0.0_PROF_WP     ! For std deviation calculation
        integer :: call_count = 0
        logical :: is_active = .false.
        real(kind=PROF_WP) :: start_time = 0.0_PROF_WP
        ! Load balance metrics (computed during finalization)
        real(kind=PROF_WP) :: mean_time = 0.0_PROF_WP
        real(kind=PROF_WP) :: std_dev = 0.0_PROF_WP
        real(kind=PROF_WP) :: load_imbalance = 0.0_PROF_WP  ! (max-min)/mean * 100
        real(kind=PROF_WP) :: efficiency = 0.0_PROF_WP      ! mean/max * 100
        integer :: participating_ranks = 0
        ! Call hierarchy tracking
        character(len=100) :: parent_name = ""
        integer :: nesting_level = 0
        ! Cumulative statistics across ranks
        real(kind=PROF_WP) :: total_time_across_ranks = 0.0_PROF_WP
        real(kind=PROF_WP) :: min_total_time = HUGE(1.0_PROF_WP)
        real(kind=PROF_WP) :: max_total_time = 0.0_PROF_WP
    end type profile_stats
    
    ! Module variables
    type(profile_stats), save :: profiles(MAX_PROFILE_SECTIONS)
    integer, save :: num_profiles = 0
    logical, save :: profiler_initialized = .false.
    logical, save :: profiler_enabled = .true.
    
    ! For nested profiling support
    integer, parameter :: MAX_CALL_STACK_DEPTH = 32
    character(len=100), save :: call_stack(MAX_CALL_STACK_DEPTH)
    integer, save :: call_stack_depth = 0
    
    ! System information
    integer, save :: total_timesteps = 0
    integer, save :: total_ranks = 0
    real(kind=PROF_WP), save :: timestep_size = 0.0_PROF_WP  ! Model timestep in seconds
    integer, save :: omp_threads = 1
    character(len=200), save :: system_info = ""
    
contains

    !=========================================================================
    ! Initialize the profiler system
    !=========================================================================
    subroutine fesom_profiler_init(enable_profiling)
        logical, intent(in), optional :: enable_profiling
        integer :: i
        
        if (present(enable_profiling)) then
            profiler_enabled = enable_profiling
        endif
        
        if (.not. profiler_enabled) return

        ! Detect OpenMP thread count
#ifdef _OPENMP
        !$ omp_threads = omp_get_max_threads()
#else
        omp_threads = 1
#endif
        
        ! Initialize all profiles
        do i = 1, MAX_PROFILE_SECTIONS
            profiles(i)%name = ""
            profiles(i)%total_time = 0.0_PROF_WP
            profiles(i)%min_time = HUGE(1.0_PROF_WP)
            profiles(i)%max_time = 0.0_PROF_WP
            profiles(i)%sum_squares = 0.0_PROF_WP
            profiles(i)%call_count = 0
            profiles(i)%is_active = .false.
            profiles(i)%start_time = 0.0_PROF_WP
            profiles(i)%mean_time = 0.0_PROF_WP
            profiles(i)%std_dev = 0.0_PROF_WP
            profiles(i)%load_imbalance = 0.0_PROF_WP
            profiles(i)%efficiency = 0.0_PROF_WP
            profiles(i)%participating_ranks = 0
        end do
        
        num_profiles = 0
        call_stack_depth = 0
        profiler_initialized = .true.
    end subroutine fesom_profiler_init

    !=========================================================================
    ! Check if profiler is enabled
    !=========================================================================
    function fesom_profiler_enabled() result(enabled)
        logical :: enabled
        enabled = profiler_enabled .and. profiler_initialized
    end function fesom_profiler_enabled

    !=========================================================================
    ! Set total timesteps for reporting
    !=========================================================================
    subroutine fesom_profiler_set_timesteps(timesteps)
        integer, intent(in) :: timesteps
        if (.not. fesom_profiler_enabled()) return
        total_timesteps = timesteps
    end subroutine fesom_profiler_set_timesteps

    !=========================================================================
    ! Set timestep size in seconds for SYPD calculation
    !=========================================================================
    subroutine fesom_profiler_set_timestep_size(dt_seconds)
        real(kind=PROF_WP), intent(in) :: dt_seconds
        if (.not. fesom_profiler_enabled()) return
        timestep_size = dt_seconds
    end subroutine fesom_profiler_set_timestep_size

    !=========================================================================
    ! Start profiling a section
    !=========================================================================
    subroutine fesom_profiler_start(section_name, source_file, line_number)
        include "mpif.h"
        character(len=*), intent(in) :: section_name
        character(len=*), intent(in), optional :: source_file
        integer, intent(in), optional :: line_number
        integer :: profile_index
        character(len=200) :: full_name
        
        if (.not. fesom_profiler_enabled()) return
        
        ! Create full name with source location if provided
        if (present(source_file) .and. present(line_number)) then
            write(full_name, '(A,A,A,I0)') trim(section_name), ' (', trim(source_file), line_number
        else
            full_name = trim(section_name)
        endif
        
        profile_index = find_or_create_profile(trim(section_name))
        
        ! Set parent relationship and nesting level ONLY on first creation
        ! Don't overwrite if already set (preserves hierarchy from first call)
        if (profiles(profile_index)%parent_name == "") then
            if (call_stack_depth > 0) then
                profiles(profile_index)%parent_name = trim(call_stack(call_stack_depth))
                profiles(profile_index)%nesting_level = call_stack_depth
            else
                profiles(profile_index)%parent_name = ""
                profiles(profile_index)%nesting_level = 0
            endif
        endif
        
        ! Add to call stack for nested profiling
        if (call_stack_depth < MAX_CALL_STACK_DEPTH) then
            call_stack_depth = call_stack_depth + 1
            call_stack(call_stack_depth) = trim(section_name)
        endif
        
        if (profiles(profile_index)%is_active) then
            ! Warning: nested call to same section - this is fine for recursive calls
        endif
        
        profiles(profile_index)%is_active = .true.
        profiles(profile_index)%start_time = MPI_Wtime()
    end subroutine fesom_profiler_start

    !=========================================================================
    ! End profiling a section
    !=========================================================================
    subroutine fesom_profiler_end(section_name, source_file, line_number)
        include "mpif.h"
        character(len=*), intent(in) :: section_name
        character(len=*), intent(in), optional :: source_file
        integer, intent(in), optional :: line_number
        integer :: profile_index
        real(kind=PROF_WP) :: end_time, elapsed_time
        
        if (.not. fesom_profiler_enabled()) return
        
        profile_index = find_profile(trim(section_name))
        if (profile_index == -1) then
            ! Error: trying to end a section that wasn't started
            return
        endif
        
        if (.not. profiles(profile_index)%is_active) then
            ! Warning: trying to end a section that's not active
            return
        endif
        
        end_time = MPI_Wtime()
        elapsed_time = end_time - profiles(profile_index)%start_time
        
        ! Update statistics
        profiles(profile_index)%total_time = profiles(profile_index)%total_time + elapsed_time
        profiles(profile_index)%min_time = min(profiles(profile_index)%min_time, elapsed_time)
        profiles(profile_index)%max_time = max(profiles(profile_index)%max_time, elapsed_time)
        profiles(profile_index)%sum_squares = profiles(profile_index)%sum_squares + elapsed_time**2
        profiles(profile_index)%call_count = profiles(profile_index)%call_count + 1
        profiles(profile_index)%is_active = .false.
        
        ! Remove from call stack
        if (call_stack_depth > 0 .and. &
            trim(call_stack(call_stack_depth)) == trim(section_name)) then
            call_stack_depth = call_stack_depth - 1
        endif
    end subroutine fesom_profiler_end

    !=========================================================================
    ! Generate comprehensive profiling report
    !=========================================================================
    subroutine fesom_profiler_report(mpi_comm, mpi_rank, output_unit)
        include "mpif.h"
        integer, intent(in) :: mpi_comm, mpi_rank
        integer, intent(in), optional :: output_unit
        integer :: unit, ierr, i, npes
        real(kind=PROF_WP), allocatable :: local_data(:), global_data(:)
        integer, allocatable :: local_counts(:), global_counts(:)
        logical :: section_has_data, file_opened
        character(len=300) :: stats_filename
        
        if (.not. fesom_profiler_enabled()) return
        
        file_opened = .false.
        unit = 6  ! stdout default
        
        if (present(output_unit)) then
            unit = output_unit
        else
            ! Only rank 0 opens the stats file
            if (mpi_rank == 0) then
                if (trim(runid) == 'fesom') then
                    stats_filename = trim(ResultPath) // 'fesom.stats'
                else
                    stats_filename = trim(ResultPath) // 'fesom.' // trim(runid) // '.stats'
                endif
                
                open(newunit=unit, file=trim(stats_filename), action='write', status='replace', iostat=ierr)
                if (ierr /= 0) then
                    write(*,*) 'Warning: Cannot open stats file ', trim(stats_filename), ', using stdout'
                    unit = 6
                else
                    file_opened = .true.
                endif
            endif
        endif
        
        call MPI_Comm_size(mpi_comm, npes, ierr)
        
        ! Prepare data for MPI reduction - separate arrays for different operations
        allocate(local_data(4 * num_profiles), global_data(4 * num_profiles))
        allocate(local_counts(num_profiles), global_counts(num_profiles))
        
        ! Pack local data: total_time, min_total_time, max_total_time, sum_squares
        do i = 1, num_profiles
            local_data(4*i-3) = profiles(i)%total_time  ! For averaging (wall clock)
            local_data(4*i-2) = profiles(i)%total_time  ! For min across ranks
            local_data(4*i-1) = profiles(i)%total_time  ! For max across ranks
            local_data(4*i)   = profiles(i)%sum_squares
            local_counts(i)   = profiles(i)%call_count
        end do
        
        ! Get averages for wall clock time calculation
        ! Use MPI_DOUBLE_PRECISION since PROF_WP is double precision
        call MPI_Allreduce(local_data(1::4), global_data(1::4), num_profiles, &
                          MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm, ierr)
        call MPI_Allreduce(local_counts, global_counts, num_profiles, &
                          MPI_INTEGER, MPI_SUM, mpi_comm, ierr)
        
        ! Get min/max total times across ranks for load balance analysis
        call MPI_Allreduce(local_data(2::4), global_data(2::4), num_profiles, &
                          MPI_DOUBLE_PRECISION, MPI_MIN, mpi_comm, ierr)
        call MPI_Allreduce(local_data(3::4), global_data(3::4), num_profiles, &
                          MPI_DOUBLE_PRECISION, MPI_MAX, mpi_comm, ierr)
        
        ! Only rank 0 prints the report
        if (mpi_rank == 0) then
            call print_detailed_report(unit, global_data, global_counts, npes)
        endif
        
        ! Close the stats file if we opened it
        if (file_opened .and. mpi_rank == 0) then
            close(unit)
        endif
        
        deallocate(local_data, global_data, local_counts, global_counts)
    end subroutine fesom_profiler_report

    !=========================================================================
    ! Print detailed profiling report
    !=========================================================================
    subroutine print_detailed_report(unit, global_data, global_counts, npes)
        integer, intent(in) :: unit, npes
        real(kind=PROF_WP), intent(in) :: global_data(:)
        integer, intent(in) :: global_counts(:)
        integer :: i
        real(kind=PROF_WP) :: total_runtime_for_percent, wall_clock_time
        logical :: found_total_for_percent
        
        ! Print header in FESOM style
        write(unit, '(A)') ''
        write(unit, '(A)') repeat('=', 90)
        write(unit, '(A)') '======== FESOM2 ENHANCED PROFILING REPORT (detailed timing statistics) ========'
        if (total_timesteps > 0) then
            write(unit, '(A,I0,A,I0,A)') '   ', npes, ' ranks, ', total_timesteps, ' timesteps'
        else
            write(unit, '(A,I0,A)') '   ', npes, ' MPI ranks'
        endif
        write(unit, '(A)') repeat('=', 90)
        
        ! Add comprehensive explanation section
        write(unit, '(A)') ''
        write(unit, '(A)') 'METRICS EXPLANATION:'
        write(unit, '(A)') '-------------------'
        write(unit, '(A)') 'Mean(s): Average wall clock time across all MPI ranks'
        write(unit, '(A)') 'Min(s)/Max(s): Fastest/slowest total time among all ranks'
        write(unit, '(A)') 'Calls: Total number of function calls across all ranks'
        write(unit, '(A)') ''
        write(unit, '(A)') 'PERCENTAGE METRICS:'
        write(unit, '(A)') '%Total: Percentage of absolute total runtime (init + runloop + finalize)'
        write(unit, '(A)') '%Parent: Percentage relative to immediate parent section'
        write(unit, '(A)') '  → Top-level: same as %Total (relative to absolute total)'
        write(unit, '(A)') '  → Nested: relative to parent (e.g., oce_mix_pres % of oce_timestep_ale)'
        write(unit, '(A)') ''
        write(unit, '(A)') 'LOAD BALANCE METRICS:'
        write(unit, '(A)') 'RngImb(%) = (Max - Min) / Mean × 100  [Range-based imbalance]'
        write(unit, '(A)') '  → How much slower is the worst rank vs average?'
        write(unit, '(A)') '  → Critical for HPC: slowest rank determines completion time'
        write(unit, '(A)') 'StdImb(%) = StdDev / Mean × 100       [Distribution-based imbalance]'
        write(unit, '(A)') '  → Overall variability across all ranks (coefficient of variation)'
        write(unit, '(A)') '  → Shows whether imbalance affects few ranks or is widespread'
        write(unit, '(A)') ''
        write(unit, '(A)') 'INTERPRETATION EXAMPLES:'
        write(unit, '(A)') '• High RngImb + Low StdImb  → One outlier rank, others well-balanced'
        write(unit, '(A)') '• High RngImb + High StdImb → General imbalance across many ranks'
        write(unit, '(A)') '• Low RngImb + High StdImb  → Multiple performance clusters, no single outlier'
        write(unit, '(A)') ''
        if (total_timesteps > 0 .and. timestep_size > 0.0_PROF_WP) then
            write(unit, '(A,F8.1,A)') 'SYPD CALCULATION: Simulated Years Per Day = (timestep_size_seconds × timesteps) / (365.25 × runtime_seconds)'
            write(unit, '(A,F8.1,A,I0,A)') '  → Based on ', timestep_size, 's timesteps, ', total_timesteps, ' steps total'
        endif
        write(unit, '(A)') repeat('=', 90)
        
        ! Print detailed header with clean formatting
        write(unit, '(A)') ''
        write(unit, '(A35,A15,A15,A15,A8,A8,A8,A10,A10)') &
            'Section Name', 'Mean(s)', 'Min(s)', 'Max(s)', 'Calls', '%Total', '%Parent', 'RngImb(%)', 'StdImb(%)'
        write(unit, '(A)') repeat('-', 124)
        
        ! Calculate absolute total runtime (sum of all top-level sections)
        total_runtime_for_percent = 0.0_PROF_WP
        
        do i = 1, num_profiles
            if (trim(profiles(i)%name) == "") cycle
            
            ! Only sum top-level sections (those with no parent)
            if (trim(profiles(i)%parent_name) == "") then
                wall_clock_time = global_data(4*i-3) / real(npes, PROF_WP)
                total_runtime_for_percent = total_runtime_for_percent + wall_clock_time
            endif
        end do
        
        ! Fallback to avoid division by zero
        if (total_runtime_for_percent <= 0.0_PROF_WP) then
            total_runtime_for_percent = 1.0_PROF_WP
        endif
        
        ! Print sections with all statistics
        call print_detailed_sections(unit, global_data, global_counts, npes, 0, "", total_runtime_for_percent, total_runtime_for_percent)
        
        write(unit, '(A)') repeat('-', 124)
        write(unit, '(A)') ''
        
        ! Load balance summary section
        write(unit, '(A)') repeat('=', 90)
        write(unit, '(A)') '========================== LOAD BALANCE ANALYSIS =========================='
        call print_load_balance_summary(unit, global_data, global_counts, npes)
        write(unit, '(A)') repeat('=', 90)
        write(unit, '(A)') ''
        
        ! Find the total runtime from runloop or main sections
        call print_benchmark_summary_with_total_runtime(unit, global_data, global_counts, npes)
    end subroutine print_detailed_report

    !=========================================================================
    ! Print detailed sections with all statistics and hierarchy
    !=========================================================================
    recursive subroutine print_detailed_sections(unit, global_data, global_counts, npes, level, parent, total_runtime, parent_runtime)
        integer, intent(in) :: unit, npes, level
        real(kind=PROF_WP), intent(in) :: global_data(:)
        integer, intent(in) :: global_counts(:)
        character(len=*), intent(in) :: parent
        real(kind=PROF_WP), intent(in) :: total_runtime
        real(kind=PROF_WP), intent(in), optional :: parent_runtime
        integer :: i, j, total_calls
        real(kind=PROF_WP) :: wall_clock_time, min_total, max_total, load_imbalance
        real(kind=PROF_WP) :: mean_time, min_call_time, max_call_time, std_dev, sum_squares, std_imbalance
        real(kind=PROF_WP) :: percent_total, percent_parent, current_parent_runtime
        character(len=100) :: display_name
        character(len=15) :: indent_prefix
        
        ! Print sections at this level
        do i = 1, num_profiles
            if (trim(profiles(i)%name) == "") cycle
            ! Only print sections that belong at this level
            if (level == 0) then
                ! At top level, only print sections with no parent
                if (trim(profiles(i)%parent_name) /= "") cycle
            else
                ! At nested levels, only print sections whose parent matches
                if (trim(profiles(i)%parent_name) /= trim(parent)) cycle
            endif
            
            ! Extract all statistics from global data
            wall_clock_time = global_data(4*i-3) / real(npes, PROF_WP)  ! Wall clock time
            min_total = global_data(4*i-2)                         ! Min total across ranks
            max_total = global_data(4*i-1)                         ! Max total across ranks
            sum_squares = global_data(4*i)                         ! Sum of squares
            total_calls = global_counts(i)
            
            if (total_calls == 0) cycle
            
            ! Calculate per-call statistics
            mean_time = wall_clock_time / real(total_calls/npes, PROF_WP)
            
            ! Calculate min/max per-call times (approximation)
            if (total_calls > 0) then
                min_call_time = min_total / real(max(1, total_calls/npes), PROF_WP)
                max_call_time = max_total / real(max(1, total_calls/npes), PROF_WP)
            else
                min_call_time = 0.0_PROF_WP
                max_call_time = 0.0_PROF_WP
            endif
            
            ! Calculate standard deviation (use mean per-call time for high-frequency calls)
            if (total_calls > npes .and. wall_clock_time > 0.0_PROF_WP) then
                if (total_calls > 100) then
                    ! For high-frequency calls, use mean-based std dev to avoid overflow
                    std_dev = (max_total - min_total) / 4.0_PROF_WP  ! Approximation: range/4 ≈ std dev
                else
                    std_dev = sqrt(abs(sum_squares/real(npes,PROF_WP) - (wall_clock_time**2/real(npes,PROF_WP))))
                endif
            else
                std_dev = 0.0_PROF_WP
            endif
            
            ! Calculate load imbalance
            if (wall_clock_time > 0.0_PROF_WP) then
                load_imbalance = (max_total - min_total) / wall_clock_time * 100.0_PROF_WP
            else
                load_imbalance = 0.0_PROF_WP
            endif
            
            ! Create indented display name with FESOM-style hierarchy
            if (level == 0) then
                display_name = trim(profiles(i)%name)
            else
                ! Use FESOM-style indentation with '>' symbols
                ! level=1 gets "  > ", level=2 gets "    > ", etc.
                write(indent_prefix, '(A)') repeat("  ", level) // "> "
                write(display_name, '(A,A)') trim(indent_prefix), trim(profiles(i)%name)
            endif
            
            ! Calculate std imbalance as percentage of mean
            if (wall_clock_time > 0.0_PROF_WP) then
                std_imbalance = (std_dev / wall_clock_time) * 100.0_PROF_WP
            else
                std_imbalance = 0.0_PROF_WP
            endif
            
            ! Calculate percentage of absolute total runtime
            if (total_runtime > 0.0_PROF_WP) then
                percent_total = (wall_clock_time / total_runtime) * 100.0_PROF_WP
            else
                percent_total = 0.0_PROF_WP
            endif
            
            ! Calculate percentage relative to parent
            if (present(parent_runtime) .and. parent_runtime > 0.0_PROF_WP) then
                percent_parent = (wall_clock_time / parent_runtime) * 100.0_PROF_WP
            else
                ! For top-level sections, parent % same as total %
                percent_parent = percent_total
            endif
            
            ! Print line with dual percentages and dual load balance metrics
            if (std_dev < 9999.0_PROF_WP) then
                write(unit, '(A35,3F15.4,I8,F8.1,F8.1,F10.1,F10.1)') &
                    display_name, wall_clock_time, min_total, max_total, &
                    total_calls, percent_total, percent_parent, load_imbalance, std_imbalance
            else
                ! Handle overflow for high-frequency calls
                write(unit, '(A35,3F15.4,I8,A8,A8,A10,A10)') &
                    display_name, wall_clock_time, min_total, max_total, &
                    total_calls, '   N/A ', '   N/A ', '     N/A  ', '     N/A  '
            endif
            
            ! Recursively print children with increased indentation
            call print_detailed_sections(unit, global_data, global_counts, npes, level+1, profiles(i)%name, total_runtime, wall_clock_time)
        end do
    end subroutine print_detailed_sections
    
    !=========================================================================
    ! Print sections in FESOM style with automatic hierarchy (simplified)
    !=========================================================================  
    recursive subroutine print_fesom_style_sections(unit, global_data, global_counts, npes, level, parent)
        integer, intent(in) :: unit, npes, level
        real(kind=PROF_WP), intent(in) :: global_data(:)
        integer, intent(in) :: global_counts(:)
        character(len=*), intent(in) :: parent
        integer :: i, j, total_calls
        real(kind=PROF_WP) :: wall_clock_time, min_total, max_total, load_imbalance
        character(len=100) :: display_name
        character(len=15) :: indent_prefix
        logical :: has_children
        
        ! Print sections at this level
        do i = 1, num_profiles
            if (trim(profiles(i)%name) == "") cycle
            ! Only print sections that belong at this level
            if (level == 0) then
                ! At top level, only print sections with no parent
                if (trim(profiles(i)%parent_name) /= "") cycle
            else
                ! At nested levels, only print sections whose parent matches
                if (trim(profiles(i)%parent_name) /= trim(parent)) cycle
            endif
            
            ! Extract statistics
            wall_clock_time = global_data(4*i-3) / real(npes, PROF_WP)
            min_total = global_data(4*i-2)
            max_total = global_data(4*i-1)
            total_calls = global_counts(i)
            
            if (total_calls == 0) cycle
            
            ! Calculate load imbalance
            if (wall_clock_time > 0.0_PROF_WP) then
                load_imbalance = (max_total - min_total) / wall_clock_time * 100.0_PROF_WP
            else
                load_imbalance = 0.0_PROF_WP
            endif
            
            ! Create indented display name
            indent_prefix = ""
            if (level == 0) then
                display_name = trim(profiles(i)%name)
            else
                write(indent_prefix, '(A)') repeat("  ", level)
                write(display_name, '(A,A)') trim(indent_prefix), trim(profiles(i)%name)
            endif
            
            ! Check if this section has children
            has_children = .false.
            do j = 1, num_profiles
                if (trim(profiles(j)%parent_name) == trim(profiles(i)%name)) then
                    has_children = .true.
                    exit
                endif
            end do
            
            ! Print in FESOM style: "name : time" with load balance info if significant
            if (load_imbalance > 15.0_PROF_WP) then
                write(unit, '(A30,A,1X,E10.3,A,F6.1,A)') &
                    display_name, ':', wall_clock_time, ' [LdImb:', load_imbalance, '%]'
            else
                write(unit, '(A30,A,1X,E10.3)') display_name, ':', wall_clock_time
            endif
            
            ! Recursively print children with increased indentation
            call print_fesom_style_sections(unit, global_data, global_counts, npes, level+1, profiles(i)%name)
            
            ! Print separator after top-level sections with children
            if (level == 0 .and. has_children) then
                write(unit, '(A)') repeat('_', 31)
            endif
        end do
    end subroutine print_fesom_style_sections
    
    !=========================================================================
    ! Print concise load balance summary
    !=========================================================================
    subroutine print_load_balance_summary(unit, global_data, global_counts, npes)
        integer, intent(in) :: unit, npes
        real(kind=PROF_WP), intent(in) :: global_data(:)
        integer, intent(in) :: global_counts(:)
        integer :: i, imbalanced_count, well_balanced_count, ignored_count, analyzed_count
        real(kind=PROF_WP) :: wall_clock_time, min_total, max_total, load_imbalance, worst_imbalance
        real(kind=PROF_WP) :: total_runtime, threshold_time
        character(len=100) :: worst_section
        logical :: found_total
        
        ! Calculate total runtime for meaningful threshold
        total_runtime = 0.0_PROF_WP
        found_total = .false.
        
        ! Search for main timing sections to determine total runtime
        do i = 1, num_profiles
            if (trim(profiles(i)%name) == "") cycle
            
            ! Look for main timing sections (same logic as benchmark summary)
            if (index(profiles(i)%name, "runloop_total") > 0 .or. &
                index(profiles(i)%name, "main_total") > 0 .or. &
                index(profiles(i)%name, "fesom_total") > 0) then
                wall_clock_time = global_data(4*i-3) / real(npes, PROF_WP)
                if (wall_clock_time > total_runtime) then
                    total_runtime = wall_clock_time
                    found_total = .true.
                endif
            endif
        end do
        
        ! If no main section found, sum up major sections
        if (.not. found_total) then
            do i = 1, num_profiles
                if (trim(profiles(i)%name) == "") cycle
                if (trim(profiles(i)%parent_name) == "") then  ! Only top-level sections
                    wall_clock_time = global_data(4*i-3) / real(npes, PROF_WP)
                    total_runtime = total_runtime + wall_clock_time
                endif
            end do
            if (total_runtime > 0.0_PROF_WP) found_total = .true.
        endif
        
        ! Set meaningful threshold: 1% of total runtime, minimum 0.01 seconds
        if (found_total .and. total_runtime > 0.0_PROF_WP) then
            threshold_time = max(0.01_PROF_WP, total_runtime * 0.01_PROF_WP)
        else
            threshold_time = 0.01_PROF_WP  ! Fallback to 0.01 seconds
        endif
        
        imbalanced_count = 0
        well_balanced_count = 0
        ignored_count = 0
        analyzed_count = 0
        worst_imbalance = 0.0_PROF_WP
        worst_section = ""
        
        do i = 1, num_profiles
            if (trim(profiles(i)%name) == "" .or. global_counts(i) == 0) cycle
            
            wall_clock_time = global_data(4*i-3) / real(npes, PROF_WP)
            min_total = global_data(4*i-2)
            max_total = global_data(4*i-1)
            
            if (wall_clock_time > threshold_time) then  ! Only analyze meaningful sections
                analyzed_count = analyzed_count + 1
                load_imbalance = (max_total - min_total) / wall_clock_time * 100.0_PROF_WP
                
                if (load_imbalance > worst_imbalance) then
                    worst_imbalance = load_imbalance
                    worst_section = trim(profiles(i)%name)
                endif
                
                if (load_imbalance > 15.0_PROF_WP) then
                    imbalanced_count = imbalanced_count + 1
                else if (load_imbalance < 5.0_PROF_WP) then
                    well_balanced_count = well_balanced_count + 1
                endif
            else
                ignored_count = ignored_count + 1
            endif
        end do
        
        write(unit, '(A,I0,A,I0,A)') '  Analyzed sections: ', analyzed_count, ' (ignored ', ignored_count, ' minor sections)'
        write(unit, '(A,I0)') '  Well-balanced sections: ', well_balanced_count
        write(unit, '(A,I0)') '  Imbalanced sections:    ', imbalanced_count
        if (worst_imbalance > 0.0_PROF_WP) then
            write(unit, '(A,F6.1,A,A,A)') '  Worst imbalance: ', worst_imbalance, '% (', trim(worst_section), ')'
        endif
        
        if (ignored_count > 0) then
            write(unit, '(A,F8.4,A)') '  Note: Sections < ', threshold_time, 's (1% of total runtime) ignored (init/finalize/minor overhead)'
        endif
    end subroutine print_load_balance_summary

    !=========================================================================
    ! Print load balance analysis
    !=========================================================================
    subroutine print_load_balance_analysis(unit, global_data, global_counts, npes)
        integer, intent(in) :: unit, npes
        real(kind=PROF_WP), intent(in) :: global_data(:)
        integer, intent(in) :: global_counts(:)
        integer :: i, imbalanced_sections, well_balanced_sections
        real(kind=PROF_WP) :: total_time, min_time, max_time, mean_time, load_imbalance
        real(kind=PROF_WP) :: worst_imbalance
        character(len=100) :: worst_section
        
        write(unit, '(A)') 'LOAD BALANCE ANALYSIS:'
        write(unit, '(A)') repeat('-', 40)
        
        imbalanced_sections = 0
        well_balanced_sections = 0
        worst_imbalance = 0.0_PROF_WP
        worst_section = ""
        
        do i = 1, num_profiles
            if (trim(profiles(i)%name) == "" .or. global_counts(i) == 0) cycle
            
            total_time = global_data(4*i-3)
            min_time = global_data(4*i-2)
            max_time = global_data(4*i-1)
            mean_time = total_time / real(npes, PROF_WP)
            
            if (mean_time > 0.001_PROF_WP) then  ! Only analyze sections with significant time
                load_imbalance = (max_time - min_time) / mean_time * 100.0_PROF_WP
                
                if (load_imbalance > worst_imbalance) then
                    worst_imbalance = load_imbalance
                    worst_section = trim(profiles(i)%name)
                endif
                
                if (load_imbalance > 15.0_PROF_WP) then
                    imbalanced_sections = imbalanced_sections + 1
                else if (load_imbalance < 5.0_PROF_WP) then
                    well_balanced_sections = well_balanced_sections + 1
                endif
            endif
        end do
        
        write(unit, '(A,I0)') '  Well-balanced sections (< 5% imbalance): ', well_balanced_sections
        write(unit, '(A,I0)') '  Imbalanced sections (> 15% imbalance):   ', imbalanced_sections
        write(unit, '(A,F6.1,A)') '  Worst load imbalance: ', worst_imbalance, '%'
        if (worst_imbalance > 0.0_PROF_WP) then
            write(unit, '(A,A)') '  Worst section: ', trim(worst_section)
        endif
        write(unit, '(A)') ''
        
        ! Recommendations
        if (imbalanced_sections > 0) then
            write(unit, '(A)') 'RECOMMENDATIONS:'
            write(unit, '(A)') '  - Consider load balancing optimization for imbalanced sections'
            write(unit, '(A)') '  - Check for uneven mesh distribution or boundary effects'
            write(unit, '(A)') '  - Consider asynchronous communication patterns'
        else
            write(unit, '(A)') 'GOOD NEWS: All major sections are well-balanced!'
        endif
        write(unit, '(A)') ''
    end subroutine print_load_balance_analysis

    !=========================================================================
    ! Print benchmark summary in FESOM style
    !=========================================================================
    subroutine print_benchmark_summary(unit, global_data, npes, total_runtime)
        integer, intent(in) :: unit, npes
        real(kind=PROF_WP), intent(in) :: global_data(:)
        real(kind=PROF_WP), intent(in) :: total_runtime
        real(kind=PROF_WP) :: sypd, simulated_time, simulated_years
        
        write(unit, '(A)') repeat('=', 90)
        write(unit, '(A)') '============================ BENCHMARK RUNTIME ============================'
        write(unit, '(A,I0)') '    Number of cores :                      ', npes
        write(unit, '(A,I0)') '    Number of OMP threads per rank :      ', omp_threads
        if (total_timesteps > 0) then
            write(unit, '(A,I0)') '    Number of timesteps :                  ', total_timesteps
        endif
        write(unit, '(A,F12.4,A)') '    Runtime for all timesteps :      ', total_runtime, ' sec'
        
        ! Calculate SYPD (Simulated Years Per Day) if we have timestep info
        if (total_timesteps > 0 .and. timestep_size > 0.0_PROF_WP .and. total_runtime > 0.0_PROF_WP) then
            simulated_time = timestep_size * real(total_timesteps, PROF_WP)  ! seconds
            simulated_years = simulated_time / (365.25_PROF_WP * 86400.0_PROF_WP)  ! years (accounting for leap years)
            sypd = simulated_years / (total_runtime / 86400.0_PROF_WP)  ! years per day
            write(unit, '(A,F12.4,A)') '    Estimated SYPD :                  ', sypd, ' years/day'
        endif
        
        write(unit, '(A)') repeat('=', 90)
    end subroutine print_benchmark_summary

    !=========================================================================
    ! Helper function to find total runtime and call benchmark summary
    !=========================================================================
    subroutine print_benchmark_summary_with_total_runtime(unit, global_data, global_counts, npes)
        integer, intent(in) :: unit, npes
        real(kind=PROF_WP), intent(in) :: global_data(:)
        integer, intent(in) :: global_counts(:)
        integer :: i
        real(kind=PROF_WP) :: total_runtime, wall_clock_time
        logical :: found_total
        
        ! Look for total runtime from main sections
        total_runtime = 0.0_PROF_WP
        found_total = .false.
        
        ! Search for common main section names
        do i = 1, num_profiles
            if (trim(profiles(i)%name) == "") cycle
            
            ! Look for main timing sections
            if (index(profiles(i)%name, "runloop_total") > 0 .or. &
                index(profiles(i)%name, "main_total") > 0 .or. &
                index(profiles(i)%name, "fesom_total") > 0) then
                wall_clock_time = global_data(4*i-3) / real(npes, PROF_WP)
                if (wall_clock_time > total_runtime) then
                    total_runtime = wall_clock_time
                    found_total = .true.
                endif
            endif
        end do
        
        ! If no main section found, sum up major sections
        if (.not. found_total) then
            do i = 1, num_profiles
                if (trim(profiles(i)%name) == "") cycle
                if (trim(profiles(i)%parent_name) == "") then  ! Only top-level sections
                    wall_clock_time = global_data(4*i-3) / real(npes, PROF_WP)
                    total_runtime = total_runtime + wall_clock_time
                endif
            end do
            if (total_runtime > 0.0_PROF_WP) found_total = .true.
        endif
        
        ! Print benchmark summary if we found a reasonable total runtime
        if (found_total .and. total_runtime > 0.001_PROF_WP) then
            call print_benchmark_summary(unit, global_data, npes, total_runtime)
        endif
    end subroutine print_benchmark_summary_with_total_runtime

    !=========================================================================
    ! Finalize profiler and generate final report
    !=========================================================================
    subroutine fesom_profiler_finalize(mpi_comm, mpi_rank)
        include "mpif.h"
        integer, intent(in) :: mpi_comm, mpi_rank
        
        if (.not. fesom_profiler_enabled()) return
        
        call fesom_profiler_report(mpi_comm, mpi_rank)
        profiler_initialized = .false.
    end subroutine fesom_profiler_finalize

    !=========================================================================
    ! Reset all profiling data
    !=========================================================================
    subroutine fesom_profiler_reset()
        integer :: i
        
        if (.not. fesom_profiler_enabled()) return
        
        do i = 1, num_profiles
            profiles(i)%total_time = 0.0_PROF_WP
            profiles(i)%min_time = HUGE(1.0_PROF_WP)
            profiles(i)%max_time = 0.0_PROF_WP
            profiles(i)%sum_squares = 0.0_PROF_WP
            profiles(i)%call_count = 0
            profiles(i)%is_active = .false.
        end do
        
        call_stack_depth = 0
    end subroutine fesom_profiler_reset

    !=========================================================================
    ! Find existing profile or create new one
    !=========================================================================
    function find_or_create_profile(name) result(index)
        character(len=*), intent(in) :: name
        integer :: index
        
        ! First try to find existing profile
        index = find_profile(name)
        if (index /= -1) return
        
        ! Create new profile
        if (num_profiles < MAX_PROFILE_SECTIONS) then
            num_profiles = num_profiles + 1
            index = num_profiles
            profiles(index)%name = trim(name)
            profiles(index)%total_time = 0.0_PROF_WP
            profiles(index)%min_time = HUGE(1.0_PROF_WP)
            profiles(index)%max_time = 0.0_PROF_WP
            profiles(index)%sum_squares = 0.0_PROF_WP
            profiles(index)%call_count = 0
            profiles(index)%is_active = .false.
        else
            ! Error: too many profiles
            index = -1
        endif
    end function find_or_create_profile

    !=========================================================================
    ! Find existing profile by name
    !=========================================================================
    function find_profile(name) result(index)
        character(len=*), intent(in) :: name
        integer :: index
        integer :: i
        
        index = -1
        do i = 1, num_profiles
            if (trim(profiles(i)%name) == trim(name)) then
                index = i
                return
            endif
        end do
    end function find_profile

end module fesom_profiler