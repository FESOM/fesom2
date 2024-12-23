program test_ifs_interface
    use mpi
    use par_kind, only: wpIFS, JPTIME
    implicit none

    ! MPI variables
    integer :: mype, npes, icomm, ierr
    
    ! Variables needed for nemogcmcoup_init
    integer :: inidate, initime, itini, itend
    real(wpIFS) :: zstp
    logical :: lwaveonly = .false.
    integer :: iatmunit = -1  ! Not used in FESOM
    logical :: lwrite = .false.

    ! Initialize MPI (simulating what IFS would do)
    call MPI_Init(ierr)
    if (ierr /= 0) then
        write(*,*) 'Error initializing MPI'
        stop 1
    endif

    ! Get basic MPI info
    call MPI_Comm_rank(MPI_COMM_WORLD, mype, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, npes, ierr)
    icomm = MPI_COMM_WORLD

    if (mype == 0) then
        write(*,*) '=== Testing FESOM initialization from IFS ==='
        write(*,*) 'Running on ', npes, ' processes'
    endif

    ! Create namfesomstep.in for testing
    !if (mype == 0) then
    !    open(9, file='namfesomstep.in', status='replace')
    !    write(9,'(A)') ' &namfesomstep'
    !    write(9,'(A)') '   substeps=4'
    !    write(9,'(A)') ' /'
    !    close(9)
    !endif
    call MPI_Barrier(icomm, ierr)

    ! Call the initialization routine
    call nemogcmcoup_init(mype, icomm, inidate, initime, itini, itend, zstp, &
                         lwaveonly, iatmunit, lwrite)

    if (mype == 0) then
        write(*,*) '=== Initialization Results ==='
        write(*,*) 'Initial date: ', inidate
        write(*,*) 'Initial time: ', initime
        write(*,*) 'Initial timestep: ', itini
        write(*,*) 'Final timestep: ', itend
        write(*,*) 'Timestep length (s): ', zstp
    endif

    ! Cleanup and finalize
    call nemogcmcoup_final()
    call MPI_Finalize(ierr)

end program test_ifs_interface 
