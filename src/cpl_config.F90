module cpl_config
  !============================================================================
  ! Runtime configuration of FESOM's coupling to external components.
  !
  ! The interface (standalone / direct / oasis28 / oasis50 / yac) is a
  ! compile-time choice, selected by one __cpl_* macro; see
  ! cmake/FesomCoupling.cmake. The external components FESOM is coupled to
  ! are configured here at run time, from config/namelist.cpl.
  !
  ! read_cpl_namelist runs before cpl_oasis3mct_init / cpl_yac_init, which
  ! both need values from it. Under FESOM_COUPLING=oasis50 that is before
  ! MPI_Init, since OASIS performs it. It therefore uses plain Fortran I/O,
  ! "error stop" rather than par_ex or MPI_Abort, and produces no output;
  ! check_cpl_config reports once a rank number is available.
  !============================================================================
  implicit none
  save
  private

  public :: read_cpl_namelist, check_cpl_config
  public :: is_coupled_to_echam, is_coupled_to_oifs
  public :: is_coupled_to_icon_a, is_coupled_to_ifs
  public :: cpl_comp_name, cpl_grid_name, cpl_config_file
  public :: compute_oasis_corners

  !____________________________________________________________________________
  ! &coupling_partner -- the external components FESOM is coupled to.
  !
  ! Interface-independent, and checked against the compiled interface by
  ! check_cpl_config. The default is the partner implied by the compiled
  ! interface, so an existing run directory without a namelist.cpl keeps
  ! working unchanged.
  !
  ! More than one may be true once a second component role (sea ice, ice
  ! sheet) exists; at most one atmosphere.
#if defined (__cpl_oasis28)
  logical :: is_coupled_to_echam  = .true.
#else
  logical :: is_coupled_to_echam  = .false.
#endif
#if defined (__cpl_oasis50)
  logical :: is_coupled_to_oifs   = .true.
#else
  logical :: is_coupled_to_oifs   = .false.
#endif
#if defined (__cpl_yac)
  logical :: is_coupled_to_icon_a = .true.
#else
  logical :: is_coupled_to_icon_a = .false.
#endif
#if defined (__cpl_direct)
  logical :: is_coupled_to_ifs    = .true.
#else
  logical :: is_coupled_to_ifs    = .false.
#endif

  namelist /coupling_partner/ is_coupled_to_echam, is_coupled_to_oifs, &
                              is_coupled_to_icon_a, is_coupled_to_ifs

  !____________________________________________________________________________
  ! Interface-specific settings. The two groups are mutually exclusive -- only
  ! the one matching the compiled interface is read -- so a single pair of
  ! variables backs both, with interface-specific defaults.
  character(len=32)  :: cpl_comp_name
  character(len=32)  :: cpl_grid_name
  character(len=256) :: cpl_config_file = 'coupling.yaml'

  ! Write grid corners to the OASIS grid files. Needed when the namcouple asks
  ! for first-order conservative remapping; not a property of the partner, so
  ! it stays a user option.
  logical :: compute_oasis_corners = .false.

#if defined (__cpl_yac)
  data cpl_comp_name /'fesom2'/
  data cpl_grid_name /'fesom_grid'/
#else
  data cpl_comp_name /'fesom'/
  data cpl_grid_name /'feom'/
#endif

  namelist /coupling_oasis/ cpl_comp_name, cpl_grid_name, &
                            compute_oasis_corners
  namelist /coupling_yac/   cpl_comp_name, cpl_grid_name, cpl_config_file

  character(len=*), parameter :: nmlfile = 'namelist.cpl'

  !> Set by read_cpl_namelist when the file is absent, so that
  !> check_cpl_config can say so once rank 0 is known.
  logical :: nmlfile_missing = .false.

contains

  !____________________________________________________________________________
  !> Read config/namelist.cpl. Called at the very start of fesom_init, before
  !> MPI_Init and before the coupler is initialised. Not called at all in a
  !> standalone build.
  !>
  !> Like the other FESOM namelists this is read by every rank from a bare
  !> relative path, with no broadcast; a missing file leaves the compiled-in
  !> defaults in place. What was read is reported by check_cpl_config, once a
  !> rank number is available.
  subroutine read_cpl_namelist()
    integer :: fileunit, istat
    logical :: file_exists

    inquire(file=nmlfile, exist=file_exists)
    if (.not. file_exists) then
       nmlfile_missing = .true.
       return
    end if

    open(newunit=fileunit, file=nmlfile, status='OLD', iostat=istat)
    if (istat /= 0) error stop 'could not open namelist.cpl'

    read(fileunit, nml=coupling_partner, iostat=istat)
    if (istat /= 0) error stop &
       'could not read &coupling_partner from namelist.cpl'

    ! Only the group belonging to the compiled interface is read. Both groups
    ! are optional: an existing file that predates them keeps working.
#if defined (__cpl_oasis)
    read(fileunit, nml=coupling_oasis, iostat=istat)
    if (istat /= 0) rewind(fileunit)
#elif defined (__cpl_yac)
    read(fileunit, nml=coupling_yac, iostat=istat)
    if (istat /= 0) rewind(fileunit)
#endif

    close(fileunit)
  end subroutine read_cpl_namelist

  !____________________________________________________________________________
  !> Validate the partner selection against the compiled interface. Called
  !> from fesom_init after setup_model, so that the deprecated
  !> compute_oasis_corners in namelist.config's &run_config can be folded in.
  subroutine check_cpl_config(deprecated_oasis_corners, comm, mype)
    logical, intent(in) :: deprecated_oasis_corners
    integer, intent(in) :: comm
    integer, intent(in) :: mype

    integer :: n_atm

    if (mype == 0 .and. nmlfile_missing) then
       write(*,*) '     could not find ', nmlfile, &
                  ', will use default values !'
    end if
    call report_cpl_config(mype)

    !__________________________________________________________________________
    ! compute_oasis_corners used to live in &run_config of namelist.config.
    ! Either location still switches it on.
    if (deprecated_oasis_corners) then
       if (mype == 0 .and. .not. compute_oasis_corners) then
          write(*,*) 'WARNING: compute_oasis_corners in &run_config of ', &
                     'namelist.config is deprecated, move it to ', &
                     '&coupling_oasis in ', nmlfile
       end if
       compute_oasis_corners = .true.
    end if

    !__________________________________________________________________________
    ! Exactly one atmosphere.
    n_atm = 0
    if (is_coupled_to_echam)  n_atm = n_atm + 1
    if (is_coupled_to_oifs)   n_atm = n_atm + 1
    if (is_coupled_to_icon_a) n_atm = n_atm + 1
    if (is_coupled_to_ifs)    n_atm = n_atm + 1

#if defined (__cpl_enabled)
    if (n_atm == 0) then
       call cpl_config_abort(comm, mype, &
            'no external atmosphere selected. Set exactly one of '// &
            'is_coupled_to_echam / _oifs / _icon_a / _ifs in '// &
            '&coupling_partner of '//nmlfile)
    end if
    if (n_atm > 1) then
       call cpl_config_abort(comm, mype, &
            'more than one external atmosphere selected in '// &
            '&coupling_partner of '//nmlfile// &
            '. FESOM couples to at most one atmosphere.')
    end if
#else
    if (n_atm > 0) then
       call cpl_config_abort(comm, mype, &
            'an external atmosphere is selected in '//nmlfile// &
            ', but this is a standalone build (FESOM_COUPLING=standalone).')
    end if
#endif

    !__________________________________________________________________________
    ! Each partner needs the interface it talks over.
#if !defined (__cpl_oasis28)
    if (is_coupled_to_echam) &
       call wrong_interface(comm, mype, 'echam', 'oasis28')
#endif
#if !defined (__cpl_oasis50)
    if (is_coupled_to_oifs) &
       call wrong_interface(comm, mype, 'oifs', 'oasis50')
#endif
#if !defined (__cpl_yac)
    if (is_coupled_to_icon_a) &
       call wrong_interface(comm, mype, 'icon_a', 'yac')
#endif
#if !defined (__cpl_direct)
    if (is_coupled_to_ifs) &
       call wrong_interface(comm, mype, 'ifs', 'direct')
#endif

    !__________________________________________________________________________
    ! Settings that only mean something for one interface.
#if !defined (__cpl_oasis)
    if (compute_oasis_corners .and. mype == 0) then
       write(*,*) 'WARNING: compute_oasis_corners is ignored, this build ', &
                  'does not use OASIS'
    end if
#endif
  end subroutine check_cpl_config

  !____________________________________________________________________________
  subroutine wrong_interface(comm, mype, partner, interface_name)
    integer,          intent(in) :: comm
    integer,          intent(in) :: mype
    character(len=*), intent(in) :: partner
    character(len=*), intent(in) :: interface_name

    call cpl_config_abort(comm, mype, &
         'is_coupled_to_'//trim(partner)//' requires a build configured '// &
         'with -DFESOM_COUPLING='//trim(interface_name))
  end subroutine wrong_interface

  !____________________________________________________________________________
  subroutine report_cpl_config(mype)
    integer, intent(in) :: mype

    if (mype /= 0) return
    write(*,*) '     coupling configuration:'
    write(*,*) '        is_coupled_to_echam  = ', is_coupled_to_echam
    write(*,*) '        is_coupled_to_oifs   = ', is_coupled_to_oifs
    write(*,*) '        is_coupled_to_icon_a = ', is_coupled_to_icon_a
    write(*,*) '        is_coupled_to_ifs    = ', is_coupled_to_ifs
    write(*,*) '        cpl_comp_name        = ', trim(cpl_comp_name)
    write(*,*) '        cpl_grid_name        = ', trim(cpl_grid_name)
#if defined (__cpl_oasis)
    write(*,*) '        compute_oasis_corners= ', compute_oasis_corners
#elif defined (__cpl_yac)
    write(*,*) '        cpl_config_file      = ', trim(cpl_config_file)
#endif
  end subroutine report_cpl_config

  !____________________________________________________________________________
  subroutine cpl_config_abort(comm, mype, message)
    use iso_fortran_env, only: error_unit
    use mpi
    integer,          intent(in) :: comm
    integer,          intent(in) :: mype
    character(len=*), intent(in) :: message

    integer :: ierror

    if (mype == 0) then
       write(error_unit,'(A)') 'ERROR: '//trim(message)
    end if
    flush(error_unit)
    call MPI_Abort(comm, 1, ierror)
  end subroutine cpl_config_abort

end module cpl_config
