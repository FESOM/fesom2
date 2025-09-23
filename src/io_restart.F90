MODULE io_RESTART
  use restart_file_group_module
  use restart_derivedtype_module
  use g_clock
  use g_config
  use o_arrays
  use g_cvmix_tke
  use g_cvmix_idemix
  use g_backscatter
  use MOD_TRACER
  use MOD_ICE
  use MOD_DYN
  use MOD_MESH
  use MOD_PARTIT
  use MOD_PARSUP
  use fortran_utils
  use mpi
#if defined(__icepack)
  use icedrv_main
#endif 
#if defined(__recom)
  use recom_glovar
  use recom_config
  use recom_ciso
#endif
  
  implicit none
  public :: restart, finalize_restart
  private

  integer,       save       :: globalstep=0 ! todo: remove this from module scope as it will mess things up if we use async read/write from the same process
  real(kind=WP)             :: ctime !current time in seconds from the beginning of the year

  type(restart_file_group) , save :: oce_files
  character(:), allocatable, save :: oce_path
  
  type(restart_file_group) , save :: ice_files
  character(:), allocatable, save :: ice_path
  
#if defined(__icepack)
  type(restart_file_group) , save, public :: icepack_files
  character(:), allocatable, save, public :: icepack_path
#endif  

#if defined(__recom)
  type(restart_file_group) , save :: bio_files
  character(:), allocatable, save :: bio_path
#endif

  character(:), allocatable, save :: raw_restart_dirpath
  character(:), allocatable, save :: raw_restart_infopath
  character(:), allocatable, save :: bin_restart_dirpath
  character(:), allocatable, save :: bin_restart_infopath
  integer, parameter :: RAW_RESTART_METADATA_RANK = 0


  contains

!--------------------------------------------------------------------------------------------
! Helper functions for constructing restart file paths
!--------------------------------------------------------------------------------------------

! Build NetCDF restart file path
pure function nc_restart_path(component, year, root_path) result(path)
  implicit none
  character(len=*), intent(in) :: component, root_path
  integer, intent(in) :: year
  character(:), allocatable :: path
  character(4) :: cyear
  
  write(cyear, '(i4)') year
  path = trim(root_path) // trim(runid) // '.' // cyear // '.' // trim(component) // '.restart.nc'
end function nc_restart_path

! Build raw restart directory path
pure function build_raw_restart_dirpath(root_path) result(path)
  implicit none
  character(len=*), intent(in) :: root_path
  character(:), allocatable :: path
  
  path = trim(root_path) // trim(runid) // '_raw_restart'
end function build_raw_restart_dirpath

! Build raw restart info file path
pure function build_raw_restart_infopath(root_path) result(path)
  implicit none
  character(len=*), intent(in) :: root_path
  character(:), allocatable :: path
  
  path = trim(root_path) // trim(runid) // '_raw_restart.info'
end function build_raw_restart_infopath

! Build binary restart directory path
pure function build_bin_restart_dirpath(root_path) result(path)
  implicit none
  character(len=*), intent(in) :: root_path
  character(:), allocatable :: path
  
  path = trim(root_path) // trim(runid) // '_bin_restart'
end function build_bin_restart_dirpath

! Build binary restart info file path
pure function build_bin_restart_infopath(root_path) result(path)
  implicit none
  character(len=*), intent(in) :: root_path
  character(:), allocatable :: path
  
  path = trim(root_path) // trim(runid) // '_bin_restart.info'
end function build_bin_restart_infopath

!
!--------------------------------------------------------------------------------------------
! ini_ocean_io initializes ocean_file datatype which contains information of all variables need to be written into 
! the ocean restart file. This is the only place need to be modified if a new variable is added!
subroutine ini_ocean_io(year, dynamics, tracers, partit, mesh)
#ifdef ENABLE_NVHPC_WORKAROUNDS
  use nvfortran_subarray_workaround_module
#endif
  integer, intent(in)       :: year
  integer                   :: j, id
  character(500)            :: longname
  character(500)            :: trname, units
  character(4)              :: cyear
  type(t_mesh), target :: mesh
  type(t_partit), intent(inout), target :: partit
  type(t_tracer), target :: tracers
  type(t_dyn), target :: dynamics
  logical, save :: has_been_called = .false.

  write(cyear,'(i4)') year
  oce_path = nc_restart_path('oce', year, RestartOutPath)

  if(has_been_called) return
  has_been_called = .true.

  !===========================================================================
  !===================== Definition part =====================================
  !===========================================================================
  !___Define the netCDF variables for 2D fields_______________________________
  !___SSH_____________________________________________________________________
  call oce_files%def_node_var('ssh', 'sea surface elevation', 'm',   dynamics%eta_n, mesh, partit)
  !___ALE related fields______________________________________________________
  call oce_files%def_node_var('hbar', 'ALE surface elevation', 'm',   mesh%hbar, mesh, partit)
!!PS   call oce_files%def_node_var('ssh_rhs', 'RHS for the elevation', '?',   ssh_rhs, mesh, partit)
  call oce_files%def_node_var('ssh_rhs_old', 'RHS for the elevation', '?',   dynamics%ssh_rhs_old, mesh, partit)
  call oce_files%def_node_var('hnode', 'nodal layer thickness', 'm',   mesh%hnode, mesh, partit)
  
  !___Define the netCDF variables for 3D fields_______________________________
#ifdef ENABLE_NVHPC_WORKAROUNDS
  dynamics_workaround => dynamics
#endif
  call oce_files%def_elem_var('u', 'zonal velocity',        'm/s', dynamics%uv(1,:,:), mesh, partit)
  call oce_files%def_elem_var('v', 'meridional velocity',   'm/s', dynamics%uv(2,:,:), mesh, partit)
  call oce_files%def_elem_var('urhs_AB', 'Adams-Bashforth for u (n-1 for AB2 and n-2 for AB3)', 'm/s', dynamics%uv_rhsAB(1,1,:,:), mesh, partit)
  call oce_files%def_elem_var('vrhs_AB', 'Adams-Bashforth for v (n-1 for AB2 and n-2 for AB3)', 'm/s', dynamics%uv_rhsAB(1,2,:,:), mesh, partit)
  if (dynamics%AB_order==3) then
       call oce_files%def_elem_var_optional('urhs_AB3', 'Adams-Bashforth for u (n-1) for AB3', 'm/s', dynamics%uv_rhsAB(2,1,:,:), mesh, partit)
       call oce_files%def_elem_var_optional('vrhs_AB3', 'Adams-Bashforth for v (n-1) for AB3', 'm/s', dynamics%uv_rhsAB(2,2,:,:), mesh, partit)
  end if
  
!___Save restart variables for TKE and IDEMIX_________________________________
!   if (trim(mix_scheme)=='cvmix_TKE' .or. trim(mix_scheme)=='cvmix_TKE+IDEMIX') then
  if (mix_scheme_nmb==5 .or. mix_scheme_nmb==56) then
        call oce_files%def_node_var_optional('tke', 'Turbulent Kinetic Energy', 'm2/s2', tke(:,:), mesh, partit)
  endif
!   if (trim(mix_scheme)=='cvmix_IDEMIX' .or. trim(mix_scheme)=='cvmix_TKE+IDEMIX') then
  if (mix_scheme_nmb==6 .or. mix_scheme_nmb==56) then
        call oce_files%def_elem_var_optional('iwe', 'Internal Wave Energy'    , 'm2/s2', iwe(:,:), mesh, partit)
  endif 
  if (dynamics%opt_visc==8) then
        call oce_files%def_elem_var_optional('uke', 'unresolved kinetic energy', 'm2/s2', uke(:,:), mesh, partit)
        call oce_files%def_elem_var_optional('uke_rhs', 'unresolved kinetic energy rhs', 'm2/s2', uke_rhs(:,:), mesh, partit)
  endif
  
  do j=1,tracers%num_tracers
     id=tracers%data(j)%ID  !MB: Avoid hard-wired tracer assignments like SELECT CASE(j)
     SELECT CASE (id) 
       CASE(1)
         trname='temp'
         longname='potential temperature'
         units='degC'
       CASE(2)
         trname='salt'
         longname='salinity'
         units='psu'
       CASE(6)
         trname='sf6'
         longname='sulfur hexafluoride'
         units='mol / m**3'
       CASE(11)
         trname='cfc11'
         longname='chlorofluorocarbon CFC-11'
         units='mol / m**3'
       CASE(12)
         trname='cfc12'
         longname='chlorofluorocarbon CFC-12'
         units='mol / m**3'
       CASE(14)
         trname='r14c'
         longname='14C / C ratio of DIC'
         units='none'
       CASE(39)
         trname='r39ar'
         longname='39Ar / Ar ratio'
         units='none'
       CASE(101)
         trname='h2o18'
         longname='h2o18 concentration'
         units='kmol/m**3'
       CASE(102)
         trname='hDo16'
         longname='hDo16 concentration'
         units='kmol/m**3'
       CASE(103)
         trname='h2o16'
         longname='h2o16 concentration'
         units='kmol/m**3'
       CASE DEFAULT
         write(trname,'(A3,i4.4)') 'tra_', j
         write(longname,'(A15,i4.4)') 'passive tracer ', j
         units='none'
     END SELECT
     if ((tracers%data(j)%ID==101) .or. (tracers%data(j)%ID==102) .or. (tracers%data(j)%ID==103)) then
        call oce_files%def_node_var_optional(trim(trname), trim(longname), trim(units), tracers%data(j)%values(:,:), mesh, partit)
     else
        call oce_files%def_node_var(trim(trname), trim(longname), trim(units), tracers%data(j)%values(:,:), mesh, partit)
     endif
     longname=trim(longname)//', Adams-Bashforth'
     if ((tracers%data(j)%ID==101) .or. (tracers%data(j)%ID==102) .or. (tracers%data(j)%ID==103)) then
        call oce_files%def_node_var_optional(trim(trname)//'_AB', trim(longname), trim(units), tracers%data(j)%valuesAB(:,:),    mesh, partit)
     else
        call oce_files%def_node_var(trim(trname)//'_AB', trim(longname), trim(units), tracers%data(j)%valuesAB(:,:),    mesh, partit)
     endif
     call oce_files%def_node_var_optional(trim(trname)//'_M1', trim(longname), trim(units), tracers%data(j)%valuesold(1,:,:), mesh, partit)
     if (tracers%data(j)%AB_order==3) &
     call oce_files%def_node_var_optional(trim(trname)//'_M2', trim(longname), trim(units), tracers%data(j)%valuesold(2,:,:), mesh, partit)
  end do
  call oce_files%def_node_var('w', 'vertical velocity', 'm/s', dynamics%w, mesh, partit)
  call oce_files%def_node_var('w_expl', 'vertical velocity', 'm/s', dynamics%w_e, mesh, partit)
  call oce_files%def_node_var('w_impl', 'vertical velocity', 'm/s', dynamics%w_i, mesh, partit)
end subroutine ini_ocean_io
!
!--------------------------------------------------------------------------------------------
! ini_ice_io initializes ice_file datatype which contains information of all variables need to be written into
! the ice restart file. This is the only place need to be modified if a new variable is added!
subroutine ini_ice_io(year, ice, partit, mesh)
  integer,      intent(in)  :: year
  character(4)              :: cyear
  type(t_mesh), intent(in) , target :: mesh
  type(t_partit), intent(inout), target :: partit
  type(t_ice), target :: ice
  logical, save :: has_been_called = .false.

  write(cyear,'(i4)') year
  ice_path = nc_restart_path('ice', year, RestartOutPath)

  if(has_been_called) return
  has_been_called = .true.

  !===========================================================================
  !===================== Definition part =====================================
  !===========================================================================
  !___Define the netCDF variables for 2D fields_______________________________
  call ice_files%def_node_var('area', 'ice concentration [0 to 1]', '%',   ice%data(1)%values(:), mesh, partit)
  call ice_files%def_node_var('hice', 'effective ice thickness',    'm',   ice%data(2)%values(:), mesh, partit)
  call ice_files%def_node_var('hsnow', 'effective snow thickness',  'm',   ice%data(3)%values(:), mesh, partit)
  call ice_files%def_node_var('uice', 'zonal velocity',             'm/s', ice%uice, mesh, partit)
  call ice_files%def_node_var('vice', 'meridional velocity',        'm',   ice%vice, mesh, partit)
#if defined (__oifs)
  call ice_files%def_node_var_optional('ice_albedo', 'ice albedo',    '-',   ice%atmcoupl%ice_alb, mesh, partit)
  call ice_files%def_node_var_optional('ice_temp', 'ice surface temperature',  'K',   ice%data(4)%values, mesh, partit)
#endif /* (__oifs) */
#if defined (__oasis)
  !---wiso-code
  if (lwiso) then
    call ice_files%def_node_var_optional('h2o18_ice', 'h2o18 concentration in sea ice', 'kmol/m**3', tr_arr_ice(:,1), mesh, partit)
    call ice_files%def_node_var_optional('hDo16_ice', 'hDo16 concentration in sea ice', 'kmol/m**3', tr_arr_ice(:,2), mesh, partit)
    call ice_files%def_node_var_optional('h2o16_ice', 'h2o16 concentration in sea ice', 'kmol/m**3', tr_arr_ice(:,3), mesh, partit)
  end if
  !---wiso-code-end
#endif

end subroutine ini_ice_io
!
!--------------------------------------------------------------------------------------------
!

! ini_bio_io initializes bid datatype which contains information of all variables need to be written into
! the bio restart file. This is the only place need to be modified if a new variable is added!
#if defined(__recom)
subroutine ini_bio_io(year, tracers, partit, mesh)
    integer,      intent(in)  :: year
    integer                   :: j
    character(500)            :: longname
    character(500)            :: trname, units
    character(4)              :: cyear

    type(t_mesh), intent(in) , target :: mesh
    type(t_partit), intent(inout), target :: partit
    type(t_tracer), target :: tracers
    logical, save :: has_been_called = .false.

    write(cyear,'(i4)') year
    bio_path = nc_restart_path('bio', year, RestartOutPath)

    if(has_been_called) return
    has_been_called = .true.

    !===========================================================================
    !===================== Definition part =====================================
    !===========================================================================
    !___Define the netCDF variables for 2D fields_______________________________
    call bio_files%def_node_var('BenN',    'Benthos Nitrogen', 'mmol/m3',   Benthos(:,1), mesh, partit);
    call bio_files%def_node_var('BenC',    'Benthos Carbon',   'mmol/m3',   Benthos(:,2), mesh, partit);
    call bio_files%def_node_var('BenSi',   'Benthos Silicate', 'mmol/m3',   Benthos(:,3), mesh, partit);
    call bio_files%def_node_var('BenCalc', 'Benthos Calcite',  'mmol/m3',   Benthos(:,4), mesh, partit);
    call bio_files%def_node_var('HPlus',   'Conc. of H-plus ions in the surface water', 'mol/kg',   GloHplus, mesh, partit);

end subroutine ini_bio_io
#endif

!--------------------------------------------------------------------------------------------
! Separate subroutine for reading restart files (initial conditions)
subroutine read_initial_conditions(which_readr, ice, dynamics, tracers, partit, mesh)
  use fortran_utils
  implicit none
  
  ! Parameters
  type(t_mesh)  , intent(inout), target :: mesh
  type(t_partit), intent(inout), target :: partit
  type(t_tracer), intent(inout), target :: tracers
  type(t_dyn)   , intent(inout), target :: dynamics
  type(t_ice)   , intent(inout), target :: ice
  integer, intent(out) :: which_readr
  
  ! Local variables
  logical :: rawfiles_exist, binfiles_exist
  integer :: mpierr
  character(:), allocatable :: read_raw_dirpath, read_raw_infopath
  character(:), allocatable :: read_bin_dirpath, read_bin_infopath
  character(:), allocatable :: read_oce_path, read_ice_path, read_bio_path
  
  ! Build paths for reading using RestartInPath
  read_raw_dirpath = build_raw_restart_dirpath(RestartInPath)//"/np"//int_to_txt(partit%npes)
  read_raw_infopath = build_raw_restart_infopath(RestartInPath)//"/np"//int_to_txt(partit%npes)//".info"
  read_bin_dirpath = build_bin_restart_dirpath(RestartInPath)//"/np"//int_to_txt(partit%npes)
  read_bin_infopath = build_bin_restart_infopath(RestartInPath)//"/np"//int_to_txt(partit%npes)//".info"
  read_oce_path = nc_restart_path('oce', yearold, RestartInPath)
  read_ice_path = nc_restart_path('ice', yearold, RestartInPath)
  read_bio_path = nc_restart_path('bio', yearold, RestartInPath)
  
  ! Initialize file groups for reading
  call ini_ocean_io(yearold, dynamics, tracers, partit, mesh)
  if (use_ice) then
#if defined(__icepack)    
      call ini_icepack_io(yearold, partit, mesh)
#else
      call ini_ice_io  (yearold, ice, partit, mesh)
#endif        
  end if     
#if defined(__recom)
  if (REcoM_restart) call ini_bio_io(yearold, tracers, partit, mesh)
#endif
  
  ! Check for raw restart files
  if(partit%mype == RAW_RESTART_METADATA_RANK) then
    inquire(file=read_raw_infopath, exist=rawfiles_exist)
  end if
  call MPI_Bcast(rawfiles_exist, 1, MPI_LOGICAL, RAW_RESTART_METADATA_RANK, partit%MPI_COMM_FESOM, mpierr)
  
  ! Check for binary restart files
  if(partit%mype == RAW_RESTART_METADATA_RANK) then
    inquire(file=read_bin_infopath, exist=binfiles_exist)
  end if
  call MPI_Bcast(binfiles_exist, 1, MPI_LOGICAL, RAW_RESTART_METADATA_RANK, partit%MPI_COMM_FESOM, mpierr)
  
  ! Read restart files in order of preference
  if(rawfiles_exist) then
    ! Read raw/core dump restart
    which_readr = 1
    ! Note: This will need to be updated once we have read functions that accept paths
    call read_all_raw_restarts(partit%MPI_COMM_FESOM, partit%mype)
    
  elseif(binfiles_exist .and. bin_restart_length_unit /= "off") then
    ! Read binary restart
    which_readr = 2
    if (use_ice) then 
        call read_all_bin_restarts(read_bin_dirpath, &
                                   partit   = partit,   &
                                   mesh     = mesh,     &
                                   ice      = ice,      &
                                   dynamics = dynamics, &
                                   tracers  = tracers   )
    else
        call read_all_bin_restarts(read_bin_dirpath, &
                                   partit   = partit,   &
                                   mesh     = mesh,     &                    
                                   dynamics = dynamics, &
                                   tracers  = tracers   )
    end if     
    
  else
    ! Read NetCDF restart
    which_readr = 0
    
    ! Read OCEAN restart
    if (partit%mype==RAW_RESTART_METADATA_RANK) print *, achar(27)//'[1;33m'//' --> read restarts from netcdf file: ocean'//achar(27)//'[0m'
    call read_restart(read_oce_path, oce_files, partit%MPI_COMM_FESOM, partit%mype)
    
    ! Read ICE/ICEPACK restart
    if (use_ice) then
#if defined(__icepack)   
        if (partit%mype==RAW_RESTART_METADATA_RANK) print *, achar(27)//'[1;33m'//' --> read restarts from netcdf file: icepack'//achar(27)//'[0m'
        call read_restart(nc_restart_path('icepack', yearold, RestartInPath), icepack_files, partit%MPI_COMM_FESOM, partit%mype)
#else            
        if (partit%mype==RAW_RESTART_METADATA_RANK) print *, achar(27)//'[1;33m'//' --> read restarts from netcdf file: ice'//achar(27)//'[0m'
        call read_restart(read_ice_path, ice_files, partit%MPI_COMM_FESOM, partit%mype)            
#endif
    end if 

#if defined(__recom)
    ! Read RECOM restarts
    if (REcoM_restart) then
        if (partit%mype==RAW_RESTART_METADATA_RANK) print *, achar(27)//'[1;33m'//' --> read restarts from netcdf file: bio'//achar(27)//'[0m'
        call read_restart(read_bio_path, bio_files, partit%MPI_COMM_FESOM, partit%mype)
    end if
#endif

    ! Immediately create raw and binary restarts after NetCDF read for backup
    if(raw_restart_length_unit /= "off") then
        call write_all_raw_restarts(0, partit%MPI_COMM_FESOM, partit%mype)
    end if
    
    if(bin_restart_length_unit /= "off") then
        call write_all_bin_restarts((/globalstep, int(ctime), yearnew/), &
                                    bin_restart_dirpath,                 &
                                    bin_restart_infopath,                &
                                    partit,                              &
                                    mesh,                                &                                        
                                    ice,                                 &
                                    dynamics,                            &
                                    tracers                              )
    end if
  end if
  
end subroutine read_initial_conditions

!--------------------------------------------------------------------------------------------
! Separate subroutine for writing restart files
subroutine write_initial_conditions(istep, nstart, ntotal, which_readr, ice, dynamics, tracers, partit, mesh)
  use fortran_utils
  implicit none
  
  ! Parameters
  integer :: istep, nstart, ntotal
  type(t_mesh)  , intent(inout), target :: mesh
  type(t_partit), intent(inout), target :: partit
  type(t_tracer), intent(inout), target :: tracers
  type(t_dyn)   , intent(inout), target :: dynamics
  type(t_ice)   , intent(inout), target :: ice
  integer, intent(in) :: which_readr
  
  ! Local variables
  logical :: is_portable_restart_write, is_raw_restart_write, is_bin_restart_write
  
  ! Skip writing on step 0
  if (istep==0) return
  
  ! Check whether restart will be written
  is_portable_restart_write = is_due(trim(restart_length_unit), restart_length, istep)
  
  ! Should write core dump restart?
  if(is_portable_restart_write .and. (raw_restart_length_unit /= "off")) then
    is_raw_restart_write = .true. ! always write a raw restart together with the portable restart
  else
#if !defined __ifsinterface
    is_raw_restart_write = is_due(trim(raw_restart_length_unit), raw_restart_length, istep)
#else
    is_raw_restart_write = is_due(trim(raw_restart_length_unit), raw_restart_length, istep) .OR. (istep==ntotal)
#endif
  end if
  
  ! Should write derived type binary restart?
  if(is_portable_restart_write .and. (bin_restart_length_unit /= "off")) then
    is_bin_restart_write = .true. ! always write a binary restart together with the portable restart
  else
    is_bin_restart_write = is_due(trim(bin_restart_length_unit), bin_restart_length, istep)
  end if

  ! Write restart files
  if(is_portable_restart_write) then
    ! Write OCEAN restart
    if (partit%mype==RAW_RESTART_METADATA_RANK) print *, achar(27)//'[1;33m'//' --> write restarts to netcdf file: ocean'//achar(27)//'[0m'
    call write_restart(oce_path, oce_files, istep)
    
    ! Write ICE/ICEPACK restart
    if(use_ice) then
#if defined(__icepack)        
        if (partit%mype==RAW_RESTART_METADATA_RANK) print *, achar(27)//'[1;33m'//' --> write restarts to netcdf file: icepack'//achar(27)//'[0m'
        call write_restart(icepack_path, icepack_files, istep)
#else
        if (partit%mype==RAW_RESTART_METADATA_RANK) print *, achar(27)//'[1;33m'//' --> write restarts to netcdf file: ice'//achar(27)//'[0m'
        call write_restart(ice_path, ice_files, istep)
#endif 
    end if

#if defined(__recom)
    ! Write RECOM restart
    if (REcoM_restart .or. use_REcoM) then
        if (partit%mype==RAW_RESTART_METADATA_RANK) print *, achar(27)//'[1;33m'//' --> write restarts to netcdf file: bio'//achar(27)//'[0m'
        call write_restart(bio_path, bio_files, istep)
    end if
#endif
  end if

  ! Write core dump
  if(is_raw_restart_write) then
    call write_all_raw_restarts(istep, partit%MPI_COMM_FESOM, partit%mype)
  end if

  ! Write derived type binary
  if(is_bin_restart_write) then
    call write_all_bin_restarts((/globalstep+istep, int(ctime), yearnew/), &
                                bin_restart_dirpath,                 &
                                bin_restart_infopath,                &
                                partit,                              &
                                mesh,                                &                                
                                ice,                                 &
                                dynamics,                            &
                                tracers                              )
  end if

  ! Update clock file to latest restart point
  if (partit%mype==0) then
    if(is_portable_restart_write .or. is_raw_restart_write .or. is_bin_restart_write) then
        write(*,*) ' --> actualize clock file to latest restart point'
        call clock_finish
    end if
  end if

end subroutine write_initial_conditions

!
!--------------------------------------------------------------------------------------------
! Original restart subroutine now acting as a wrapper for backward compatibility
subroutine restart(istep, nstart, ntotal, l_read, which_readr, ice, dynamics, tracers, partit, mesh)

  use fortran_utils

  implicit none
  ! Wrapper subroutine for backward compatibility
  ! Calls appropriate read/write subroutines based on l_read flag

  integer :: istep, nstart, ntotal
  logical :: l_read
  type(t_mesh)  , intent(inout), target :: mesh
  type(t_partit), intent(inout), target :: partit
  type(t_tracer), intent(inout), target :: tracers
  type(t_dyn)   , intent(inout), target :: dynamics
  type(t_ice)   , intent(inout), target :: ice
  integer, intent(out):: which_readr
  
  ! Local variables for initialization
  logical, save :: initialized_raw = .false.
  logical, save :: initialized_bin = .false.
  integer :: mpierr
  
  !_____________________________________________________________________________
  ! initialize directory for core dump restart 
  if(.not. initialized_raw) then
    initialized_raw = .true.
    raw_restart_dirpath  = build_raw_restart_dirpath(RestartOutPath)//"/np"//int_to_txt(partit%npes)
    raw_restart_infopath = build_raw_restart_infopath(RestartOutPath)//"/np"//int_to_txt(partit%npes)//".info"
    if(raw_restart_length_unit /= "off") then
      if(partit%mype == RAW_RESTART_METADATA_RANK) then
        ! execute_command_line with mkdir sometimes fails, use a custom implementation around mkdir from C instead
        call mkdir(build_raw_restart_dirpath(RestartOutPath)) ! we have no mkdir -p, create the intermediate dirs separately
        call mkdir(raw_restart_dirpath)
      end if
      call MPI_Barrier(partit%MPI_COMM_FESOM, mpierr) ! make sure the dir has been created before we continue...
    end if
  end if

  !_____________________________________________________________________________
  ! initialize directory for derived type binary restart
  if(.not. initialized_bin) then
    initialized_bin = .true.
    bin_restart_dirpath  = build_bin_restart_dirpath(RestartOutPath)//"/np"//int_to_txt(partit%npes)
    bin_restart_infopath = build_bin_restart_infopath(RestartOutPath)//"/np"//int_to_txt(partit%npes)//".info"
    if(bin_restart_length_unit /= "off") then
        if(partit%mype == RAW_RESTART_METADATA_RANK) then
            ! execute_command_line with mkdir sometimes fails, use a custom implementation around mkdir from C instead
            call mkdir(build_bin_restart_dirpath(RestartOutPath)) ! we have no mkdir -p, create the intermediate dirs separately
            call mkdir(bin_restart_dirpath)
        end if
        call MPI_Barrier(partit%MPI_COMM_FESOM, mpierr) ! make sure the dir has been created before we continue...
    end if
  end if
  
  !_____________________________________________________________________________
  ! Initialize file groups for writing (when not reading)
  if (.not. l_read) then
    call ini_ocean_io(yearnew, dynamics, tracers, partit, mesh)
    if (use_ice) then
#if defined(__icepack)
        call ini_icepack_io(yearnew, partit, mesh)
#else        
        call ini_ice_io(yearnew, ice, partit, mesh)        
#endif        
    end if     
#if defined(__recom)
    if (use_REcoM) call ini_bio_io  (yearnew, tracers, partit, mesh)
#endif

  else
    call ini_ocean_io(yearold, dynamics, tracers, partit, mesh)
    if (use_ice) then
#if defined(__icepack)    
        call ini_icepack_io(yearold, partit, mesh)
#else
        call ini_ice_io  (yearold, ice, partit, mesh)
#endif        
    end if     
#if defined(__recom)
    if (REcoM_restart) call ini_bio_io(yearold, tracers, partit, mesh)
#endif
  end if ! --> if (.not. l_read) then

  !_____________________________________________________________________________
  ! Call appropriate subroutines based on operation type
  if (l_read) then
    ! Read initial conditions
    call read_initial_conditions(which_readr, ice, dynamics, tracers, partit, mesh)
  else
    ! Write restart files
    call write_initial_conditions(istep, nstart, ntotal, which_readr, ice, dynamics, tracers, partit, mesh)
  end if


end subroutine restart
!
!
!_______________________________________________________________________________
subroutine write_restart(path, filegroup, istep)
  use fortran_utils
  character(len=*), intent(in) :: path
  type(restart_file_group), intent(inout) :: filegroup
  integer,  intent(in)          :: istep
  ! EO parameters
  integer cstep
  integer i
  character(:), allocatable :: dirpath
  character(:), allocatable :: filepath
  logical file_exists
  
  cstep = globalstep+istep
  
  do i=1, filegroup%nfiles
    call filegroup%files(i)%join() ! join the previous write (if required)

    if(filegroup%files(i)%is_iorank()) then
      if(filegroup%files(i)%is_attached()) call filegroup%files(i)%close_file() ! close the file from previous write
            
      dirpath = path(1:len(path)-3) ! chop of the ".nc" suffix
      filepath = dirpath//"/"//filegroup%files(i)%varname//".nc"
      if(filegroup%files(i)%path == "" .or. (.not. filegroup%files(i)%must_exist_on_read)) then
        ! the path to an existing restart file is not set in read_restart if we had a restart from a raw restart
        ! OR we might have skipped the file when reading restarts and it does not exist at all
        inquire(file=filepath, exist=file_exists)
        if(file_exists) then
          filegroup%files(i)%path = filepath
        else if(.not. filegroup%files(i)%must_exist_on_read) then
          filegroup%files(i)%path = ""
        end if
      end if
      if(filegroup%files(i)%path .ne. filepath) then
        ! execute_command_line with mkdir sometimes fails, use a custom implementation around mkdir from C instead
        call mkdir(dirpath)
        filegroup%files(i)%path = filepath
        call filegroup%files(i)%open_write_create(filegroup%files(i)%path)
      else
        call filegroup%files(i)%open_write_append(filegroup%files(i)%path) ! todo: keep the file open between writes
      end if

      write(*,*) 'writing restart record ', filegroup%files(i)%rec_count()+1, ' to ', filegroup%files(i)%path
      call filegroup%files(i)%write_var(filegroup%files(i)%iter_varindex, [filegroup%files(i)%rec_count()+1], [1], [cstep])
      ! todo: write time via the fesom_file_type
      call filegroup%files(i)%write_var(filegroup%files(i)%time_varindex(), [filegroup%files(i)%rec_count()+1], [1], [ctime])
    end if

    call filegroup%files(i)%async_gather_and_write_variables()
  end do
  
end subroutine write_restart
!
!
!_______________________________________________________________________________
subroutine write_all_raw_restarts(istep, mpicomm, mype)
  integer,  intent(in):: istep
  integer, intent(in) :: mpicomm
  integer, intent(in) :: mype
  ! EO parameters
  integer cstep
  integer fileunit

  open(newunit = fileunit, file = raw_restart_dirpath//'/'//mpirank_to_txt(mpicomm)//'.dump', form = 'unformatted')
  call write_raw_restart_group(oce_files, fileunit)
  if(use_ice) call write_raw_restart_group(ice_files, fileunit)
#if defined(__recom)
  call write_raw_restart_group(bio_files, fileunit)
#endif
  close(fileunit)

  if(mype == RAW_RESTART_METADATA_RANK) then
    print *,"writing raw restart to "//raw_restart_dirpath
    ! store metadata about the raw restart
    cstep = globalstep+istep
    open(newunit = fileunit, file = raw_restart_infopath)
    write(fileunit, '(g0)') cstep
    write(fileunit, '(g0)') ctime
    write(fileunit, '(2(g0))') "! year: ",yearnew
    write(fileunit, '(3(g0))') "! oce: ", oce_files%nfiles, " variables"
    if(use_ice) write(fileunit, '(3(g0))') "! ice: ", ice_files%nfiles, " variables"
#if defined(__recom)
    write(fileunit, '(3(g0))') "! bio: ", bio_files%nfiles, " variables"
#endif
    close(fileunit)
  end if
end subroutine write_all_raw_restarts
!
!
!_______________________________________________________________________________
subroutine write_raw_restart_group(filegroup, fileunit)
  type(restart_file_group), intent(inout) :: filegroup
  integer, intent(in) :: fileunit
  ! EO parameters
  integer i
  
  do i=1, filegroup%nfiles
    call filegroup%files(i)%write_variables_raw(fileunit)
  end do
end subroutine write_raw_restart_group
! ! !
! ! !
! ! !_______________________________________________________________________________
! ! subroutine write_all_bin_restarts(istep, ice, dynamics, tracers, partit, mesh)
! !     integer, intent(in) :: istep
! !     type(t_ice)   , target, intent(in) :: ice
! !     type(t_dyn)   , target, intent(in) :: dynamics
! !     type(t_tracer), target, intent(in) :: tracers
! !     type(t_partit), target, intent(in) :: partit
! !     type(t_mesh)  , target, intent(in) :: mesh
! !     
! !     ! EO parameters
! !     integer cstep
! !     integer fileunit, fileunit_i
! !     
! !     !___________________________________________________________________________
! !     ! write info file
! !     if(partit%mype == RAW_RESTART_METADATA_RANK) then
! !         print *, achar(27)//'[1;33m'//' --> writing derived type binary restarts to '//bin_restart_dirpath//achar(27)//'[0m'
! !         ! store metadata about the raw restart
! !         cstep = globalstep+istep
! !         fileunit_i = 299
! !         open(newunit = fileunit_i, file = bin_restart_infopath)
! !         write(fileunit_i, '(g0)') cstep
! !         write(fileunit_i, '(g0)') ctime
! !         write(fileunit_i, '(2(g0))') "! year: ",yearnew
! !     end if 
! !     
! !     !___________________________________________________________________________
! !     ! mesh derived type 
! !     fileunit = partit%mype+300
! !     open(newunit = fileunit, &
! !         file     = bin_restart_dirpath//'/'//'t_mesh.'//mpirank_to_txt(partit%MPI_COMM_FESOM), &
! !         status   = 'replace', &
! !         form     = 'unformatted')
! !     write(fileunit) mesh
! !     close(fileunit)
! !     if(partit%mype == RAW_RESTART_METADATA_RANK) then
! !         write(fileunit_i, '(1(g0))') "!   t_mesh"
! !         print *, achar(27)//'[33m'//'     > write derived type t_mesh'//achar(27)//'[0m'
! !     end if     
! !     
! !     !___________________________________________________________________________
! !     ! partit derived type 
! !     fileunit = partit%mype+300
! !     open(newunit = fileunit, &
! !         file     = bin_restart_dirpath//'/'//'t_partit.'//mpirank_to_txt(partit%MPI_COMM_FESOM), &
! !         status   = 'replace', &
! !         form     = 'unformatted')
! !     write(fileunit) partit
! !     close(fileunit)
! !     if(partit%mype == RAW_RESTART_METADATA_RANK) then 
! !         write(fileunit_i, '(1(g0))') "!   t_partit"
! !         print *, achar(27)//'[33m'//'     > write derived type t_partit'//achar(27)//'[0m'
! !     end if 
! !     
! !     !___________________________________________________________________________
! !     ! tracer derived type 
! !     fileunit = partit%mype+300
! !     open(newunit = fileunit, &
! !         file     = bin_restart_dirpath//'/'//'t_tracer.'//mpirank_to_txt(partit%MPI_COMM_FESOM), &
! !         status   = 'replace', &
! !         form     = 'unformatted')
! !     write(fileunit) tracers  
! !     close(fileunit)
! !     if(partit%mype == RAW_RESTART_METADATA_RANK) then 
! !         write(fileunit_i, '(1(g0))') "!   t_tracer"
! !         print *, achar(27)//'[33m'//'     > write derived type t_tracer'//achar(27)//'[0m'
! !     end if     
! !     
! !     !___________________________________________________________________________
! !     ! dynamics derived type 
! !     fileunit = partit%mype+300
! !     open(newunit = fileunit, &
! !         file     = bin_restart_dirpath//'/'//'t_dynamics.'//mpirank_to_txt(partit%MPI_COMM_FESOM), &
! !         status   = 'replace', &
! !         form     = 'unformatted')
! !     write(fileunit) dynamics
! !     close(fileunit)
! !     if(partit%mype == RAW_RESTART_METADATA_RANK) then 
! !         write(fileunit_i, '(1(g0))') "!   t_dynamics"
! !         print *, achar(27)//'[33m'//'     > write derived type t_dynamics'//achar(27)//'[0m'
! !     end if     
! !     
! !     !___________________________________________________________________________
! !     ! ice derived type 
! !     if (use_ice) then
! !         fileunit = partit%mype+300
! !         open(newunit = fileunit, &
! !             file    = bin_restart_dirpath//'/'//'t_ice.'//mpirank_to_txt(partit%MPI_COMM_FESOM), &
! !             status  = 'replace', &
! !             form    = 'unformatted')
! !         write(fileunit) ice
! !         close(fileunit)
! !         if(partit%mype == RAW_RESTART_METADATA_RANK) then 
! !             write(fileunit_i, '(1(g0))') "!   t_ice"
! !             print *, achar(27)//'[33m'//'     > write derived type t_ice'//achar(27)//'[0m'
! !         end if     
! !     end if 
! !     
! !     !___________________________________________________________________________
! !     if(partit%mype == RAW_RESTART_METADATA_RANK) close(fileunit_i)
! ! 
! ! end subroutine
! ! !
! ! !
! ! !_______________________________________________________________________________
! ! subroutine read_all_bin_restarts(path_in, ice, dynamics, tracers, partit, mesh)
! !     implicit none 
! !     
! !     ! do optional here for the usage with dwarfs, since there only specific derived  
! !     ! types will be needed
! !     character(len=*), intent(in)                    :: path_in
! !     type(t_ice)   , intent(inout), target, optional :: ice
! !     type(t_dyn)   , intent(inout), target, optional :: dynamics
! !     type(t_tracer), intent(inout), target, optional :: tracers
! !     type(t_partit), intent(inout), target, optional :: partit
! !     type(t_mesh)  , intent(inout), target, optional :: mesh
! !     integer fileunit
! !     
! !     if (partit%mype==RAW_RESTART_METADATA_RANK) print *, achar(27)//'[1;33m'//' --> read restarts from derived type binary'//achar(27)//'[0m'
! !     
! !     !___________________________________________________________________________
! !     ! mesh derived type 
! !     if (present(mesh)) then
! !         fileunit = partit%mype+300
! !         open(newunit = fileunit, &
! !             file     = trim(path_in)//'/'//'t_mesh.'//mpirank_to_txt(partit%MPI_COMM_FESOM), &
! !             status   = 'old', &
! !             form     = 'unformatted')
! !         read(fileunit) mesh
! !         close(fileunit)
! !         if (partit%mype==RAW_RESTART_METADATA_RANK) print *, achar(27)//'[33m'//'     > read derived type t_mesh'//achar(27)//'[0m'
! !     end if
! !     
! !     !___________________________________________________________________________
! !     ! partit derived type 
! !     if (present(partit)) then
! !         fileunit = partit%mype+300
! !         open(newunit = fileunit, &
! !             file     = trim(path_in)//'/'//'t_partit.'//mpirank_to_txt(partit%MPI_COMM_FESOM), &
! !             status   = 'old', &
! !             form     = 'unformatted')
! !         read(fileunit) partit
! !         close(fileunit)
! !         if (partit%mype==RAW_RESTART_METADATA_RANK) print *, achar(27)//'[33m'//'     > read derived type t_partit'//achar(27)//'[0m'
! !     end if 
! !     
! !     !___________________________________________________________________________
! !     ! tracer derived type     
! !     if (present(tracers)) then
! !         fileunit = partit%mype+300
! !         open(newunit = fileunit, &
! !             file     = trim(path_in)//'/'//'t_tracer.'//mpirank_to_txt(partit%MPI_COMM_FESOM), &
! !             status   = 'old', &
! !             form     = 'unformatted')
! !         read(fileunit) tracers  
! !         close(fileunit)
! !         if (partit%mype==RAW_RESTART_METADATA_RANK) print *, achar(27)//'[33m'//'     > read derived type t_tracer'//achar(27)//'[0m'
! !     end if 
! !     
! !     !___________________________________________________________________________
! !     ! dynamics derived type 
! !     if (present(dynamics)) then
! !         fileunit = partit%mype+300
! !         open(newunit = fileunit, &
! !             file     = trim(path_in)//'/'//'t_dynamics.'//mpirank_to_txt(partit%MPI_COMM_FESOM), &
! !             status   = 'old', &
! !             form     = 'unformatted')
! !         read(fileunit) dynamics
! !         close(fileunit)
! !         if (partit%mype==RAW_RESTART_METADATA_RANK) print *, achar(27)//'[33m'//'     > read derived type t_dynamics'//achar(27)//'[0m'
! !     end if 
! !     
! !     !___________________________________________________________________________
! !     ! ice derived type 
! !     if (present(ice)) then
! !         fileunit = partit%mype+300
! !         open(newunit = fileunit, &
! !             file    = trim(path_in)//'/'//'t_ice.'//mpirank_to_txt(partit%MPI_COMM_FESOM), &
! !             status  = 'old', &
! !             form    = 'unformatted')
! !         read(fileunit) ice
! !         close(fileunit)
! !         if (partit%mype==RAW_RESTART_METADATA_RANK) print *, achar(27)//'[33m'//'     > read derived type t_ice'//achar(27)//'[0m'
! !     end if 
! ! end subroutine
!
!
!_______________________________________________________________________________
subroutine read_all_raw_restarts(mpicomm, mype)
  integer, intent(in) :: mpicomm
  integer, intent(in) :: mype
  ! EO parameters
  integer rstep
  real(kind=WP) rtime
  integer fileunit
  integer status
  integer mpierr

  if(mype == RAW_RESTART_METADATA_RANK) then
    ! read metadata info for the raw restart
    open(newunit = fileunit, status = 'old', iostat = status, file = raw_restart_infopath)
    if(status == 0) then
      read(fileunit,*) rstep
      read(fileunit,*) rtime
      close(fileunit)
    else
      print *,"can not open ",raw_restart_infopath
      stop 1
    end if
    
    ! compare the restart time with our actual time
    if(int(ctime) /= int(rtime)) then
      print *, "raw restart time ",rtime,"does not match current clock time",ctime
      stop 1
    end if
    globalstep = rstep
    print *,"reading raw restart from "//raw_restart_dirpath
  end if
  ! sync globalstep with the other processes to let all processes writing portable restart files know the globalstep
  call MPI_Bcast(globalstep, 1, MPI_INTEGER, RAW_RESTART_METADATA_RANK, mpicomm, mpierr)

  open(newunit = fileunit, status = 'old', iostat = status, file = raw_restart_dirpath//'/'//mpirank_to_txt(mpicomm)//'.dump', form = 'unformatted')
  if(status == 0) then
    call read_raw_restart_group(oce_files, fileunit)
    if(use_ice) call read_raw_restart_group(ice_files, fileunit)
#if defined(__recom)
    call read_raw_restart_group(bio_files, fileunit)
#endif
    close(fileunit)
  else
    print *,"can not open ",raw_restart_dirpath//'/'//mpirank_to_txt(mpicomm)//'.dump'
    stop 1
  end if
end subroutine read_all_raw_restarts
!
!
!_______________________________________________________________________________
subroutine read_raw_restart_group(filegroup, fileunit)
  type(restart_file_group), intent(inout) :: filegroup
  integer, intent(in) :: fileunit
  ! EO parameters
  integer i
  
  do i=1, filegroup%nfiles
    call filegroup%files(i)%read_variables_raw(fileunit)
  end do  
end subroutine read_raw_restart_group
!
!
!_______________________________________________________________________________
! join remaining threads and close all open files
subroutine finalize_restart()
  integer i

  ! join all previous writes
  ! close all restart files

  do i=1, oce_files%nfiles
    call oce_files%files(i)%join()
    if(oce_files%files(i)%is_iorank()) then
      if(oce_files%files(i)%is_attached()) call oce_files%files(i)%close_file()
    end if
  end do

  if(use_ice) then
    do i=1, ice_files%nfiles
      call ice_files%files(i)%join()
      if(ice_files%files(i)%is_iorank()) then
        if(ice_files%files(i)%is_attached()) call ice_files%files(i)%close_file()
      end if
    end do
  end if
#if defined(__recom)
  do i=1, bio_files%nfiles
    call bio_files%files(i)%join()
    if(bio_files%files(i)%is_iorank()) then
      if(bio_files%files(i)%is_attached()) call bio_files%files(i)%close_file()
    end if
  end do
#endif
end subroutine finalize_restart
!
!
!_______________________________________________________________________________
subroutine read_restart(path, filegroup, mpicomm, mype)
  character(len=*), intent(in) :: path
  type(restart_file_group), intent(inout) :: filegroup
  integer, intent(in) :: mpicomm
  integer, intent(in) :: mype
  ! EO parameters
  real(kind=WP) rtime
  integer i
  character(:), allocatable :: dirpath
  integer mpistatus(MPI_STATUS_SIZE)
  logical file_exists
  logical, allocatable :: skip_file(:)
  integer current_iorank_snd, current_iorank_rcv
  integer max_globalstep
  integer mpierr
  
  allocate(skip_file(filegroup%nfiles))
  skip_file = .false.
  
  do i=1, filegroup%nfiles
    current_iorank_snd = 0
    current_iorank_rcv = 0
    if( filegroup%files(i)%is_iorank() ) then
      dirpath = path(1:len(path)-3) ! chop of the ".nc" suffix
      if(filegroup%files(i)%path .ne. dirpath//"/"//filegroup%files(i)%varname//".nc") then
        filegroup%files(i)%path = dirpath//"/"//filegroup%files(i)%varname//".nc"

        ! determine if the file should be skipped
        if(.not. filegroup%files(i)%must_exist_on_read) then
          current_iorank_snd = mype
          inquire(file=filegroup%files(i)%path, exist=file_exists)
          if(.not. file_exists) skip_file(i) = .true.
        end if

        if(.not. skip_file(i)) then
#ifndef DISABLE_PARALLEL_RESTART_READ
          write(*,*) 'reading restart PARALLEL for ', filegroup%files(i)%varname, ' at ', filegroup%files(i)%path
#else
          write(*,*) 'reading restart SEQUENTIAL for ', filegroup%files(i)%varname, ' at ', filegroup%files(i)%path
#endif
        else
#ifndef DISABLE_PARALLEL_RESTART_READ
          write(*,*) 'skipping reading restart PARALLEL for ', filegroup%files(i)%varname, ' at ', filegroup%files(i)%path
#else
          write(*,*) 'skipping reading restart SEQUENTIAL for ', filegroup%files(i)%varname, ' at ', filegroup%files(i)%path
#endif
        end if
        
        if(.not. skip_file(i)) call filegroup%files(i)%open_read(filegroup%files(i)%path) ! do we need to bother with read-only access?
        ! todo: print a reasonable error message if the file does not exist
      end if      
    end if

    ! iorank already knows if we skip the file, tell the others
    if(.not. filegroup%files(i)%must_exist_on_read) then
      call MPI_Allreduce(current_iorank_snd, current_iorank_rcv, 1, MPI_INTEGER, MPI_SUM, mpicomm, mpierr)
      call MPI_Bcast(skip_file(i), 1, MPI_LOGICAL, current_iorank_rcv, mpicomm, mpierr)
    end if      

    if(.not. skip_file(i)) call filegroup%files(i)%async_read_and_scatter_variables()
#ifndef DISABLE_PARALLEL_RESTART_READ
  end do
  
  do i=1, filegroup%nfiles
#endif
    if(skip_file(i)) cycle
    call filegroup%files(i)%join()

    if(filegroup%files(i)%is_iorank()) then
      write(*,*) 'restart from record ', filegroup%files(i)%rec_count(), ' of ', filegroup%files(i)%rec_count(), filegroup%files(i)%path

      ! read the last entry from the iter variable
      call filegroup%files(i)%read_var1(filegroup%files(i)%iter_varindex, [filegroup%files(i)%rec_count()], globalstep)

      ! read the last entry from the time variable
      call filegroup%files(i)%read_var1(filegroup%files(i)%time_varindex(), [filegroup%files(i)%rec_count()], rtime)
      call filegroup%files(i)%close_file()

     if (int(ctime)/=int(rtime)) then
        write(*,*) 'Reading restart: timestamps in restart and in clock files do not match for ', filegroup%files(i)%varname, ' at ', filegroup%files(i)%path
        write(*,*) 'restart/ times are:', ctime, rtime
        write(*,*) 'the model will stop!'
        call par_ex(mpicomm, mype, 1)
      end if
    end if
  end do

  ! sync globalstep with processes which may have skipped a restart upon reading and thus need to know the globalstep when writing their restart 
  if( any(skip_file .eqv. .true.) ) then
    call MPI_Allreduce(globalstep, max_globalstep, 1, MPI_INTEGER, MPI_MAX, mpicomm, mpierr)
    globalstep = max_globalstep
  end if

  ! sync globalstep with the process responsible for raw restart metadata
  if(filegroup%nfiles >= 1) then
    ! use the first restart I/O process to send the globalstep
    if( filegroup%files(1)%is_iorank() .and. (mype .ne. RAW_RESTART_METADATA_RANK)) then
      call MPI_Send(globalstep, 1, MPI_INTEGER, RAW_RESTART_METADATA_RANK, 42, mpicomm, mpierr)
    else if((mype == RAW_RESTART_METADATA_RANK) .and. (.not. filegroup%files(1)%is_iorank())) then
      call MPI_Recv(globalstep, 1, MPI_INTEGER, MPI_ANY_SOURCE, 42, mpicomm, mpistatus, mpierr)
    end if
  end if
end subroutine read_restart
!
!
!_______________________________________________________________________________
  function is_due(unit, frequency, istep) result(d)
    character(len=*), intent(in) :: unit
    integer, intent(in) :: frequency
    integer, intent(in) :: istep
    logical d
    ! EO parameters
    d = .false.
    
    if(unit.eq.'y') then
      call annual_event(d)
    else if(unit.eq.'m') then 
      call monthly_event(d) 
    else if(unit.eq.'d') then
      call daily_event(d, frequency)
    else if(unit.eq.'h') then
      call hourly_event(d, frequency)
    else if(unit.eq.'s') then
      call step_event(d, istep, frequency)
    else if(unit.eq.'off') then
      d = .false.
    else
      write(*,*) 'You did not specify a supported outputflag.'
      write(*,*) 'The program will stop to give you opportunity to do it.'
      stop 1
    stop
    end if
  end function
! !
! !
! !_______________________________________________________________________________
!   function mpirank_to_txt(mpicomm) result(txt)
!     use fortran_utils
!     integer, intent(in) :: mpicomm
!     character(:), allocatable :: txt
!     ! EO parameters
!     integer mype
!     integer npes
!     integer mpierr
!   
!     call MPI_Comm_Rank(mpicomm, mype, mpierr)
!     call MPI_Comm_Size(mpicomm, npes, mpierr)
!     txt = int_to_txt_pad(mype,int(log10(real(npes)))+1) ! pad to the width of the number of processes
!   end function
!!PS --> move this function also to fortran_utils.F90 

end module
