!
! Module of model configuration parameters + routines to set model configuration
! It combines gen_modules_config and gen_setup_model of FESOM release 
! in a single module. 
! S. Danilov, 5.04.2012
!
! ======================================================================
module g_config
  use o_param
  implicit none
  save

  ! *** Modelname ***
  character(5)           :: runid='test1'       ! a model/setup name
  namelist /modelname/ runid

  ! *** time step ***
  integer                :: step_per_day=72     ! number of steps per day
  integer                :: run_length=1	! run length
  character              :: run_length_unit='y' ! unit: y, d, s
  namelist /timestep/ step_per_day, run_length, run_length_unit
  
  ! *** Paths for all in and out ***
  character(100)         :: MeshPath='./mesh/'
  character(100)         :: OpbndPath='./opbnd/'
  character(100)         :: ClimateDataPath='./hydrography/'
  character(100)         :: ForcingDataPath='./forcing/'
  character(100)         :: TideForcingPath='./tide_forcing/'
  character(100)         :: ResultPath='./result/'
  namelist /paths/  MeshPath, OpbndPath, ClimateDataPath, ForcingDataPath, &
       TideForcingPath, ResultPath

  ! *** ocean climatology data name ***
  character(100)         :: OceClimaDataName='annual_woa01_ts.out'
  logical                :: use_prepared_init_ice=.false.  !initialize ice externally
  namelist /initialization/ OceClimaDataName, use_prepared_init_ice

  ! *** in out ***
  character*4            :: restartflag='last'     !restart from record,'#','last'
  integer                :: output_length=1        !valid for d,h,s
  character              :: output_length_unit='m' !output period: y, m, d, h, s 
  integer                :: logfile_outfreq=1      ! logfile info. outp. freq., # steps
  logical                :: use_means=.true.       !Mean or snapshot
  integer                :: restart_length=1
  character              :: restart_length_unit='m'
  integer                :: output_offset=32, restart_offset=32

  namelist /inout/ restartflag, output_length, output_length_unit, restart_length, restart_length_unit, &
	logfile_outfreq, use_means, output_offset, restart_offset

  ! *** mesh ***
  integer                :: grid_type=1 	! z-level, 2 sigma, 3 sigma + z-level
  logical                :: use_ALE=.false.     ! switch on/off ALE
  character(20)          :: which_ALE='zlevel' ! 'zlevel', 'zstar', 'zstar-weighted', 'ztilde'
  logical                :: use_partial_cell=.false.  ! use partial bottom cell configuration  
! for zlevel: layer thickness should not become smaller than min_hnode of 
  ! original layer thickness. If it happens switch from zelvel to local zstar
  real(kind=WP)          :: min_hnode=0.25 
  ! for zlevel: in case min_hnode criteria is reached over how many level should 
  ! ssh change be distributed
  integer                :: lzstar_lev=3
  ! maximal pressure from ice felt by the ocean
  real(kind=8)           :: max_ice_loading=5.0

  namelist /mesh_def/ grid_type, use_ALE, which_ALE, use_partial_cell, min_hnode, lzstar_lev, max_ice_loading

  ! *** model geometry
  logical                :: cartesian=.false.
  logical                :: fplane=.false.
  logical                :: betaplane=.false.      ! not used
  real(kind=WP)          :: f_fplane=-1.4e-4       ! [1/s] not used
  real(kind=WP)          :: beta_betaplane=2.0e-11 ! [1/s/m] not used
  real(kind=WP)          :: cyclic_length=360.     ! [degree]
                                         ! if there is no cyclicity just use 
	                                 ! value which is larger than maximum 
			                 ! triangle size.
  logical                :: rotated_grid=.true. ! not used
  logical                :: force_rotation=.true.
  real(kind=WP)          :: alphaEuler=50. 	! [degree] Euler angles, convention:
  real(kind=WP)          :: betaEuler=15.  	! first around z, then around new x,
  real(kind=WP)		 :: gammaEuler=-90.	! then around new z.
  				                ! Set to zeros to work with
						! geographical coordinates
  namelist /geometry/  cartesian, fplane, betaplane, f_fplane, beta_betaplane, &
       cyclic_length, rotated_grid, alphaEuler, betaEuler, gammaEuler, force_rotation

  ! *** fleap_year ***
  logical                       :: include_fleapyear=.false.
  namelist /calendar/ include_fleapyear
  
   ! *** machine ***
  integer                :: system_arch = 1    ! XD1 2(byte), HLRN 1(word)
  integer                :: n_levels = 1       ! Number of levels for hierarchic partitioning
  integer, dimension(10) :: n_part = RESHAPE((/0/), (/10/), (/0/)) ! Number of partitions on each hierarchy level
  namelist /machine/ system_arch, n_levels, n_part
  
  ! *** configuration***
  logical                       :: use_sw_pene=.true.
  logical                       :: use_ice=.false.  
  logical 						:: use_floatice = .false.
  logical                       :: toy_ocean=.false. ! Ersatz forcing has
                                                   ! to be supplied
  namelist /run_config/ use_ice,use_floatice, use_sw_pene, toy_ocean
  
  ! *** others ***
  real(kind=8)             	:: dt
  integer                       :: save_count_mean, save_count_restart
  logical                       :: r_restart
end module g_config

