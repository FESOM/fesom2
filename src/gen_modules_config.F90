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
  character(100)         :: ClimateDataPath='./hydrography/'
  character(100)         :: TideForcingPath='./tide_forcing/'
  character(1000)        :: ResultPath='./result/'
  namelist /paths/  MeshPath, ClimateDataPath, &
       TideForcingPath, ResultPath

  ! *** restart_log ***
  integer                :: logfile_outfreq=1      ! logfile info. outp. freq., # steps
  integer                :: restart_length=1
  character              :: restart_length_unit='m'
  
  namelist /restart_log/   restart_length, restart_length_unit, logfile_outfreq

  ! *** ale_def ***
  character(20)          :: which_ALE='linfs' ! 'zlevel', 'zstar', 'zstar-weighted', 'ztilde'
  logical                :: use_partial_cell=.false.  ! use partial bottom cell configuration  
  ! for zlevel: layer thickness should not become smaller than min_hnode of 
  ! original layer thickness. If it happens switch from zelvel to local zstar
  real(kind=WP)          :: min_hnode=0.5 
  ! for zlevel: in case min_hnode criteria is reached over how many level should 
  ! ssh change be distributed
  integer                :: lzstar_lev=4
  ! maximal pressure from ice felt by the ocean
  real(kind=WP)          :: max_ice_loading=5.0

  namelist /ale_def/ which_ALE, use_partial_cell, min_hnode, lzstar_lev, max_ice_loading

  ! *** model geometry ***
  logical                :: cartesian=.false.
  logical                :: fplane=.false.
  real(kind=WP)          :: cyclic_length=360. ! [degree]
  logical                :: rotated_grid=.true. ! not used
  logical                :: force_rotation=.true.
  real(kind=WP)          :: alphaEuler=50. 	! [degree] Euler angles, convention:
  real(kind=WP)          :: betaEuler=15.  	! first around z, then around new x,
  real(kind=WP)          :: gammaEuler=-90.	! then around new z.
  				                ! Set to zeros to work with
						! geographical coordinates
  namelist /geometry/  cartesian, fplane, &
       cyclic_length, rotated_grid, alphaEuler, betaEuler, gammaEuler, force_rotation

  ! *** fleap_year ***
  logical                       :: include_fleapyear=.false.
  namelist /calendar/ include_fleapyear
  
   ! *** machine ***
  integer                :: n_levels = 1       ! Number of levels for hierarchic partitioning
  integer, dimension(10) :: n_part = RESHAPE((/0/), (/10/), (/0/)) ! Number of partitions on each hierarchy level
  namelist /machine/ n_levels, n_part
  
  ! *** configuration***
  logical                       :: use_sw_pene=.true.
  logical                       :: use_ice=.false.  
  logical 						:: use_floatice = .false.
  logical                       :: toy_ocean=.false. ! Ersatz forcing has to be supplied
  character(100)                :: which_toy="soufflet" 
  logical                       :: flag_debug=.false.
  namelist /run_config/ use_ice,use_floatice, use_sw_pene, toy_ocean, which_toy, flag_debug
  
  ! *** others ***
  real(kind=WP)            	:: dt
  integer                       :: save_count_mean, save_count_restart
  logical                       :: r_restart
  real(kind=WP)             	:: rtime_ice=0.0, rtime_tot=0.0
  real(kind=WP)                 :: rtime_oce=0.0, rtime_oce_dyn=0.0, rtime_oce_dynssh=0.0,  rtime_oce_solvessh=0.0
  real(kind=WP)                 :: rtime_oce_solvetra=0.0, rtime_oce_GMRedi=0.0, rtime_oce_mixpres=0.0
  real(kind=WP)             	:: dummy=1.e10
  
  
end module g_config

