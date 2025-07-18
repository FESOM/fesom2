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
  !_____________________________________________________________________________
  ! *** Modelname ***
  character(10)           :: runid='test1'       ! a model/setup name
  namelist /modelname/ runid
  
  !_____________________________________________________________________________
  ! *** time step ***
  integer                :: step_per_day=72     ! number of steps per day
  integer                :: run_length=1	! run length
  character              :: run_length_unit='y' ! unit: y, d, s
  namelist /timestep/ step_per_day, run_length, run_length_unit
  
  !_____________________________________________________________________________
  ! *** Paths for all in and out ***
! kh 01.03.21 paths in test environments can easily become longer than 100 characters (the former value) 
  character(MAX_PATH)        :: MeshPath='./mesh/'
  character(MAX_PATH)        :: ClimateDataPath='./hydrography/'
  character(MAX_PATH)        :: TideForcingPath='./tide_forcing/'
  character(MAX_PATH)        :: ResultPath='./result/'
  character(20)              :: MeshId='NONE'
  namelist /paths/  MeshPath, ClimateDataPath, &
       TideForcingPath, ResultPath, MeshId
       
  !_____________________________________________________________________________
  ! *** restart_log ***
  integer                :: logfile_outfreq=1      ! logfile info. outp. freq., # steps
  integer                :: restart_length=1
  character(3)           :: restart_length_unit='m'
  integer                :: raw_restart_length=1
  character(3)           :: raw_restart_length_unit='m'
  integer                :: bin_restart_length=1
  character(3)           :: bin_restart_length_unit='m'
  
  namelist /restart_log/   restart_length    , restart_length_unit, & 
                           raw_restart_length, raw_restart_length_unit, &
                           bin_restart_length, bin_restart_length_unit, &
                           logfile_outfreq
  !_____________________________________________________________________________
  ! *** ale_def ***
  ! Which ALE case to use : 'linfs', 'zlevel', 'zstar'
  character(20)          :: which_ALE='linfs' 
  
  ! use partial bottom cell configuration .true./.false. 
  logical                :: use_partial_cell=.false.  
  
  ! if a thicker partial bottom layer thickness is more realistic than always 
  ! apply it, BUT when a thinner partial bottom layer thickness is more realistic than 
  ! only apply it when the initial full cell bottom layer thickness is above the 
  ! treshhold partial_cell_tresh to not allow already thin layers to become even 
  ! thinner. e.g 
  ! partial_cell_tresh=10 --> thinner partial bottom cells will be only applied for initial 
  ! bottom layer thicknesses larger than 10m 
  real(kind=WP)          :: partial_cell_thresh=0.0_WP 
  
  ! for zlevel: layer thickness should not become smaller than min_hnode of 
  ! original layer thickness. If it happens switch from zelvel to local zstar
  real(kind=WP)          :: min_hnode=0.5 
  
  ! for zlevel: in case min_hnode criteria is reached over how many level should 
  ! ssh change be distributed
  integer                :: lzstar_lev=4
  
  ! maximal pressure from ice felt by the ocean
  real(kind=WP)          :: max_ice_loading=5.0

  namelist /ale_def/ which_ALE, use_partial_cell, partial_cell_thresh, min_hnode, lzstar_lev, max_ice_loading

  !_____________________________________________________________________________
  ! *** model geometry ***
  logical                :: cartesian=.false.
  logical                :: fplane=.false.
  real(kind=WP)          :: cyclic_length=360. ! [degree]
  logical                :: rotated_grid=.true. ! not used
  logical                :: force_rotation=.true.
  real(kind=WP)          :: alphaEuler=50.  ! [degree] Euler angles, convention:
  real(kind=WP)          :: betaEuler=15.   ! first around z, then around new x,
  real(kind=WP)          :: gammaEuler=-90. ! then around new z.
                                            ! Set to zeros to work with
                                            ! geographical coordinates
  integer                :: thers_zbar_lev=5     ! minimum number of levels to be                                            
  character(len=5)       :: which_depth_n2e='mean'    
  
  logical                :: use_depthonelem =.false.
  character(len=10)      :: use_depthfile='aux3d'   ! 'aux3d', 'depth@'        
  logical                :: use_cavityonelem=.false.
  
  namelist /geometry/   cartesian, fplane, &
                        cyclic_length, rotated_grid, force_rotation, &
                        alphaEuler, betaEuler, gammaEuler, &
                        which_depth_n2e, use_depthonelem, use_cavityonelem, use_depthfile

  !_____________________________________________________________________________
  ! *** fleap_year ***
  logical                       :: include_fleapyear=.false.
  logical                       :: use_flpyrcheck   =.true.
  namelist /calendar/ include_fleapyear, use_flpyrcheck
  
  !_____________________________________________________________________________
  ! *** machine ***
  integer                :: n_levels = 1       ! Number of levels for hierarchic partitioning
  integer, dimension(10) :: n_part = RESHAPE((/0/), (/10/), (/0/)) ! Number of partitions on each hierarchy level
  namelist /machine/ n_levels, n_part
  
  !_____________________________________________________________________________
  ! *** configuration***
  logical                       :: use_sw_pene=.true.
  logical                       :: use_ice=.false.  
                                                   ! to be supplied
  ! *** icebergs ***
  logical                       :: use_icebergs=.false.
  logical                       :: turn_off_hf=.false.
  logical                       :: turn_off_fw=.false.
  logical                       :: use_icesheet_coupling=.false.  
  logical                       :: lbalance_fw=.true.
  integer                       :: cell_saturation=2 ! 0=no cell saturation, 1=one additional iceberg allowed, 2=no daddtional iceberg allowed
  logical                       :: lmin_latent_hf=.true.
  logical                       :: lverbose_icb=.false.  
  integer                       :: ib_num=0
  integer                       :: steps_per_ib_step=8

! kh 02.02.21
! ib_async_mode == 0: original sequential behavior for both ice sections (for testing purposes, creating reference results etc.)
! ib_async_mode == 1: OpenMP code active to overlapped computations in first (ocean ice) and second (icebergs) parallel section
! ib_async_mode == 2: OpenMP code active but computations still serialized via spinlock (for testing purposes)
  integer                       :: ib_async_mode=0
  integer                       :: thread_support_level_required=3 ! 2 = MPI_THREAD_SERIALIZED, 3 = MPI_THREAD_MULTIPLE

  namelist /icebergs/   use_icebergs, turn_off_hf, turn_off_fw, use_icesheet_coupling, lbalance_fw, cell_saturation, lmin_latent_hf, &
                        ib_num, steps_per_ib_step, ib_async_mode, thread_support_level_required, lverbose_icb

!wiso-code!!!
  logical                       :: lwiso  =.false.  ! enable isotope?
!wiso-code!!!
  logical                       :: use_floatice = .false.
  logical                       :: use_cavity = .false. ! switch on/off cavity usage
  logical                       :: use_cavity_partial_cell = .false. ! switch on/off cavity usage
  logical                       :: use_cavity_fw2press = .true. ! switch on/off cavity+zstar input of freshwater leads to increase in pressure
  real(kind=WP)                 :: cavity_partial_cell_thresh=0.0_WP ! same as partial_cell_tresh but for surface
  logical                       :: toy_ocean=.false. ! Ersatz forcing has to be supplied
  character(100)                :: which_toy="soufflet" 
  logical                       :: flag_debug=.false.    ! prints name of actual subroutine he is in 
  logical                       :: flag_warn_cflz=.true. ! switches off cflz warning
  logical                       :: use_transit=.false.    ! switches off transient tracers
  logical                       :: compute_oasis_corners=.false. ! switches on corner calculation for 1st order conserv remapping 
  namelist /run_config/ use_ice,use_floatice, use_sw_pene, use_cavity, & 
                        use_cavity_partial_cell, cavity_partial_cell_thresh, &
                        use_cavity_fw2press, toy_ocean, which_toy, flag_debug, flag_warn_cflz, lwiso, &
                        use_transit, compute_oasis_corners
  
  !_____________________________________________________________________________
  ! *** others ***
  real(kind=WP)                 :: dt
  integer                       :: save_count_mean, save_count_restart
  logical                       :: r_restart
  real(kind=WP)                 :: rtime_ice=0.0, rtime_tot=0.0
  real(kind=WP)                 :: rtime_oce=0.0, rtime_oce_dyn=0.0, rtime_oce_dynssh=0.0,  rtime_oce_solvessh=0.0
  real(kind=WP)                 :: rtime_oce_solvetra=0.0, rtime_oce_GMRedi=0.0, rtime_oce_mixpres=0.0
  real(kind=WP)                 :: dummy=1.e10
  
  
end module g_config
