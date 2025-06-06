
! Modules of cell-vertex ocean model
! S. Danilov, 2012 (sergey.danilov@awi.de)
! SI units are used

!==========================================================
MODULE o_PARAM
integer, parameter            :: WP=8        ! Working precision
integer, parameter            :: MAX_PATH=4096 ! Maximum file path length
integer		                  :: mstep
real(kind=WP), parameter      :: pi=3.14159265358979
real(kind=WP), parameter      :: rad=pi/180.0_WP
real(kind=WP), parameter      :: density_0=1030.0_WP
real(kind=WP), parameter      :: density_0_r=1.0_WP/density_0 ! [m^3/kg]
real(kind=WP), parameter      :: g=9.81_WP
real(kind=WP), parameter      :: r_earth=6367500.0_WP
real(kind=WP), parameter      :: omega=2*pi/(3600.0_WP*24.0_WP)
real(kind=WP), parameter      :: vcpw=4.2e6   ![J/m^3/K] water heat cap
real(kind=WP), parameter      :: inv_vcpw = 1._WP / vcpw  ! inverse, to replace divide by multiply
real(kind=WP), parameter      :: small=1.0e-8 !small value
integer                       :: state_equation = 1     !1 - full equation of state, 0 - linear equation of state

real(kind=WP)                 :: C_d= 0.0025_WP ! Bottom drag coefficient
real(kind=WP)	              :: kappa=0.4      !von Karman's constant
real(kind=WP)                 :: mix_coeff_PP=0.01_WP   ! mixing coef for PP scheme
real(kind=WP)                 :: A_ver=0.001_WP ! Vertical harm. visc.
real(kind=WP)                 :: K_hor=10._WP
real(kind=WP)                 :: K_ver=0.00001_WP
real(kind=WP)                 :: scale_area=2.0e8
real(kind=WP)                 :: surf_relax_T= 0.0_WP
real(kind=WP)                 :: surf_relax_S= 10.0_WP/(60*3600.0_WP*24)
logical                       :: balance_salt_water =.true.
real(kind=WP)                 :: clim_relax= 1.0_WP/(10*3600.0_WP*24)
real(kind=WP)                 :: clim_decay, clim_growth
                                 ! set to 0.0 if no relaxation
logical                       :: ref_sss_local=.false.
real(kind=WP)                 :: ref_sss=34.7
logical                       :: Fer_GM =.false.  !flag for Ferrari et al. (2010) GM scheme
real(kind=WP)                 :: K_GM_max = 3000.
real(kind=WP)                 :: K_GM_min = 2.0
integer                       :: K_GM_bvref = 2 ! 0...surface, 1...bottom mixlay, 2...mean over mixlay
real(kind=WP)                 :: K_GM_resscalorder = 2.0
real(kind=WP)                 :: K_GM_rampmax = 40.0 ! Resol >K_GM_rampmax[km] GM full
real(kind=WP)                 :: K_GM_rampmin = 30.0 ! Resol <K_GM_rampmin[km] GM off
logical                       :: scaling_Ferreira   =.true.
logical                       :: scaling_Rossby     =.false.
logical                       :: scaling_resolution =.true.
logical                       :: scaling_FESOM14    =.false.
logical                       :: Redi               =.false.  !flag for Redi scheme

real(kind=WP)                 :: visc_sh_limit=5.0e-3      !for KPP, max visc due to shear instability
real(kind=WP)                 :: diff_sh_limit=5.0e-3      !for KPP, max diff due to shear instability
logical                       :: Kv0_const=.true.		    !use Kv0 varying with depth and latitude
logical                       :: double_diffusion=.false.  !for KPP,dd switch
                                 ! KPP parametrization
character(25)                 :: mix_scheme     ='KPP'	   !'KPP','PP'
integer                       :: mix_scheme_nmb = 1       ! choosen in oce_setup_step --> replace string by int comparison
real(KIND=WP)                 :: Ricr   = 0.3_WP  ! critical bulk Richardson Number
real(KIND=WP)                 :: concv  = 1.6_WP  ! constant for pure convection (eqn. 23) (Large 1.5-1.6; MOM default 1.8)

logical                       :: hbl_diag =.false.        ! writen boundary layer depth
logical                       :: use_global_tides=.false. ! tidal potential will be computed and used in the SSH gradient computation
! Time stepping
! real(kind=WP)                 :: alpha=1.0_WP, theta=1.0_WP ! implicitness for
real(kind=WP)                 :: alpha=1.0_WP, theta=1.0_WP ! implicitness for
                                                 ! elevation and divergence
real(kind=WP)                 :: epsilon=0.1_WP  ! AB2 offset
! Tracers

logical                       :: SPP=.false.

integer                       :: acc_vl = 64

TYPE tracer_source3d_type
    integer                             :: locID
    integer                             :: ID
    integer, allocatable, dimension(:)  :: ind2
END TYPE tracer_source3d_type

type(tracer_source3d_type), &
    allocatable, dimension(:)  :: ptracers_restore
integer                        :: ptracers_restore_total=0

!---wiso-code
! add water isotope parameters
real(kind=WP), dimension(3)    :: wiso_smow = (/2005.2e-6_WP, 155.76e-6_WP, 1.0_WP/)  ! water isotope SMOW values
integer, dimension(3)          :: index_wiso_tracers = (/-1, -1, -1/)  ! water isotope index in all tracers
!---wiso-code-end
!---age-code-begin
integer                        :: index_age_tracer = -1 ! water age tracer index in all tracers
!---age-code-end

! Momentum
!!PS logical                       :: free_slip=.false.
!!PS                                 ! false=no slip
!!PS integer                       :: mom_adv=2
                                ! 1 vector control volumes, p1 velocities
				! 2 scalar control volumes
				! 3 vector invariant

logical                       :: open_b=.false.   ! Reserved

!_______________________________________________________________________________
!--> mixing enhancement than can be applied via subroutine mo_convect(mesh)
!    additionally to every mixing scheme i.e. KPP, PP, cvmix_KPP, cvmix_PP, cvmix_TKE

! Switch for Monin-Obukov TB04 mixing --> can be additionally applied for all mixing schemes
! --> definetely recommented for KPP
logical                       :: use_momix     = .true. !.false. !Monin-Obukhov -> TB04 mixing on/off
real(kind=WP)                 :: momix_lat     = -50.0_WP ! latitudinal treshhold to apply mo_on <lat
real(kind=WP)                 :: momix_kv      = 0.01   ! for PP/KPP, mixing coefficient within MO length

! Switch for enhanced vertical mixing in case of instable stratification --> enhanced
! convection
logical                       :: use_instabmix = .true.
real(kind=WP)                 :: instabmix_kv  = 0.1

! Switch for enhanced wind mixing --> nasty trick from pp mixing in FESOM1.4
logical                       :: use_windmix   = .false.
real(kind=WP)                 :: windmix_kv    = 1.e-3
integer                       :: windmix_nl    = 2

!_______________________________________________________________________________
! use non-constant reference density if .false. density_ref=density_0
logical                       :: use_density_ref   = .false.
real(kind=WP)                 :: density_ref_T     = 2.0_WP
real(kind=WP)                 :: density_ref_S     = 34.0_WP

!_______________________________________________________________________________
! use k-profile nonlocal fluxes
logical                       :: use_kpp_nonlclflx = .false.

!_______________________________________________________________________________
! *** active tracer cutoff
logical          :: limit_salinity=.true.         !set an allowed range for salinity
real(kind=WP)    :: salinity_min=5.0              !minimal salinity
real(kind=WP)    :: coeff_limit_salinity=0.0023   !m/s, coefficient to restore s to s_min

  namelist /tracer_cutoff/ limit_salinity, salinity_min, coeff_limit_salinity

! *** others ***
 real(kind=WP)                        :: time_sum=0.0 ! for runtime estimate

!___________________________________________
! Pressure Gradient Force  calculation (pgf)
! calculation of pgf either:
! only linfs:
! > 'nemo'         ... like NEMO (interpolate to elemental depth, inter-/extrapolation)
! linfs, zlevel, zstar:
! > 'shchepetkin'  ... based on density jacobian
! > 'cubicspline'  ... like in FESOM1.4
! > 'easypgf'      ... interpolate pressure on elemental depth
character(20)                  :: which_pgf='shchepetkin'


 NAMELIST /oce_dyn/ state_equation, C_d, A_ver, &
                    scale_area, SPP,&
                    Fer_GM, K_GM_max, K_GM_min, K_GM_bvref, K_GM_resscalorder, K_GM_rampmax, K_GM_rampmin, &
                    scaling_Ferreira, scaling_Rossby, scaling_resolution, scaling_FESOM14, &
                    Redi, visc_sh_limit, mix_scheme, Ricr, concv, which_pgf, alpha, theta, use_density_ref

 NAMELIST /tracer_phys/ diff_sh_limit, Kv0_const, double_diffusion, K_ver, K_hor, surf_relax_T, surf_relax_S, &
            balance_salt_water, clim_relax, ref_sss_local, ref_sss, &
            use_momix, momix_lat, momix_kv, &
            use_instabmix, instabmix_kv, &
            use_windmix, windmix_kv, windmix_nl, &
            use_kpp_nonlclflx

END MODULE o_PARAM
!==========================================================
MODULE o_ARRAYS
USE o_PARAM
IMPLICIT NONE
! Arrays are described in subroutine array_setup
!real(kind=WP), allocatable         :: UV_ib(:,:,:) ! kh 08.03.21 additional array for asynchronous iceberg computations
real(kind=WP), allocatable         :: hpressure(:,:)
real(kind=WP), allocatable         :: stress_surf(:,:)
real(kind=WP), allocatable         :: stress_node_surf(:,:)
REAL(kind=WP), ALLOCATABLE         :: stress_atmoce_x(:)
REAL(kind=WP), ALLOCATABLE         :: stress_atmoce_y(:)
real(kind=WP), allocatable         :: heat_flux(:), Tsurf(:) 
real(kind=WP), allocatable         :: heat_flux_in(:) !to keep the unmodified (by SW penetration etc.) heat flux 
real(kind=WP), allocatable         :: Tsurf_ib(:) ! kh 15.03.21 additional array for asynchronous iceberg computations
real(kind=WP), allocatable    :: water_flux(:), Ssurf(:)
real(kind=WP), allocatable    :: Ssurf_ib(:) ! kh 15.03.21 additional array for asynchronous iceberg computations
real(kind=WP), allocatable    :: virtual_salt(:), relax_salt(:)
real(kind=WP), allocatable    :: Tclim(:,:), Sclim(:,:)
real(kind=WP), allocatable    :: t_star(:), qsr_c(:)

!--------------
! LA: add iceberg tracer arrays 2023-02-08
!--------------
real(kind=WP), allocatable    :: Tclim_ib(:,:), Sclim_ib(:,:)
!!PS real(kind=WP), allocatable    :: Visc(:,:)
real(kind=WP), allocatable    :: Tsurf_t(:,:), Ssurf_t(:,:)
real(kind=WP), allocatable    :: tau_x_t(:,:), tau_y_t(:,:)
real(kind=WP), allocatable    :: heat_flux_t(:,:), heat_rel_t(:,:), heat_rel(:)
!!PS real(kind=WP), allocatable    :: coriolis(:), coriolis_node(:)
real(kind=WP), allocatable    :: relax2clim(:)
real(kind=WP), allocatable    :: MLD1(:), MLD2(:), MLD3(:)
integer,       allocatable    :: MLD1_ind(:), MLD2_ind(:), MLD3_ind(:)
real(kind=WP), allocatable    :: ssh_gp(:)
!Tracer gradients&RHS
real(kind=WP), allocatable :: tr_xy(:,:,:)
real(kind=WP), allocatable :: tr_z(:,:)

#if defined(__recom)
real(kind=WP), allocatable    :: dtr_bf(:,:), str_bf(:,:)
real(kind=WP), allocatable    :: vert_sink(:,:)
#endif

!Viscosity and diff coefs
real(kind=WP), allocatable,dimension(:,:)   :: Av,Kv
real(kind=WP), allocatable,dimension(:,:,:) :: Kv_double
real(kind=WP), allocatable,dimension(:)     :: Kv0
!Velocities interpolated to nodes
!!PS real(kind=WP), allocatable,dimension(:,:,:)   :: Unode

! Auxiliary arrays to store Redi-GM fields
real(kind=WP), allocatable,dimension(:,:,:) :: neutral_slope
real(kind=WP), allocatable,dimension(:,:,:) :: slope_tapered
real(kind=WP), allocatable,dimension(:,:,:) :: sigma_xy
real(kind=WP), allocatable,dimension(:,:)   :: sw_beta, sw_alpha
real(kind=WP), allocatable,dimension(:)     :: dens_flux
!real(kind=WP), allocatable,dimension(:,:,:) :: tsh, tsv, tsh_nodes
!real(kind=WP), allocatable,dimension(:,:)   :: hd_flux,vd_flux
!Isoneutral diffusivities (or xy diffusivities if Redi=.false)
real(kind=WP), allocatable :: Ki(:,:)

! --> auxiliary array to store an intermediate part of the rhs computations.
real(kind=WP), allocatable,dimension(:)     :: ssh_rhs_old !, ssh_rhs_old2 !PS
real(kind=WP)                               :: is_nonlinfs

!_______________________________________________________________________________
! Arrays added for pressure gradient force calculation
real(kind=WP), allocatable,dimension(:,:)   :: density_m_rho0
real(kind=WP), allocatable,dimension(:,:)   :: density_m_rho0_slev
real(kind=WP), allocatable,dimension(:,:)   :: density_ref
real(kind=WP), allocatable,dimension(:,:)   :: density_dmoc
real(kind=WP), allocatable,dimension(:,:)   :: pgf_x, pgf_y

!_______________________________________________________________________________
!!PS ! dummy arrays
real(kind=WP), allocatable,dimension(:,:)   :: dum_3d_n !, dum_3d_e
!!PS real(kind=WP), allocatable,dimension(:)     :: dum_2d_n, dum_2d_e

!_______________________________________________________________________________
!Monin-Obukhov correction
real(kind=WP),allocatable :: mo(:,:),mixlength(:)
!GM_stuff
real(kind=WP),allocatable :: bvfreq(:,:),mixlay_dep(:),bv_ref(:)

real(kind=WP), target, allocatable    :: fer_c(:), fer_scal(:), fer_K(:,:), fer_gamma(:,:,:)

real(kind=WP),         allocatable    :: ice_rejected_salt(:)

!---wiso-code
real(kind=WP), allocatable    :: tr_arr_ice(:,:)      !---wiso-code: add sea ice isotope tracers
real(kind=WP), allocatable    :: wiso_flux_oce(:,:)   !---wiso-code: add isotope fluxes over open water
real(kind=WP), allocatable    :: wiso_flux_ice(:,:)   !---wiso-code: add isotope fluxes over sea ice
!---wiso-code-end

END MODULE o_ARRAYS
!==========================================================
