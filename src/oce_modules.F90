
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
real(kind=WP)                 :: gamma0=0.01! [m/s], gamma0*len*dt is the background viscosity
real(kind=WP)                 :: gamma1=0.1!  [non dim.], or computation of the flow aware viscosity
real(kind=WP)                 :: gamma2=10.!  [s/m],      is only used in easy backscatter option
real(kind=WP)                 :: Div_c  =1.0_WP !modified Leith viscosity weight
real(kind=WP)                 :: Leith_c=1.0_WP	!Leith viscosity weight. It needs vorticity!
real(kind=WP)                 :: easy_bs_return=1.0 !backscatter option only (how much to return)
real(kind=WP)                 :: A_ver=0.001_WP ! Vertical harm. visc.
integer                       :: visc_option=5
logical                       :: uke_scaling=.true.
real(kind=WP)                 :: uke_scaling_factor=1._WP
real(kind=WP)		          :: rosb_dis=1._WP
integer                       :: smooth_back=2
integer                       :: smooth_dis=2
integer                       :: smooth_back_tend=4
real(kind=WP)		          :: K_back=600._WP
real(kind=WP)                 :: c_back=0.1_8
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
logical                       :: i_vert_diff= .true.
logical                       :: i_vert_visc= .true.
character(20)                 :: tra_adv_ver, tra_adv_hor, tra_adv_lim
real(kind=WP)                 :: tra_adv_ph, tra_adv_pv
logical                       :: w_split  =.false.
real(kind=WP)                 :: w_max_cfl=1.e-5_WP

logical                       :: SPP=.false.

TYPE tracer_source3d_type
    integer                             :: locID
    integer                             :: ID
    integer, allocatable, dimension(:)  :: ind2
END TYPE tracer_source3d_type

integer	                       :: num_tracers=2
integer, dimension(100)        :: tracer_ID  = RESHAPE((/0, 1/), (/100/), (/0/)) ! ID for each tracer for treating the initialization and surface boundary condition
                                                                                 ! 0=temp, 1=salt etc.
type(tracer_source3d_type), &
    allocatable, dimension(:)  :: ptracers_restore
integer                        :: ptracers_restore_total=0


! Momentum
logical                       :: free_slip=.false.
                                ! false=no slip 
integer                       :: mom_adv=2
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

! bharmonic diffusion for tracers. We recommend to use this option in very high resolution runs (Redi is generally off there).
logical                       :: smooth_bh_tra = .false.
real(kind=WP)                 :: gamma0_tra    = 0.0005
real(kind=WP)                 :: gamma1_tra    = 0.0125
real(kind=WP)                 :: gamma2_tra    = 0.
!_______________________________________________________________________________
! use non-constant reference density if .false. density_ref=density_0
logical                       :: use_density_ref   = .false.
real(kind=WP)                 :: density_ref_T     = 2.0_WP
real(kind=WP)                 :: density_ref_S     = 34.0_WP

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


 NAMELIST /oce_dyn/ state_equation, C_d, A_ver, gamma0, gamma1, gamma2, Leith_c, Div_c, easy_bs_return, &
                    scale_area, mom_adv, free_slip, i_vert_visc, w_split, w_max_cfl, SPP,&
                    Fer_GM, K_GM_max, K_GM_min, K_GM_bvref, K_GM_resscalorder, K_GM_rampmax, K_GM_rampmin, & 
                    scaling_Ferreira, scaling_Rossby, scaling_resolution, scaling_FESOM14, & 
                    Redi, visc_sh_limit, mix_scheme, Ricr, concv, which_pgf, visc_option, alpha, theta, use_density_ref, &
                    K_back, c_back, uke_scaling, uke_scaling_factor, smooth_back, smooth_dis, &
                    smooth_back_tend, rosb_dis

 NAMELIST /oce_tra/ diff_sh_limit, Kv0_const, double_diffusion, K_ver, K_hor, surf_relax_T, surf_relax_S, &
            balance_salt_water, clim_relax, ref_sss_local, ref_sss, i_vert_diff, tra_adv_ver, tra_adv_hor, &
            tra_adv_lim, tra_adv_ph, tra_adv_pv, num_tracers, tracer_ID, &
            use_momix, momix_lat, momix_kv, &
            use_instabmix, instabmix_kv, &
            use_windmix, windmix_kv, windmix_nl, &
            smooth_bh_tra, gamma0_tra, gamma1_tra, gamma2_tra
            
END MODULE o_PARAM  
!==========================================================

!==========================================================
MODULE o_MESH
USE o_PARAM
USE, intrinsic :: ISO_FORTRAN_ENV
! All variables used to keep the mesh structure +
! auxiliary variables involved in implementation 
! of open boundaries and advection schemes
!
! The fct part
real(kind=WP),allocatable,dimension(:,:)      :: fct_LO          ! Low-order solution
real(kind=WP),allocatable,dimension(:,:)      :: adv_flux_hor    ! Antidif. horiz. contrib. from edges / backup for iterafive fct scheme
real(kind=WP),allocatable,dimension(:,:)      :: adv_flux_ver    ! Antidif. vert. fluxes from nodes    / backup for iterafive fct scheme

real(kind=WP),allocatable,dimension(:,:)      :: fct_ttf_max,fct_ttf_min
real(kind=WP),allocatable,dimension(:,:)      :: fct_plus,fct_minus
! Quadratic reconstruction part
integer,allocatable,dimension(:)              :: nn_num, nboundary_lay
real(kind=WP),allocatable,dimension(:,:,:)    :: quad_int_mat, quad_int_coef
integer,allocatable,dimension(:,:)            :: nn_pos
! MUSCL type reconstruction
integer,allocatable,dimension(:,:)            :: edge_up_dn_tri
real(kind=WP),allocatable,dimension(:,:,:)    :: edge_up_dn_grad
end module o_MESH
!==========================================================

!==========================================================
MODULE o_ARRAYS
USE o_PARAM
IMPLICIT NONE
! Arrays are described in subroutine array_setup  
real(kind=WP), allocatable, target :: Wvel(:,:), Wvel_e(:,:), Wvel_i(:,:)
real(kind=WP), allocatable         :: UV(:,:,:)
real(kind=WP), allocatable         :: UV_rhs(:,:,:), UV_rhsAB(:,:,:)
real(kind=WP), allocatable         :: uke(:,:), v_back(:,:), uke_back(:,:), uke_dis(:,:), uke_dif(:,:) 
real(kind=WP), allocatable         :: uke_rhs(:,:), uke_rhs_old(:,:)
real(kind=WP), allocatable         :: UV_dis_tend(:,:,:), UV_back_tend(:,:,:), UV_total_tend(:,:,:), UV_dis_tend_node(:,:,:)
real(kind=WP), allocatable         :: UV_dis_posdef_b2(:,:), UV_dis_posdef(:,:), UV_back_posdef(:,:)
real(kind=WP), allocatable         :: eta_n(:), d_eta(:)
real(kind=WP), allocatable         :: ssh_rhs(:), hpressure(:,:)
real(kind=WP), allocatable         :: CFL_z(:,:)
real(kind=WP), allocatable         :: stress_surf(:,:)
real(kind=WP), allocatable         :: stress_node_surf(:,:)
REAL(kind=WP), ALLOCATABLE         :: stress_atmoce_x(:)
REAL(kind=WP), ALLOCATABLE         :: stress_atmoce_y(:)
real(kind=WP), allocatable         :: T_rhs(:,:) 
real(kind=WP), allocatable         :: heat_flux(:), Tsurf(:) 
real(kind=WP), allocatable         :: heat_flux_in(:) !to keep the unmodified (by SW penetration etc.) heat flux 
real(kind=WP), allocatable         :: S_rhs(:,:)
real(kind=WP), allocatable         :: tr_arr(:,:,:),tr_arr_old(:,:,:)
real(kind=WP), allocatable         :: del_ttf(:,:)
real(kind=WP), allocatable         :: del_ttf_advhoriz(:,:),del_ttf_advvert(:,:) !!PS ,del_ttf_diff(:,:)

real(kind=WP), allocatable    :: dtr_bf(:,:), str_bf(:,:) ! OG, jh -> for fluxes from benthos
real(kind=WP), allocatable    :: vert_sink(:,:)   ! OG -> vertical sinking
real(kind=WP), allocatable    :: nss(:,:)   ! OG -> vertical sinking

real(kind=WP), allocatable    :: water_flux(:), Ssurf(:)
real(kind=WP), allocatable    :: virtual_salt(:), relax_salt(:)
real(kind=WP), allocatable    :: Tclim(:,:), Sclim(:,:)
real(kind=WP), allocatable    :: Visc(:,:)
real(kind=WP), allocatable    :: Tsurf_t(:,:), Ssurf_t(:,:)
real(kind=WP), allocatable    :: tau_x_t(:,:), tau_y_t(:,:) 
real(kind=WP), allocatable    :: heat_flux_t(:,:), heat_rel_t(:,:), heat_rel(:) 
real(kind=WP), allocatable    :: coriolis(:), coriolis_node(:)
real(kind=WP), allocatable    :: relax2clim(:)
real(kind=WP), allocatable    :: MLD1(:), MLD2(:)
integer,       allocatable    :: MLD1_ind(:), MLD2_ind(:)
real(kind=WP), allocatable    :: ssh_gp(:)
! Passive and age tracers
real(kind=WP), allocatable    :: tracer(:,:,:), tracer_rhs(:,:,:)   
!Tracer gradients&RHS      
real(kind=WP), allocatable :: ttrhs(:,:)
real(kind=WP), allocatable :: tr_xy(:,:,:)
real(kind=WP), allocatable :: tr_z(:,:)

! Auxiliary arrays for vector-invariant form of momentum advection
real(kind=WP), allocatable,dimension(:,:)   :: vorticity

!Viscosity and diff coefs
real(kind=WP), allocatable,dimension(:,:)   :: Av,Kv
real(kind=WP), allocatable,dimension(:,:,:) :: Kv_double
real(kind=WP), allocatable,dimension(:)     :: Kv0
!Velocities interpolated to nodes
real(kind=WP), allocatable,dimension(:,:,:)   :: Unode

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

!_______________________________________________________________________________
! Arrays added for ALE implementation:
! --> layer thinkness at node and depthlayer for t=n and t=n+1
real(kind=WP), allocatable,dimension(:,:)   :: hnode, hnode_new, zbar_3d_n, Z_3d_n

! --> layer thinkness at elements, interpolated from hnode
real(kind=WP), allocatable,dimension(:,:)   :: helem

! --> thinkness of bottom elem (important for partial cells)
real(kind=WP), allocatable,dimension(:)     :: bottom_elem_thickness 
real(kind=WP), allocatable,dimension(:)     :: bottom_node_thickness 

! --> The increment of total fluid depth on elements. It is used to update the matrix
real(kind=WP), allocatable,dimension(:)     :: dhe

! --> hbar, hbar_old: correspond to the elevation, but on semi-integer time steps.
real(kind=WP), allocatable,dimension(:)     :: hbar, hbar_old

! --> auxiliary array to store an intermediate part of the rhs computations.
real(kind=WP), allocatable,dimension(:)     :: ssh_rhs_old !, ssh_rhs_old2 !PS

! --> auxiliary array to store depth of layers and depth of mid level due to changing 
!     layer thinkness at every node
real(kind=WP), allocatable,dimension(:)     :: zbar_n, Z_n

! new bottom depth at node and element due to partial cells
real(kind=WP), allocatable,dimension(:)     :: zbar_n_bot
real(kind=WP), allocatable,dimension(:)     :: zbar_e_bot

! new depth of cavity-ocean interface at node and element due to partial cells
real(kind=WP), allocatable,dimension(:)     :: zbar_n_srf
real(kind=WP), allocatable,dimension(:)     :: zbar_e_srf

! --> multiplication factor for surface boundary condition in 
!     diff_ver_part_impl_ale(tr_num) between linfs -->=0.0 and noninfs 
!     (zlevel,zstar...) --> = 1.0
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

real(kind=WP),         allocatable    :: fer_UV(:,:,:), fer_wvel(:,:)
real(kind=WP), target, allocatable    :: fer_c(:), fer_scal(:), fer_K(:,:), fer_gamma(:,:,:)

real(kind=WP),         allocatable    :: ice_rejected_salt(:)

!_______________________________________________________________________________
! in case ldiag_DVD=.true. --> calculate discrete variance decay (DVD)
real(kind=WP), allocatable    :: tr_dvd_horiz(:,:,:),tr_dvd_vert(:,:,:)

END MODULE o_ARRAYS
!==========================================================
