! Modules of cell-vertex ocean model
! S. Danilov, 2012 (sergey.danilov@awi.de)
! SI units are used
  
!==========================================================
MODULE o_PARAM
integer, parameter            :: WP=8        ! Working precision
integer		                  :: mstep
real(kind=WP), parameter      :: pi=3.14159265358979
real(kind=WP), parameter      :: rad=pi/180.0_WP
real(kind=WP), parameter      :: density_0=1030.0_WP
real(kind=WP), parameter      :: density_0_r=1.0_WP/density_0 ! [m^3/kg]         
real(kind=WP), parameter      :: g=9.81_WP
real(kind=WP), parameter      :: r_earth=6367500.0_WP
real(kind=WP), parameter      :: omega=2*pi/(3600.0_WP*24.0_WP)
real(kind=WP), parameter      :: vcpw=4.2e6   ![J/m^3/K] water heat cap
real(kind=WP), parameter      :: small=1.0e-8 !small value

real(kind=WP)                 :: C_d= 0.0025_WP ! Bottom drag coefficient
real(kind=WP)	              :: kappa=0.4      !von Karman's constant
real(kind=WP)                 :: mix_coeff_PP=0.01_WP   ! mixing coef for PP scheme
real(kind=WP)                 :: A_hor=100.0_WP		! Horizontal harm. visc.    
real(kind=WP)                 :: A_hor_max=1500.0_WP	! Maximum viscosity allowed (to limit Smag and Leith
							! contributions when they are too large
real(kind=WP)                 :: Leith_c=0.		! Leith viscosity. It needs vorticity, which is only computed for 
							! the vector invariant form of momentum advection (mom_adv=4)
real(kind=WP)                 :: tau_c=0.4		! Controls the strength of filters. Should be about 0.4
real(kind=WP)                 :: A_ver=0.001_WP ! Vertical harm. visc.
real(kind=WP)                 :: Div_c=0.5_WP
real(kind=WP)                 :: Smag_c=0.0_WP  ! 0.2   ! (C/pi)^2
real(kind=WP)                 :: Abh0=8.0e12    ! 
logical                       :: laplacian=.false.
logical                       :: biharmonic=.true.
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
real(kind=WP)                 :: K_GM=1000.
logical			      :: scaling_Ferreira   =.true.
logical			      :: scaling_Rossby     =.false.
logical			      :: scaling_resolution =.true.
logical			      :: scaling_FESOM14    =.false.
logical			      :: Redi               =.false.  !flag for Redi scheme

real(kind=WP)                 :: visc_sh_limit=5.0e-3      !for KPP, max visc due to shear instability
real(kind=WP)                 :: diff_sh_limit=5.0e-3      !for KPP, max diff due to shear instability
logical                       :: Kv0_const=.true.		    !use Kv0 varying with depth and latitude 
logical                       :: double_diffusion=.false.  !for KPP,dd switch
                                 ! KPP parametrization
 character(5)                 :: mix_scheme='KPP'	   !'KPP','PP'
real(KIND=WP)                 :: Ricr   = 0.3_WP  ! critical bulk Richardson Number
real(KIND=WP)                 :: concv  = 1.6_WP  ! constant for pure convection (eqn. 23) (Large 1.5-1.6; MOM default 1.8)

logical                       :: hbl_diag =.false.   ! writen boundary layer depth

! Time stepping                               
! real(kind=WP)                 :: alpha=1.0_WP, theta=1.0_WP ! implicitness for
real(kind=WP)                 :: alpha=1.0_WP, theta=1.0_WP ! implicitness for
                                                 ! elevation and divergence
real(kind=WP)                 :: epsilon=0.1_WP  ! AB2 offset 
! Tracers
logical                       :: i_vert_diff= .true.
logical                       :: i_vert_visc= .true.
integer                       :: tracer_adv=1

logical                       :: w_split  =.false.
real(kind=WP)                 :: w_exp_max=1.e-5_WP

				! 1 MUSCL
				! 2 MUSCL-FCT
integer	                       :: num_tracers=2
! Momentum
logical                       :: free_slip=.false.
                                ! false=no slip 
integer                       :: mom_adv=2
                                ! 1 vector control volumes, p1 velocities
				! 2 scalar control volumes  
				! 3 vector invariant 

logical                       :: open_b=.false.   ! Reserved    

logical                       :: mo_on=.false. !Monin-Obukhov
real*8 :: modiff=0.01                   !for PP, mixing coefficient within MO length

  ! *** active tracer cutoff
logical          :: limit_salinity=.true.         !set an allowed range for salinity
real(kind=WP)    :: salinity_min=5.0              !minimal salinity 
real(kind=WP)    :: coeff_limit_salinity=0.0023   !m/s, coefficient to restore s to s_min

  namelist /tracer_cutoff/ limit_salinity, salinity_min, coeff_limit_salinity

! *** others ***
 integer                       :: num_tracer
 real*8                        :: time_sum=0.0 ! for runtime estimate


 NAMELIST /oce_dyn/ C_d, A_ver, laplacian, A_hor, A_hor_max, Leith_c, tau_c, Div_c, Smag_c, &
                    biharmonic, Abh0, scale_area, mom_adv, free_slip, i_vert_visc, w_split, w_exp_max, &
                    Fer_GM, K_GM, scaling_Ferreira, scaling_Rossby, scaling_resolution, scaling_FESOM14, Redi, visc_sh_limit, mix_scheme, Ricr, concv
 NAMELIST /oce_tra/ diff_sh_limit, Kv0_const, double_diffusion, K_ver, K_hor, surf_relax_T, surf_relax_S, balance_salt_water, clim_relax, &
		    ref_sss_local, ref_sss, i_vert_diff, &
		    tracer_adv
END MODULE o_PARAM  
!==========================================================

!==========================================================
MODULE o_MESH
USE o_PARAM
! All variables used to keep the mesh structure +
! auxiliary variables involved in implementation 
! of open boundaries and advection schemes
! 

integer, parameter                         :: MAX_ADJACENT=32 ! Max allowed number of adjacent nodes
integer                                    ::   nod2D      ! the number of 2D nodes
real(kind=WP)                              ::   ocean_area
real(kind=WP), allocatable, dimension(:,:) ::   coord_nod2D, geo_coord_nod2D
integer                                    ::   edge2D     ! the number of 2D edges
integer                                    ::   edge2D_in  
                                              ! the number of internal 2D edges
integer                                    ::   elem2D     ! the number of 2D elements
integer, allocatable, dimension(:,:)       ::   elem2D_nodes
                                              ! elem2D_nodes(:,n) lists
				              ! 3 nodes of element n   
integer, allocatable, dimension(:,:)       ::   edges
                                              ! edge(:,n) lists 2 nodes
				              ! edge n
integer, allocatable, dimension(:,:)       ::   edge_tri
                                              ! edge_tri(:,n) lists 2 
				              ! elements containing edge n
				              ! The first one is to left 
				              ! of the line directed
				              ! to the second node
integer, allocatable, dimension(:,:)       ::   elem_edges
                                              ! elem_edges(:,n) are edges of 
                                              ! element n.  
real(kind=WP), allocatable, dimension(:)   ::   elem_area
real(kind=WP), allocatable, dimension(:,:) ::   edge_dxdy, edge_cross_dxdy
real(kind=WP), allocatable, dimension(:)   ::   elem_cos, metric_factor
integer,allocatable,dimension(:,:)         ::   elem_neighbors
integer,allocatable,dimension(:,:)         ::   nod_in_elem2D
integer,allocatable,dimension(:)           ::   nod_in_elem2D_num
real(kind=WP),allocatable,dimension(:)     ::   depth
                                              ! depth(n) is the depths at 
				              ! node n 
real(kind=WP),allocatable,dimension(:,:)    ::   gradient_vec 
                                              ! Coefficients of linear reconstruction
					      ! of velocities on elements
real(kind=WP),allocatable,dimension(:,:)    ::   gradient_sca
                                              ! Coefficients to compute
					      ! gradient of scalars on elements
! Vertical structure             
integer                                    :: nl
real(kind=8), allocatable, dimension(:)    :: zbar, Z,elem_depth
integer, allocatable, dimension(:)         :: nlevels, nlevels_nod2D
real(kind=8), allocatable, dimension(:,:)  :: area, area_inv
real(kind=WP), allocatable, dimension(:)   :: mesh_resolution


  type sparse_matrix 
     integer :: nza
     integer :: dim
     real(kind=8), allocatable, dimension(:)      :: values
     integer(KIND=4), allocatable,   dimension(:) :: colind
     integer(KIND=4), allocatable,   dimension(:) :: rowptr
  end type sparse_matrix
! Elevation stiffness matrix
type(sparse_matrix)                           :: ssh_stiff

! Auxiliary arrays. They are not related to mesh structure, but are 
! kept here because they are just used for temporary storage in computations

! Open boundary:
integer                                       :: ob_num  ! number of OB fragments

TYPE ob_type
    integer      :: len
    integer, allocatable, dimension(:)       :: list
END TYPE ob_type

TYPE ob_rhs_type
    integer      :: len
    real(kind=WP), allocatable, dimension(:) :: list
END TYPE ob_rhs_type

type(ob_type), allocatable                    ::  ob_info(:)
type(ob_rhs_type), allocatable                ::  ob_2rhs(:)
!
! The fct part
integer                                       :: fct_iter=1
real(kind=WP),allocatable,dimension(:,:)      :: fct_LO, fct_HO           ! Low-order solution
real(kind=WP),allocatable,dimension(:,:)      :: fct_adf_h, fct_adf_h2    ! Antidif. horiz. contrib. from edges / backup for iterafive fct scheme
real(kind=WP),allocatable,dimension(:,:)      :: fct_adf_v, fct_adf_v2    ! Antidif. vert. fluxes from nodes    / backup for iterafive fct scheme

real(kind=WP),allocatable,dimension(:,:)      :: fct_ttf_max,fct_ttf_min
real(kind=WP),allocatable,dimension(:,:)      :: fct_plus,fct_minus
! Quadratic reconstruction part
integer,allocatable,dimension(:)              :: nlevels_nod2D_min, nn_num, nboundary_lay
real(kind=WP),allocatable,dimension(:,:,:)    :: quad_int_mat, quad_int_coef
integer,allocatable,dimension(:,:)            :: nn_pos
! MUSCL type reconstruction
integer,allocatable,dimension(:,:)            :: edge_up_dn_tri
real(kind=WP),allocatable,dimension(:,:,:)    :: edge_up_dn_grad

#if defined (__oasis)
  real(kind=8), allocatable, dimension(:)      :: lump2d_south, lump2d_north  
  integer, allocatable, dimension(:)           :: ind_south, ind_north    
#endif  
end module o_MESH
!==========================================================

!==========================================================
MODULE o_ARRAYS
USE o_PARAM
IMPLICIT NONE
! Arrays are described in subroutine array_setup  
real(kind=WP), allocatable    :: UV(:,:,:)
real(kind=WP), allocatable    :: UV_rhs(:,:,:), UV_rhsAB(:,:,:)
real(kind=WP), allocatable    :: eta_n(:), d_eta(:)
real(kind=WP), allocatable    :: ssh_rhs(:), Wvel(:,:), hpressure(:,:)
real(kind=WP), allocatable    :: Wvel_e(:,:), Wvel_i(:,:)
real(kind=WP), allocatable    :: CFL_z(:,:)
real(kind=WP), allocatable    :: stress_surf(:,:)
real(kind=WP), allocatable    :: T_rhs(:,:) 
real(kind=WP), allocatable    :: heat_flux(:), Tsurf(:) 
real(kind=WP), allocatable    :: heat_flux_old(:), Tsurf_old(:)  !PS
real(kind=WP), allocatable    :: S_rhs(:,:)
real(kind=WP), allocatable    :: tr_arr(:,:,:),tr_arr_old(:,:,:)
real(kind=WP), allocatable    :: del_ttf(:,:)

real(kind=WP), allocatable    :: water_flux(:), Ssurf(:)
real(kind=WP), allocatable    :: virtual_salt(:), relax_salt(:)
real(kind=WP), allocatable    :: water_flux_old(:), Ssurf_old(:) !PS
real(kind=WP), allocatable    :: Tclim(:,:), Sclim(:,:)
real(kind=WP), allocatable    :: Visc(:,:)
real(kind=WP), allocatable    :: Tsurf_t(:,:), Ssurf_t(:,:)
real(kind=WP), allocatable    :: tau_x_t(:,:), tau_y_t(:,:) 
real(kind=WP), allocatable    :: heat_flux_t(:,:), heat_rel_t(:,:), heat_rel(:) 
real(kind=WP), allocatable    :: coriolis(:), coriolis_node(:)
real(kind=WP), allocatable    :: relax2clim(:)
real(kind=WP), allocatable    :: MLD1(:), MLD2(:)
integer,       allocatable    :: MLD1_ind(:), MLD2_ind(:)

! Passive and age tracers
 character(4), allocatable    :: prog_tracer_name(:)
real(kind=WP), allocatable    :: tracer(:,:,:), tracer_rhs(:,:,:)   
!Tracer gradients&RHS      
real(kind=8), allocatable :: ttrhs(:,:)
real(kind=8), allocatable :: tr_xy(:,:,:)
real(kind=8), allocatable :: tr_z(:,:)

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
!real(kind=WP), allocatable,dimension(:,:,:) :: tsh, tsv, tsh_nodes
!real(kind=WP), allocatable,dimension(:,:)   :: hd_flux,vd_flux
!Isoneutral diffusivities (or xy diffusivities if Redi=.false)
real(kind=WP), allocatable :: Ki(:)

!_______________________________________________________________________________
! Arrays added for ALE implementation:
! --> layer thinkness at node and depthlayer for t=n and t=n+1
real(kind=WP), allocatable,dimension(:,:)   :: hnode, hnode_new, zbar_3d_n, Z_3d_n

! --> layer thinkness at elements, interpolated from hnode
real(kind=WP), allocatable,dimension(:,:)   :: helem

! --> thinkness of bottom elem (important for partial cells)
real(kind=WP), allocatable,dimension(:)     :: bottom_elem_thickness 

! --> The increment of total fluid depth on elements. It is used to update the matrix
real(kind=WP), allocatable,dimension(:)     :: dhe

! --> hbar, hbar_old: correspond to the elevation, but on semi-integer time steps.
real(kind=WP), allocatable,dimension(:)     :: hbar, hbar_old

! --> auxiliary array to store an intermediate part of the rhs computations.
real(kind=WP), allocatable,dimension(:)     :: ssh_rhs_old !, ssh_rhs_old2 !PS

! --> auxiliary array to store depth of layers and depth of mid level due to changing 
!     layer thinkness at every node
real(kind=WP), allocatable,dimension(:)     :: zbar_n, Z_n

! --> multiplication factor for surface boundary condition in 
!     diff_ver_part_impl_ale(tr_num) between linfs -->=0.0 and noninfs 
!     (zlevel,zstar...) --> = 1.0
real(kind=WP)                               :: is_nonlinfs

!_______________________________________________________________________________
!Monin-Obukhov correction
real*8,allocatable :: mo(:,:),mixlength(:)
!GM_stuff
real*8,allocatable :: bvfreq(:,:),mixlay_dep(:),bv_ref(:)

real(kind=WP), allocatable    :: fer_UV(:,:,:), fer_wvel(:,:)
real(kind=WP), target, allocatable    :: fer_c(:), fer_K(:), fer_gamma(:,:,:)
END MODULE o_ARRAYS
!==========================================================
