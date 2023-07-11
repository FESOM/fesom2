
module cvmix_tke
!! This module contains the main computations of diffusivities based on
!! TKE (following Gaspar'90)  with the calculation of the mixing length following (Blanke, B., P. Delecluse)
!!
!! @see  Gaspar, P., Y. Grégoris, and J.-M. Lefevre
!!       J. Geophys. Res., 95(C9), 16179–16193, doi:10.1029/JC095iC09p16179.
!!
!! @see  Blanke, B., P. Delecluse
!!       J. Phys. Oceanogr., 23, 1363–1388. doi: http://dx.doi.org/10.1175/1520-0485(1993)023<1363:VOTTAO>2.0.CO;2
!!
!! @author Nils Brueggemann, University of Hamburg / MPIMET
!! @author Hannah Kleppin, MPIMET/University of Hamburg
!! @author Oliver Gutjahr, MPIMET
!!
!!

use cvmix_kinds_and_types,    only : cvmix_r8,                     &
                                      CVMIX_OVERWRITE_OLD_VAL,     &
                                      CVMIX_SUM_OLD_AND_NEW_VALS,  &
                                      CVMIX_MAX_OLD_AND_NEW_VALS,  &
                                      cvmix_data_type,             &
                                      cvmix_global_params_type

use cvmix_utils,              only : cvmix_update_tke, solve_tridiag


implicit none
private 
save


!public member functions

public :: init_tke
public :: cvmix_coeffs_tke
public :: put_tke


!=================================================================================
!---------------------------------------------------------------------------------
! Interface to call the TKE parameterization! 
!---------------------------------------------------------------------------------

interface cvmix_coeffs_tke
  module procedure integrate_tke              ! calculation if prognostic TKE equation
  module procedure tke_wrap                   ! necessary to handle old/new values and to hand over user_defined constants
end interface cvmix_coeffs_tke 

!---------------------------------------------------------------------------------
! Interface to put values to TKE variables
!---------------------------------------------------------------------------------

interface put_tke
  module procedure cvmix_tke_put_tke_int
  module procedure cvmix_tke_put_tke_real
  module procedure cvmix_tke_put_tke_logical
end interface put_tke

!=================================================================================

! types for TKE
type, public :: tke_type
private

real(cvmix_r8)       ::  &
  c_k                   ,& ! 
  c_eps                 ,& ! dissipation parameter
  cd                    ,& ! 
  alpha_tke             ,& ! 
  mxl_min               ,& ! minimum value for mixing length
  kappaM_min            ,& ! minimum value for Kappa momentum
  kappaM_max            ,& ! maximum value for Kappa momentum
  tke_surf_min          ,& ! minimum value for surface TKE 
  tke_min                  ! minimum value for TKE, necessary to set this value when 
                           ! run without IDEMIX, since there are no sources for TKE in the deep ocean otherwise

integer               :: &
  tke_mxl_choice        ,& ! choice of calculation of mixing length;
  handle_old_vals          ! Flag for what to do with old values of Vmix_vars%TKE,tke_diss,KappaM_iface,KappaH_iface

logical                :: &
  only_tke               ,&
  use_ubound_dirichlet   ,&
  use_lbound_dirichlet                           

end type tke_type

type(tke_type), target :: tke_constants_saved 

 contains
 
!=================================================================================

subroutine init_tke(c_k, c_eps, cd, alpha_tke, mxl_min, KappaM_min, KappaM_max, &
                    tke_mxl_choice, use_ubound_dirichlet, use_lbound_dirichlet, &
                    handle_old_vals, only_tke, tke_min, tke_surf_min, &
                    tke_userdef_constants)

! This subroutine sets user or default values for TKE parameters

real(cvmix_r8),optional, intent(in)           ::  & 
  c_k                                            ,&
  c_eps                                          ,&
  cd                                             ,&
  alpha_tke                                      ,&
  mxl_min                                        ,&
  KappaM_min                                     ,& 
  KappaM_max                                     ,&
  tke_surf_min                                   ,&
  tke_min

integer, intent(in),optional                   :: &
  tke_mxl_choice                                 ,& 
  handle_old_vals

logical, intent(in), optional                  :: &
  only_tke                                       ,&
  use_ubound_dirichlet                           ,&
  use_lbound_dirichlet                           

type(tke_type), intent(inout),target, optional :: &
  tke_userdef_constants

! FIXME: not sure about the allowed ranges for TKE parameters
if (present(c_k)) then
  if(c_k.lt. 0.05d0 .or. c_k .gt. 0.3d0) then
    print*, "ERROR:c_k can only be allowed_range"
    stop 1
  end if
  call put_tke('c_k', c_k, tke_userdef_constants)
else
  call put_tke('c_k',0.1_cvmix_r8 , tke_userdef_constants)
end if

if (present(c_eps)) then
  if(c_eps.lt. 0.5d0 .or. c_eps .gt. 1.d0) then
    print*, "ERROR:c_eps can only be allowed_range"
    stop 1
  end if
  call put_tke('c_eps', c_eps, tke_userdef_constants)
else
  call put_tke('c_eps', 0.7_cvmix_r8, tke_userdef_constants)
end if

if (present(cd)) then
  if(cd.lt. 0.1d0 .or. cd .gt. 30.d0) then
    print*, "ERROR:cd can only be allowed_range"
    stop 1
  end if
  call put_tke('cd', cd, tke_userdef_constants)
else
  call put_tke('cd', 3.75_cvmix_r8, tke_userdef_constants)
end if

if (present(alpha_tke)) then
  if(alpha_tke.lt. 1.d0 .or. alpha_tke .gt. 30.d0) then
    print*, "ERROR:alpha_tke can only be allowed_range"
    stop 1
  end if
  call put_tke('alpha_tke', alpha_tke, tke_userdef_constants)
else
  call put_tke('alpha_tke', 30._cvmix_r8, tke_userdef_constants)
end if

if (present(mxl_min)) then
  if(mxl_min.lt. 1.d-12 .or. mxl_min .gt. 0.4d0) then
    print*, "ERROR:mxl_min can only be allowed_range"
    stop 1
  end if
  call put_tke('mxl_min', mxl_min, tke_userdef_constants)
else
  call put_tke('mxl_min', 1._cvmix_r8-8, tke_userdef_constants)
end if

if (present(KappaM_min)) then
  if(KappaM_min.lt. 0.d0 .or. KappaM_min .gt. 1.d0) then
    print*, "ERROR:KappaM_min can only be allowed_range"
    stop 1
  end if
  call put_tke('kappaM_min', KappaM_min, tke_userdef_constants)
else
  call put_tke('kappaM_min', 0._cvmix_r8, tke_userdef_constants)
end if

if (present(KappaM_max)) then
  if(KappaM_max.lt. 1.d0 .or. KappaM_max .gt. 100.d0) then
    print*, "ERROR:kappaM_max can only be allowed_range"
    stop 1
  end if
  call put_tke('kappaM_max', KappaM_max, tke_userdef_constants)
else
  call put_tke('kappaM_max', 100._cvmix_r8, tke_userdef_constants)
end if

if (present(tke_mxl_choice)) then
  if(tke_mxl_choice.lt. 1 .or. tke_mxl_choice .gt. 2 ) then
    print*, "ERROR:tke_mxl_choice can only be 1 or 2"
    stop 1
  end if
  call put_tke('tke_mxl_choice', tke_mxl_choice, tke_userdef_constants)
else
  call put_tke('tke_mxl_choice', 2, tke_userdef_constants)
end if

if (present(handle_old_vals)) then
  if(handle_old_vals.lt. 1 .or. handle_old_vals.gt. 3 ) then
    print*, "ERROR:handle_old_vals can only be 1 to 3"
    stop 1
  end if
  call put_tke('handle_old_vals', handle_old_vals, tke_userdef_constants)
else
  call put_tke('handle_old_vals', 1, tke_userdef_constants)
end if

if (present(tke_min)) then
  if(tke_min.lt. 1.d-7 .or. tke_min.gt. 1.d-4 ) then
    print*, "ERROR:tke_min can only be 10^-7 to 10^-4"
    stop 1
  end if
  call put_tke('tke_min', tke_min, tke_userdef_constants)
else
  call put_tke('tke_min', 1._cvmix_r8-6, tke_userdef_constants)
end if

if (present(tke_surf_min)) then
  if(tke_surf_min.lt. 1.d-7 .or. tke_surf_min.gt. 1.d-2 ) then
    print*, "ERROR:tke_surf_min can only be 10^-7 to 10^-4"
    stop 1
  end if
  call put_tke('tke_surf_min', tke_surf_min, tke_userdef_constants)
else
  call put_tke('tke_surf_min', 1._cvmix_r8-4, tke_userdef_constants)
end if

if (present(use_ubound_dirichlet)) then
  call put_tke('use_ubound_dirichlet', use_ubound_dirichlet, tke_userdef_constants)
else
  call put_tke('use_ubound_dirichlet', .false., tke_userdef_constants)
end if

if (present(use_lbound_dirichlet)) then
  call put_tke('use_lbound_dirichlet', use_lbound_dirichlet, tke_userdef_constants)
else
  call put_tke('use_lbound_dirichlet', .false., tke_userdef_constants)
end if

if (present(only_tke)) then

  call put_tke('only_tke', only_tke, tke_userdef_constants)
else
  call put_tke('only_tke', .true., tke_userdef_constants)
end if

end subroutine init_tke

!=================================================================================

subroutine tke_wrap(Vmix_vars, Vmix_params, tke_userdef_constants)

! This subroutine is necessary to handle old/new values and to hand over the TKE parameters set in previous subroutine
! This subroutine should be called from calling ocean model or driver

type(tke_type), intent(in), optional, target ::  &
  tke_userdef_constants                            !

type(cvmix_data_type), intent(inout)         ::  & 
  Vmix_vars                                        ! 

type(cvmix_global_params_type), intent(in)   ::  &
  Vmix_params                                      !

real(cvmix_r8), dimension(Vmix_vars%nlev)  ::  &
  cvmix_int_1                                       ,& !
  cvmix_int_2                                       ,& !
  cvmix_int_3                                       ,& !
  tke_Tbpr                                          ,&
  tke_Tspr                                          ,&
  tke_Tdif                                          ,&
  tke_Tdis                                          ,&
  tke_Twin                                          ,&
  tke_Tiwf                                          ,&
  tke_Tbck                                          ,&
  tke_Ttot                                          ,&
  !tke                                               ,&
  tke_Lmix                                          ,&
  tke_Pr                                            ,&
  new_KappaM                                    ,& !
  new_KappaH                                    ,& !
  new_tke                                       ,& !
  new_tke_diss                          

integer                                       :: &
  nlev                                          ,& !
  max_nlev                                         !
! nils
integer :: i, j, tstep_count

type(tke_type), pointer                       :: &
  tke_constants_in                                 !

! FIXME: nils: this routine is never called. 
! This is an emergency stop in case th routine is called.
write(*,*) 'i am wrapping'
stop

tke_constants_in => tke_constants_saved

if (present(tke_userdef_constants)) then
  tke_constants_in => tke_userdef_constants
end if

nlev = Vmix_vars%nlev
max_nlev = Vmix_vars%max_nlev

! call to actual computation of TKE parameterization
call cvmix_coeffs_tke( &
                 ! parameter
                 dzw          = Vmix_vars%dzw,                                    &
                 dzt          = Vmix_vars%dzt,                                    &
                 nlev         = nlev,                                             &
                 max_nlev     = max_nlev,                                         &
                 dtime        = Vmix_vars%dtime,                                  & 
                 rho_ref      = Vmix_vars%rho_ref,                                & 
                 grav         = Vmix_params%Gravity,                              & 
                 ! essentials
                 tke_old      = Vmix_vars%tke,            & ! in
                 tke_new      = new_tke,                  & ! out
                 KappaM_out   = new_KappaM,               & ! out
                 KappaH_out   = new_KappaH,               & ! out
                 Ssqr         = Vmix_vars%Ssqr_iface,                             &
                 Nsqr         = Vmix_vars%Nsqr_iface,                             &
                 old_kappaM   = Vmix_vars%KappaM_iface,                           &
                 old_KappaH   = Vmix_vars%KappaH_iface,                           &
                 ! FIXME: nils: better calc IDEMIX Ri directly in ! CVMIX/IDEMIX
                 alpha_c      = Vmix_vars%alpha_c,                                &
                 E_iw         = Vmix_vars%E_iw,                                   &
                 ! forcing
                 forc_tke_surf= Vmix_vars%forc_tke_surf,                          &
                 forc_rho_surf= Vmix_vars%forc_rho_surf,                          &
                 bottom_fric  = Vmix_vars%bottom_fric,                            &
                 iw_diss      = Vmix_vars%iw_diss,                                &
                 ! diagnostics
                 tke_Tbpr     = tke_Tbpr,                                         &
                 tke_Tspr     = tke_Tspr,                                         &
                 tke_Tdif     = tke_Tdif,                                         &
                 tke_Tdis     = tke_Tdis,                                         &
                 tke_Twin     = tke_Twin,                                         &
                 tke_Tiwf     = tke_Tiwf,                                         &
                 tke_Tbck     = tke_Tbck,                                         &
                 tke_Ttot     = tke_Ttot,                                         &
                 tke_Lmix     = tke_Lmix,                                         &
                 tke_Pr       = tke_Pr,                                           &
                 ! debugging
                 cvmix_int_1   = cvmix_int_1,             &
                 cvmix_int_2   = cvmix_int_2,             &
                 cvmix_int_3   = cvmix_int_3,             &
                 i = i , &
                 j = j , &
                 tstep_count = tstep_count , &
                 !tke_diss_out = new_tke_diss,                                &
!                 nlev         = Vmix_vars%nlev,                                   &
!                 max_nlev     = Vmix_vars%max_nlev,                               &
                 !old_tke_diss = Vmix_vars%tke_diss,                               &
                 !tke          = tke,                                              &
                 !Kappa_GM     = Vmix_vars%Kappa_GM,                               &
                 tke_userdef_constants = tke_userdef_constants)


!!update Vmix_vars to new values
!call cvmix_update_tke(tke_constants_in%handle_old_vals,                           &
!!                     Vmix_vars%nlev,                                              &
!                     nlev,                                                        &  
!                     !tke_diss_out = Vmix_vars%tke_diss,                           &
!                     new_tke_diss = new_tke_diss,                                 &
!                     tke_new      = Vmix_vars%tke,                                &
!                     new_tke      = new_tke,                                      &
!                     KappaM_out   = Vmix_vars%KappaM_iface,                       &
!                     new_KappaM   = new_KappaM,                                   &
!                     KappaH_out   = Vmix_vars%KappaH_iface,                       &
!                     new_KappaH   = new_KappaH)

end subroutine tke_wrap

!=================================================================================

subroutine integrate_tke( &
                         i,                    & ! FIXME: nils: for debuging
                         j,                    & ! FIXME: nils: for debuging
                         tstep_count ,         & ! FIXME: nils: for debuging
                         !tke_diss_out,         & ! FIXME: nils: today: delete?
                         tke_old,              & ! FIXME: nils: today: rename?
                         tke_new,              & ! FIXME: nils: today: rename?
                         KappaM_out,           &
                         KappaH_out,           &
                         cvmix_int_1,          & ! FIXME: nils: for debuging
                         cvmix_int_2,          & ! FIXME: nils: for debuging
                         cvmix_int_3,          & ! FIXME: nils: for debuging
                         dzw,                  &
                         dzt,                  &
                         nlev,                 &
                         max_nlev,             &
                         !old_tke_diss,         & ! FIXME: nils: today: delete?
                         Ssqr,                 &
                         Nsqr,                 & 
                         tke_Tbpr,             & ! diagnostic
                         tke_Tspr,             & ! diagnostic
                         tke_Tdif,             & ! diagnostic
                         tke_Tdis,             & ! diagnostic
                         tke_Twin,             & ! diagnostic
                         tke_Tiwf,             & ! diagnostic
                         tke_Tbck,             & ! diagnostic
                         tke_Ttot,             & ! diagnostic
                         !tke,                  &
                         tke_Lmix,             & ! diagnostic
                         tke_Pr,               & ! diagnostic
                         forc_tke_surf,        &
                         E_iw,                 &
                         dtime,                &
                         bottom_fric,          &
                         old_KappaM,           &
                         old_KappaH,           &
                         iw_diss,              & ! FIXME: nils: rename?
                         forc_rho_surf,        &
                         !Kappa_GM,             & ! FIXME: nils: today: delete?
                         rho_ref,              & ! FIXME: today: put to initialize
                         grav,                 & ! FIXME: today: put to initialize
                         alpha_c,              & ! FIXME: today: put to initialize
                         tke_userdef_constants)
  
  type(tke_type), intent(in), optional, target                 :: &
    tke_userdef_constants
  
  integer,intent(in)                                           :: &
    nlev                                                         ,& !
    max_nlev                                                      !
  integer,intent(in)                                           :: &
    i, j, tstep_count
  
  ! OLD values 
  real(cvmix_r8), dimension(nlev+1), intent(in)                :: & 
    tke_old                                                      ,& !
    !old_tke_diss                                                 ,& !
    old_KappaM                                                   ,& !
    old_KappaH                                                   ,& !
    dzt                                                             !
  real(cvmix_r8), dimension(nlev+1), intent(in)             :: & 
    Ssqr                                                         ,& !
    Nsqr                                                         
   
  real(cvmix_r8), dimension(nlev), intent(in)                  :: &
    dzw                                                             !
   
  ! IDEMIX variables, if run coupled iw_diss is added as forcing to TKE
  real(cvmix_r8), dimension(max_nlev), intent(in), optional  :: &
    E_iw                                                         ,& !  
    alpha_c                                                      ,& !
    iw_diss                                                         !
  
  real(cvmix_r8), intent(in)                                   :: & 
    bottom_fric                                                  ,& !
    forc_rho_surf                                                ,& ! FIXME: nils: what is this?
    rho_ref                                                      ,& !
    dtime                                                        ,& ! time step
    grav                                                            ! gravity constant
  
  ! FIXME: nils: why inout? why not only in?
  real(cvmix_r8), intent(in)                                   :: & 
    forc_tke_surf
              
  !real(cvmix_r8),dimension(nlev+1), intent(in), optional       :: &
  !  Kappa_GM                                                        !
  
  ! NEW values
  real(cvmix_r8), dimension(nlev+1), intent(out)             :: &
    tke_new                                                      ,& ! 
    !tke_diss_out                                                 ,& !
    KappaM_out                                                   ,& !
    KappaH_out
  
  ! diagnostics
  real(cvmix_r8), dimension(nlev+1) ::                &
     tke_Tbpr                                                     ,&
     tke_Tspr                                                     ,&
     tke_Tdif                                                     ,&
     tke_Tdis                                                     ,&
     tke_Twin                                                     ,&
     tke_Tiwf                                                     ,&
     tke_Tbck                                                     ,&
     tke_Ttot                                                     ,&
     !tke                                                          ,&
     tke_Lmix                                                     ,&
     tke_Pr                                                       !,&
  real(cvmix_r8), dimension(nlev+1), intent(out) ::                &
     cvmix_int_1                                                  ,&
     cvmix_int_2                                                  ,&
     cvmix_int_3                                                
  
  ! local variables
  real(cvmix_r8), dimension(nlev+1)                            :: &
    tke_unrest                                                   ,& ! copy of tke before restorring to background value
    tke_upd                                                      ,& ! copy of tke before in which surface/bottom values given by Dirichlet boundary conditions
    mxl                                                          ,& ! mixing length scale (m)
    sqrttke                                                      ,& ! square root of TKE (m/s)
    prandtl                                                      ,& ! Prandtl number
    Rinum                                                        ,& ! Richardson number 
    K_diss_v                                                     ,& ! shear production of TKE (m^2/s^3)
    P_diss_v                                                     ,& ! buoyancy production of TKE (m^2/s^3)
    forc                                                            ! combined forcing for TKE (m^2/s^3)
  
  real(cvmix_r8) :: tke_surf, tke_bott
  
  real(cvmix_r8)                                               :: &
    ! FIXME: nils: today: delete tke_surf_corr?
    !tke_surf_corr                                                ,& ! correction of surface density according to surface buoyancy flux
    alpha_tke                                                    ,& ! {30}
    c_eps                                                        ,& ! {0.7}
    cd                                                           ,& ! (3.75}
    KappaM_max                                                   ,& ! 
    mxl_min                                                      ,& ! {1e-8}
    c_k                                                          ,& ! {0.1}
    tke_surf_min                                                 ,& ! {1e-4}
    tke_min                                                         ! {1e-6}
  integer :: tke_mxl_choice

  logical :: only_tke, use_ubound_dirichlet, use_lbound_dirichlet
  
  real(cvmix_r8)                                               :: &
    zzw                                                          ,& ! depth of interface k 
    depth                                                        ,&! total water depth
    diff_surf_forc                                               ,&
    diff_bott_forc
    
  ! input to tri-diagonal solver
  real(cvmix_r8), dimension(nlev+1)                            :: &
    a_dif                                                        ,& !
    b_dif                                                        ,& !
    c_dif                                                        ,& !
    a_tri                                                        ,& !
    b_tri                                                        ,& !
    c_tri                                                        ,& !
    d_tri                                                        ,& !
    ke                                                              !  diffusivity for tke
  
  integer :: k, kk, kp1
  
  ! FIXME: nils: Where are these values set?
  real(cvmix_r8) :: kappaM_min
  integer :: recnum
  
  type(tke_type), pointer :: tke_constants_in
  
  ! FIXME: nils: What happens here?
  tke_constants_in => tke_constants_saved
  if (present(tke_userdef_constants)) then
   tke_constants_in => tke_userdef_constants
  end if

  ! FIXME: nils: What should we do with height of last grid box dzt(nlev+1)?
  !              This should not be as thick as the distance to the next tracer
  !              point (which is a dry point). 
  !              Be careful if you divide by 0.5 here. Maybe later we use ddpo
  !              or something like this that might include already this factor.

  ! initialize diagnostics
  tke_Tbpr = 0.0
  tke_Tspr = 0.0
  tke_Tdif = 0.0
  tke_Tdis = 0.0
  tke_Twin = 0.0
  tke_Tiwf = 0.0
  tke_Tbck = 0.0
  tke_Ttot = 0.0
 
  !---------------------------------------------------------------------------------
  ! set tke_constants locally
  !---------------------------------------------------------------------------------
 
  alpha_tke  = tke_constants_in%alpha_tke
  c_eps      = tke_constants_in%c_eps
  cd         = tke_constants_in%cd
  KappaM_max = tke_constants_in%KappaM_max
  mxl_min    = tke_constants_in%mxl_min
  c_k        = tke_constants_in%c_k
  tke_min    = tke_constants_in%tke_min
  tke_surf_min   = tke_constants_in%tke_surf_min
  tke_mxl_choice = tke_constants_in%tke_mxl_choice
  only_tke = tke_constants_in%only_tke
  use_ubound_dirichlet = tke_constants_in%use_ubound_dirichlet
  use_lbound_dirichlet = tke_constants_in%use_lbound_dirichlet

  ! FIXME: nils: Is kappaM_min ever used?
  ! FIXME: use kappaM_min from namelist
  ! FIXME: where is kappaM_min used?
  ! FIXME: start with small k in all KappaM_min 
  kappaM_min = 0.0
  !KappaM_out=max(tke_userdef_constants%kappaM_min,KappaM_out)
  
  !---------------------------------------------------------------------------------
  ! Part 1: calculate mixing length scale
  !---------------------------------------------------------------------------------
  sqrttke = sqrt(max(0d0,tke_old))
 
  ! turbulent mixing length
  mxl = sqrt(2D0)*sqrttke/sqrt(max(1d-12,Nsqr))

  ! constrain mixing length scale as in MITgcm
  if (tke_mxl_choice==2) then 
    !FIXME: What should we do at the surface and bottom?
    mxl(1) = 0.d0
    mxl(nlev+1) = 0.d0
    do k=2,nlev
      mxl(k) = min(mxl(k), mxl(k-1)+dzw(k-1))
    enddo
    mxl(nlev) = min(mxl(nlev), mxl_min+dzw(nlev))
    do k=nlev-1,2,-1
      mxl(k) = min(mxl(k), mxl(k+1)+dzw(k))
    enddo
    mxl= max(mxl,mxl_min)
  ! bounded by the distance to surface/bottom
  elseif (tke_mxl_choice==3) then
    depth = sum(dzw(1:nlev))
    do k=2,nlev+1
     zzw = sum(dzw(1:k-1))
     mxl(k) = min(zzw,mxl(k),depth-zzw)
    enddo
    mxl(1) = mxl(2)
    mxl= max(mxl,mxl_min)
  else
    write(*,*) 'Wrong choice of tke_mxl_choice. Aborting...'
    stop
  endif
     
  !---------------------------------------------------------------------------------
  ! Part 2: calculate diffusivities
  !---------------------------------------------------------------------------------
  ! see. Blanke and Delecluse 1993, eq. 2.25
  KappaM_out = min(KappaM_max,c_k*mxl*sqrttke)
  
  ! Richardson number --> see. Blanke and Delecluse 1993, eq. 2.18
  Rinum = Nsqr/max(Ssqr,1d-12)
 
  ! FIXME: nils: Check this later if IDEMIX is coupled.
  ! FIXME: nils: Why E_iw**2 and not dissipation with mixed time level?
  ! FIXME: nils: Why not passing Rinum as Rinum_idemix to tke scheme?
  if (.not.only_tke) then  !IDEMIX is on
    Rinum = min(Rinum,KappaM_out*Nsqr/max(1d-12,alpha_c*E_iw**2))
  end if
  
  ! Richardson number dependent expression is used for the P_rt (Osborn 1980, 
  ! Crawford 1982) --> see. Blanke and Delecluse 1993, eq. 2.33 & 2.34
  ! Why 6.6*... ???
  prandtl=max(1d0,min(10d0,6.6*Rinum))
  
  ! see. Blanke and Delecluse 1993, eq. 2.26
  KappaH_out=KappaM_out/prandtl

  !---------------------------------------------------------------------------------
  ! Part 3: tke forcing
  !---------------------------------------------------------------------------------
  ! initialize forcing
  forc = 0.0

  ! --- forcing by shear and buoycancy production
  K_diss_v   = Ssqr*KappaM_out 
  P_diss_v   = Nsqr*KappaH_out
  ! FIXME: nils: Is forc_rho_surf set somewhere?
  ! FIXME: nils: What does forc_rho_surf mean?
  P_diss_v(1) = -forc_rho_surf*grav/rho_ref
  forc = forc + K_diss_v - P_diss_v
   
  ! --- forcing by internal wave dissipation
  if (.not.only_tke) then
    forc = forc + iw_diss
  endif

  !---------------------------------------------------------------------------------
  ! Part 4: vertical diffusion and dissipation is solved implicitely 
  !---------------------------------------------------------------------------------
  ke = 0.d0
  do k = 1, nlev
    kp1 = min(k+1,nlev)
    kk  = max(k,2) 
    ke(k) = alpha_tke*0.5*(KappaM_out(kp1)+KappaM_out(kk))
  enddo

  !--- c is lower diagonal of matrix
  do k=1,nlev
    !c_dif(k) = delta(k+1)/dzt(k)
    c_dif(k) = ke(k)/( dzt(k)*dzw(k) )
  enddo
  c_dif(nlev+1) = 0.d0 ! not part of the diffusion matrix, thus value is arbitrary
 
  !--- b is main diagonal of matrix
  do k=2,nlev
    !b_dif(k) = delta(k)/dzt(k)+delta(k+1)/dzt(k)
    b_dif(k) = ke(k-1)/( dzt(k)*dzw(k-1) ) + ke(k)/( dzt(k)*dzw(k) )
  enddo
 
  !--- a is upper diagonal of matrix
  do k=2,nlev+1
    !a_dif(k) = delta(k)/dzt(k)
    a_dif(k) = ke(k-1)/( dzt(k)*dzw(k-1) )
  enddo
  a_dif(1) = 0.d0 ! not part of the diffusion matrix, thus value is arbitrary

  ! copy tke_old
  tke_upd = tke_old
 
  ! upper boundary condition
  if (use_ubound_dirichlet) then
    sqrttke(1)      = 0.d0 ! to suppres dissipation for k=1
    forc(1)         = 0.d0 ! to suppres forcing for k=1
    tke_surf        = max(tke_surf_min, cd*forc_tke_surf)
    tke_upd(1)      = tke_surf
    ! add diffusive part that depends on tke_surf to forcing
    diff_surf_forc  = a_dif(2)*tke_surf
    forc(2)         = forc(2)+diff_surf_forc
    a_dif(2)        = 0.d0 ! and set matrix element to zero
    b_dif(1)        = 0.d0 ! 0 line in matrix for k=1
    c_dif(1)        = 0.d0 ! 0 line in matrix for k=1
  else
    ! add wind forcing
    forc(1) = forc(1) + (cd*forc_tke_surf**(3./2.))/(dzt(1))
    b_dif(1)        = ke(1)/( dzt(1)*dzw(1) )
    diff_surf_forc  = 0.0
  endif

  ! lower boundary condition
  if (use_lbound_dirichlet) then
    sqrttke(nlev+1) = 0.d0 ! to suppres dissipation for k=nlev+1
    forc(nlev+1)    = 0.d0 ! to suppres forcing for k=nlev+1
    ! FIXME: make tke_bott dependend on bottom friction?
    tke_bott        = tke_min
    !tke_old(nlev+1) = tke_bott
    tke_upd(nlev+1) = tke_bott
    ! add diffusive part that depends on tke_bott to forcing
    diff_bott_forc  = c_dif(nlev)*tke_bott
    forc(nlev)      = forc(nlev)+diff_bott_forc
    c_dif(nlev)     = 0.d0 ! and set matrix element to zero
    b_dif(nlev+1)   = 0.d0 ! 0 line in matrix for k=nlev+1
    a_dif(nlev+1)   = 0.d0 ! 0 line in matrix for k=nlev+1   
  else
    b_dif(nlev+1)   = ke(nlev)/( dzt(nlev+1)*dzw(nlev) )
    diff_bott_forc  = 0.0
  endif

  !--- construct tridiagonal matrix to solve diffusion and dissipation implicitely
  a_tri = -dtime*a_dif
  b_tri = 1+dtime*b_dif
  b_tri(2:nlev) = b_tri(2:nlev) + dtime*c_eps*sqrttke(2:nlev)/mxl(2:nlev)
  c_tri = -dtime*c_dif
 
  !--- d is r.h.s. of implicite equation (d: new tke with only explicite tendencies included)
  d_tri(1:nlev+1)  = tke_upd(1:nlev+1) + dtime*forc(1:nlev+1)
 
  ! solve the tri-diag matrix
  call solve_tridiag(a_tri, b_tri, c_tri, d_tri, tke_new, nlev+1)
 
  ! --- diagnose implicite tendencies (only for diagnostics)
  ! vertical diffusion of TKE
  do k=2,nlev
    tke_Tdif(k) = a_dif(k)*tke_new(k-1) - b_dif(k)*tke_new(k) + c_dif(k)*tke_new(k+1)
  enddo
  tke_Tdif(1) = - b_dif(1)*tke_new(1) + c_dif(1)*tke_new(2)
  tke_Tdif(nlev+1) = a_dif(nlev+1)*tke_new(nlev) - b_dif(nlev+1)*tke_new(nlev+1)
  tke_Tdif(2) = tke_Tdif(2) + diff_surf_forc
  tke_Tdif(nlev) = tke_Tdif(nlev) + diff_bott_forc

  ! flux out of first box due to diffusion with Dirichlet boundary value of TKE
  ! (tke_surf=tke_upd(1)) and TKE of box below (tke_new(2))
  if (use_ubound_dirichlet) then
    tke_Tdif(1) = - ke(1)/dzw(1)/dzt(1) &
                    * (tke_surf-tke_new(2))
  endif 
  if (use_lbound_dirichlet) then
    k = nlev+1
    tke_Tdif(k) = ke(k-1)/dzw(k-1)/dzt(k) &
                    * (tke_new(k-1)-tke_bott)
  endif 
 
  ! dissipation of TKE
  tke_Tdis = 0.d0
  tke_Tdis(2:nlev) = -c_eps/mxl(2:nlev)*sqrttke(2:nlev)*tke_new(2:nlev)
  !tke_diss_out(1:nlev+1) = c_eps/mxl(1:nlev+1)*sqrttke(1:nlev+1)*tke_new(1:nlev)

  !---------------------------------------------------------------------------------
  ! Part 5: reset tke to bounding values
  !---------------------------------------------------------------------------------
  ! copy of unrestored tke to diagnose energy input by restoring
  tke_unrest = tke_new

  ! FIXME: nils: today: delete
  ! add TKE if surface density flux drains TKE in uppermost box
  !tke_surf_corr = 0.0
 
  ! restrict values of TKE to tke_min, if IDEMIX is not used
  if (only_tke) then
    tke_new = MAX(tke_new, tke_min)
  end if 
 
  !---------------------------------------------------------------------------------
  ! Part 6: Assign diagnostic variables
  !---------------------------------------------------------------------------------
  ! tke_Ttot =   tke_Tbpr + tke_Tspr + tke_Tdif + tke_Tdis 
  !            + tke_Twin + tke_Tiwf
  tke_Tbpr = -P_diss_v
  tke_Tspr = K_diss_v
  !tke_Tdif is set above
  !tke_Tdis = -tke_diss_out
  tke_Tbck = (tke_new-tke_unrest)/dtime
  if (use_ubound_dirichlet) then
    tke_Twin(1) = (tke_new(1)-tke_old(1))/dtime - tke_Tdif(1) 
    tke_Tbck(1) = 0.0
  else
    !tke_Twin(1) = forc_tke_surf/(dzt(1)) 
    tke_Twin(1) = (cd*forc_tke_surf**(3./2.))/(dzt(1))
  endif
  ! FIXME: Find better name for tke_Twin either tke_Tbou or use tke_Tsur
  ! tke_Tbot
  if (use_lbound_dirichlet) then
    tke_Twin(nlev+1) = (tke_new(nlev+1)-tke_old(nlev+1))/dtime - tke_Tdif(nlev+1) 
    tke_Tbck(nlev+1) = 0.0
  else
    !FIXME: no flux condition so far, add bottom friction later
    tke_Twin(nlev+1) = 0.0
  endif

  tke_Tiwf = iw_diss
  tke_Ttot = (tke_new-tke_old)/dtime
  !tke = tke_new
  tke_Lmix = mxl
  tke_Pr = prandtl
   
  ! -----------------------------------------------
  ! the rest is for debugging
  ! -----------------------------------------------
  cvmix_int_1 = KappaM_out
  cvmix_int_2 = 0.0
  cvmix_int_2(1) = tke_surf
  cvmix_int_3 = Nsqr
  !cvmix_int_1 = forc
  !cvmix_int_2 = Nsqr
  !cvmix_int_3 = Ssqr

  if (.false.) then
  !  write(*,*) 'i = ', i, 'j = ', j, 'tstep_count = ', tstep_count
  if (i==45 .and. j==10) then
  !if (i==45 .and. j==10 .and. tstep_count==10) then
  ! -----------------------------------------------
 
    write(*,*) '================================================================================'
    write(*,*) 'i = ', i, 'j = ', j, 'tstep_count = ', tstep_count
    write(*,*) 'nlev = ', nlev
    write(*,*) 'dtime = ', dtime
    write(*,*) 'dzt = ', dzt
    write(*,*) 'dzw = ', dzw
    write(*,*) 'Nsqr = ', Nsqr
    write(*,*) 'Ssqr = ', Ssqr
    !write(*,*) 'tho = ', tho(i,j,1:nlev)
    !write(*,*) 'sao = ', sao(i,j,1:nlev)
    !write(*,*) 'bottom_fric = ', bottom_fric
    !write(*,*) 'forc_tke_surf = ', forc_tke_surf
    write(*,*) 'sqrttke = ', sqrttke
    write(*,*) 'mxl = ', mxl
    write(*,*) 'KappaM_out = ', KappaM_out
    write(*,*) 'KappaH_out = ', KappaH_out
    write(*,*) 'forc = ', forc
    !write(*,*) 'Rinum = ', Rinum
    write(*,*) 'prandtl = ', prandtl
    !write(*,*) 'checkpoint d_tri'
    !write(*,*) 'K_diss_v = ', K_diss_v
    !write(*,*) 'P_diss_v = ', P_diss_v
    !write(*,*) 'delta = ', delta
    write(*,*) 'ke = ', ke
    write(*,*) 'a_tri = ', a_tri
    write(*,*) 'b_tri = ', b_tri
    write(*,*) 'c_tri = ', c_tri
    write(*,*) 'd_tri = ', d_tri
    !write(*,*) 'tke_old = ', tke_old
    write(*,*) 'tke_new = ', tke_new
    write(*,*) 'tke_Tbpr = ', tke_Tbpr
    write(*,*) 'tke_Tspr = ', tke_Tspr
    write(*,*) 'tke_Tdif = ', tke_Tdif
    write(*,*) 'tke_Tdis = ', tke_Tdis
    write(*,*) 'tke_Twin = ', tke_Twin
    write(*,*) 'tke_Tiwf = ', tke_Tiwf
    write(*,*) 'tke_Ttot = ', tke_Ttot
    write(*,*) 'tke_Ttot - tke_Tsum = ', &
      tke_Ttot-(tke_Tbpr+tke_Tspr+tke_Tdif+tke_Tdis+tke_Twin+tke_Tiwf)
    !write(*,*) 'dzw = ', dzw
    !write(*,*) 'dzt = ', dzt
    ! FIXME: partial bottom cells!!
    ! namelist parameters
    write(*,*) 'c_k = ', c_k
    write(*,*) 'c_eps = ', c_eps
    write(*,*) 'alpha_tke = ', alpha_tke
    write(*,*) 'mxl_min = ', mxl_min
    write(*,*) 'kappaM_min = ', kappaM_min
    write(*,*) 'kappaM_max = ', kappaM_max
    ! FIXME: Make tke_mxl_choice available!
    !write(*,*) 'tke_mxl_choice = ', tke_mxl_choice
    !write(*,*) 'cd = ', cd
    write(*,*) 'tke_min = ', tke_min
    write(*,*) 'tke_surf_min = ', tke_surf_min
    write(*,*) 'only_tke = ', only_tke
    write(*,*) 'use_ubound_dirichlet = ', use_ubound_dirichlet
    write(*,*) 'use_lbound_dirichlet = ', use_lbound_dirichlet
    !write(*,*) 'tke(nlev) = ', tke(nlev), 'tke(nlev+1) = ', tke(nlev+1)
    !write(*,*) 'tke(nlev+2) = ', tke(nlev+2)
    write(*,*) '================================================================================'
  !end if
  !if (i==45 .and. j==10 .and. tstep_count==10) then
    !stop
  end if ! if (i==, j==, tstep==)
  end if ! if (.true./.false.)
end subroutine integrate_tke

!=================================================================================

subroutine cvmix_tke_put_tke_int(varname,val,tke_userdef_constants)
!This subroutine puts integer values to TKE variables
!IN
    character(len=*),           intent(in) :: varname
    integer,                    intent(in) :: val
!OUT   
    type(tke_type), intent(inout), target, optional:: tke_userdef_constants
    type(tke_type), pointer :: tke_constants_out

  tke_constants_out=>tke_constants_saved
  if (present(tke_userdef_constants)) then
  tke_constants_out=> tke_userdef_constants
 end if

  select case(trim(varname))

    case('handle_old_vals')
    tke_constants_out%handle_old_vals=val
    case ('tke_mxl_choice') 
    tke_constants_out%tke_mxl_choice=val 
  end select
    
end subroutine cvmix_tke_put_tke_int

!=================================================================================

subroutine cvmix_tke_put_tke_logical(varname,val,tke_userdef_constants)
!This subroutine puts logicals to TKE variables
!IN
    character(len=*),           intent(in) :: varname
    logical,                    intent(in) :: val
!OUT   
    type(tke_type), intent(inout), target, optional:: tke_userdef_constants
    type(tke_type), pointer :: tke_constants_out

  tke_constants_out=>tke_constants_saved
  if (present(tke_userdef_constants)) then
  tke_constants_out=> tke_userdef_constants
 end if

  select case(trim(varname))

    case('only_tke')
      tke_constants_out%only_tke=val  
    case('use_ubound_dirichlet')
      tke_constants_out%use_ubound_dirichlet=val  
    case('use_lbound_dirichlet')
      tke_constants_out%use_lbound_dirichlet=val  
    case DEFAULT
      print*, "ERROR:", trim(varname), " not a valid choice"
      stop 1
  end select
    !!enable_GM etc can go in here
end subroutine cvmix_tke_put_tke_logical

!=================================================================================

subroutine cvmix_tke_put_tke_real(varname,val,tke_userdef_constants)
!This subroutine puts real values to TKE variables
!IN
    character(len=*),           intent(in) :: varname
    real(cvmix_r8),             intent(in) :: val
!OUT   
    type(tke_type), intent(inout), target, optional:: tke_userdef_constants
    type(tke_type), pointer :: tke_constants_out


  tke_constants_out=>tke_constants_saved
  if (present(tke_userdef_constants)) then
  tke_constants_out=> tke_userdef_constants
  end if

  select case(trim(varname))

    case('c_k') 
      tke_constants_out%c_k= val
    case('c_eps') 
      tke_constants_out%c_eps= val
    case('cd') 
      tke_constants_out%cd= val
    case('alpha_tke') 
      tke_constants_out%alpha_tke = val
    case('mxl_min') 
      tke_constants_out%mxl_min = val
    case('kappaM_min')
      tke_constants_out%kappaM_min = val
    case('kappaM_max')
      tke_constants_out%kappaM_max = val
    case('tke_min')
      tke_constants_out%tke_min = val    
    case('tke_surf_min')
      tke_constants_out%tke_surf_min = val    
    case DEFAULT
      print*, "ERROR:", trim(varname), " not a valid choice"
      stop 1

  end select

end subroutine cvmix_tke_put_tke_real

!=================================================================================


end module cvmix_tke 
