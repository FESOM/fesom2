module cvmix_idemix2
! This module contains the main computations of the IDEMIX 2 parameterization (described in "A Global Model for the Diapycnal
! Diffusivity Induced by Internal Gravity Waves", Olbers&Eden 2013) of Internal wave energy and its dissipation
!
! EdenOlbers2014: 
! --------------
! Eden, C. and Olbers, D. (2014): An Energy Compartment Model for Propagation, 
! Nonlinear Interaction, and Dissipation of Internal Gravity Waves , Journal of 
! Physical Oceanography, 44 (8), pp. 2093-2106 . doi: 10.1175/JPO-D-13-0224.1 

use cvmix_kinds_and_types,    only : cvmix_r8,                     &
                                     CVMIX_OVERWRITE_OLD_VAL,     &
                                     CVMIX_SUM_OLD_AND_NEW_VALS,  &
                                     CVMIX_MAX_OLD_AND_NEW_VALS,  &
                                     cvmix_data_type,             &
                                     cvmix_PI,                    & 
                                     cvmix_global_params_type
use cvmix_kinds_and_types_addon,    only : cvmix_data_type_addon
use cvmix_utils_addon,              only : cvmix_update_tke, solve_tridiag
implicit none

!_______________________________________________________________________________
private 
save
!public member functions
public :: cvmix_init_idemix2
public :: cvmix_compute_groupvel_idemix2
! public :: calc_idemix2_v0
! public :: cvmix_coeffs_idemix2
public :: gofx2, hofx1, hofx2


!_______________________________________________________________________________
! Interface to call the IDEMIX parameterization
! interface cvmix_coeffs_idemix2
!     module procedure integrate_idemix2  ! calculation ! FIXME: rename in cvmix_coeffs_low..
!     !module procedure idemix_wrap       ! necessary to handle old/new values and to hand over user_defined constants
! end interface cvmix_coeffs_idemix2

interface cvmix_init_idemix2
    module procedure init_idemix2  
end interface cvmix_init_idemix2

interface cvmix_compute_groupvel_idemix2
    module procedure compute_groupvel_idemix2  
end interface cvmix_compute_groupvel_idemix2

interface idemix2_put
    module procedure vmix_tke_put_idemix2_int
    module procedure vmix_tke_put_idemix2_real
    module procedure vmix_tke_put_idemix2_bool
end interface idemix2_put

!_______________________________________________________________________________
! types for Idemix
type, public :: idemix2_type
    private
    real(cvmix_r8) ::     &
    tau_v               , & ! time scale for vertical symmetrisation (sec)
    tau_h               , & ! time scale for horizontal symmetrisation, only necessary for lateral diffusion (sec)
    gamma               , & ! constant of order one derived from the shape of the spectrum in m space (dimensionless)
    jstar               , & ! spectral bandwidth in modes (dimensionless)
    mu0                     ! dissipation parameter (dimensionless)
    
    logical        ::     &
    enable_M2           , &
    enable_niw          , &
    enable_superbee_adv , &
    enable_AB_timestep  , &
    enable_hor_diffusion
    
    integer        ::     &
    handle_old_vals     , &
    nfbin               , &
    hor_diffusion_niter
    
end type idemix2_type
type(idemix2_type), target :: idemix2_constants_saved 

!_______________________________________________________________________________
contains



!
!
!_______________________________________________________________________________
subroutine init_idemix2(tau_v               , & ! time scale for vertical symmetrisation
                        tau_h               , & ! time scale for horizontal symmetrisation
                        gamma               , & ! const. derived from the shape of the spectrum in m spac
                        jstar               , & ! mode number scale
                        mu0                 , & ! dissipation parameter (dimensionless)
                        nfbin               , & ! number of spectral bins
                        enable_M2           , &
                        enable_niw          , & 
                        enable_superbee_adv , &
                        enable_AB_timestep  , &
                        enable_hor_diffusion, &
                        hor_diffusion_niter , & 
                        handle_old_vals, idemix2_userdef_constants)

    ! This subroutine sets user or default values for IDEMIX parameters
    real(cvmix_r8), optional, intent(in) :: tau_v, tau_h, gamma, jstar, mu0
    logical       , optional, intent(in) :: enable_M2, enable_niw, &
                                            enable_superbee_adv, enable_AB_timestep, &
                                            enable_hor_diffusion
    integer       , optional, intent(in) :: nfbin, hor_diffusion_niter
    integer       , optional, intent(in) :: handle_old_vals
    type(idemix2_type), target, optional, intent(inout) :: idemix2_userdef_constants
    

    ! FIXME: not sure about the allowed ranges for idemix parameters, default values confirm with pyOM testcases
    !___ tau_v _________________________________________________________________ 
    if (present(tau_v)) then
        if(tau_v.lt.1.d0*86400.0 .or. tau_v .gt. 100.d0*86400.0) then
            print*, "ERROR:tau_v can only be allowed_range"
            stop 1
        end if
        call idemix2_put('tau_v', tau_v, idemix2_userdef_constants)
    else
        call idemix2_put('tau_v',1._cvmix_r8*86400.0 , idemix2_userdef_constants)
    end if

    !___ tau_h _________________________________________________________________
    if (present(tau_h)) then
        if(tau_h.lt. 0.01*864000. .or. tau_h .gt. 100.*86400.) then
            print*, "ERROR:tau_h can only be allowed_range"
            stop 1
        end if
        call idemix2_put('tau_h', tau_h, idemix2_userdef_constants)
    else
        call idemix2_put('tau_h', 15._cvmix_r8*86400.0, idemix2_userdef_constants)
    end if

    !___ gamma _________________________________________________________________
    if (present(gamma)) then
        if(gamma.lt. 1.d0 .or. gamma .gt. 3.d0) then
            print*, "ERROR:gamma can only be allowed_range"
            stop 1
        end if
        call idemix2_put('gamma', gamma, idemix2_userdef_constants)
    else
        call idemix2_put('gamma', 1.57_cvmix_r8, idemix2_userdef_constants)
    end if

    !___ jstar _________________________________________________________________
    if (present(jstar)) then
        if(jstar.lt. 5.d0 .or. jstar .gt. 15.d0) then
            print*, "ERROR:jstar can only be allowed_range"
            stop 1
        end if
        call idemix2_put('jstar', jstar, idemix2_userdef_constants)
    else
        call idemix2_put('jstar', 10._cvmix_r8 , idemix2_userdef_constants)
    end if

    !___ mu0 ___________________________________________________________________
    if (present(mu0)) then
        if(mu0.lt. 0.d0 .or. mu0 .gt. 3.d0) then
            print*, "ERROR: mu0 can only be allowed_range"
            stop 1
        end if
        call idemix2_put('mu0', mu0, idemix2_userdef_constants)
    else
        call idemix2_put('mu0', 4._cvmix_r8/3.0 , idemix2_userdef_constants)
    end if
    
    !___ superbee_adv __________________________________________________________
    if (present(nfbin)) then
        call idemix2_put('nfbin', nfbin, idemix2_userdef_constants)
    else
        call idemix2_put('nfbin', 52 , idemix2_userdef_constants)
    end if
    
    !___ enable_M2 _____________________________________________________________
    if (present(enable_M2)) then
        call idemix2_put('enable_M2', enable_M2, idemix2_userdef_constants)
    else
        call idemix2_put('enable_M2', .true. , idemix2_userdef_constants)
    end if
    
    !___ enable_niw ____________________________________________________________
    if (present(enable_niw)) then
        call idemix2_put('enable_niw', enable_niw, idemix2_userdef_constants)
    else
        call idemix2_put('enable_niw', .true. , idemix2_userdef_constants)
    end if
    
    !___ superbee_adv __________________________________________________________
    if (present(enable_superbee_adv)) then
        call idemix2_put('enable_superbee_adv', enable_superbee_adv, idemix2_userdef_constants)
    else
        call idemix2_put('enable_superbee_adv', .true. , idemix2_userdef_constants)
    end if
    
    !___ superbee_adv __________________________________________________________
    if (present(enable_AB_timestep)) then
        call idemix2_put('enable_AB_timestep', enable_AB_timestep, idemix2_userdef_constants)
    else
        call idemix2_put('enable_AB_timestep', .true. , idemix2_userdef_constants)
    end if
    
    !___ homogenous diffusion idemix1 functionality ____________________________
    if (present(enable_hor_diffusion)) then
        call idemix2_put('enable_hor_diffusion', enable_hor_diffusion, idemix2_userdef_constants)
    else
        call idemix2_put('enable_hor_diffusion', .false. , idemix2_userdef_constants)
    end if
    
    if (present(hor_diffusion_niter)) then
        call idemix2_put('hor_diffusion_niter', hor_diffusion_niter, idemix2_userdef_constants)
    else
        call idemix2_put('hor_diffusion_niter', 5 , idemix2_userdef_constants)
    end if
    
    !___ handle_old_vals________________________________________________________
    if (present(handle_old_vals)) then
        if(handle_old_vals.lt. 1 .or. handle_old_vals.gt. 3 ) then
            print*, "ERROR:handle_old_vals can only be 1 to 3"
            stop 1
        end if
        call idemix2_put('handle_old_vals', handle_old_vals, idemix2_userdef_constants)
    else
        call idemix2_put('handle_old_vals', 1, idemix2_userdef_constants)
    end if
end subroutine init_idemix2


!
!
!_______________________________________________________________________________
! Compute idemix representative vertical (v0) and horizontal (c0) group velocites 
! as well as the enery dissipation coefficient alpha_c
subroutine compute_groupvel_idemix2(nlev        , &
                                 nfbin       , &   
                                 dtime       , &
                                 coriolis    , & 
                                 grady_coriol, &
                                 coslat      , &
                                 Nsqr        , & 
                                 omega_M2    , &
                                 omega_niw   , & 
                                 cn          , &
                                 cn_gradx    , &
                                 cn_grady    , &
                                 phit        , &
                                 phiu        , &
                                 alpha_c     , & 
                                 c0          , & 
                                 v0          , &
                                 cg_M2       , &
                                 cg_niw      , &
                                 u_M2        , &
                                 v_M2        , &
                                 w_M2        , &
                                 u_niw       , &
                                 v_niw       , &
                                 w_niw       , &
                                 idemix2_const_userdef &
                                )
    !___Input___________________________________________________________________
    type(idemix2_type), intent(in) , optional, target :: idemix2_const_userdef
    integer           , intent(in)                    :: nlev           , &
                                                         nfbin
    real(cvmix_r8)    , intent(in)                    :: dtime          , &
                                                         coriolis       , &
                                                         grady_coriol   , &
                                                         coslat         , &
                                                         cn             , &
                                                         cn_gradx       , &
                                                         cn_grady       , &
                                                         omega_M2       , &
                                                         omega_niw
    real(cvmix_r8)    , intent(in) , dimension(nlev+1):: Nsqr 
    real(cvmix_r8)    , intent(in) , dimension(nfbin) :: phit           , &
                                                         phiu 
    !___Output__________________________________________________________________
    real(cvmix_r8)    , intent(out), dimension(nlev)  :: alpha_c        , &
                                                         c0             , &
                                                         v0
    real(cvmix_r8)    , intent(out)                   :: cg_M2          , &
                                                         cg_niw
    real(cvmix_r8)    , intent(out), dimension(nfbin) :: u_M2           , &
                                                         v_M2           , & 
                                                         w_M2           , & 
                                                         u_niw          , &
                                                         v_niw          , &
                                                         w_niw
    !___Local___________________________________________________________________
    integer                                           :: di, fbin_i
    real(cvmix_r8)                                    :: fxa, intNz, cstar, &
                                                         kdot_x_M2, kdot_x_niw, &
                                                         kdot_y_M2, kdot_y_niw
    type(idemix2_type), pointer                       :: idemix2_const_in
    
    ! do pointer into save variable or into user defined input variable 
    idemix2_const_in => idemix2_constants_saved
    if (present(idemix2_const_userdef)) then
        idemix2_const_in => idemix2_const_userdef
    end if
    
    !___________________________________________________________________________
    ! calculate cstar from OE13 Eq. (13) --> vertical integrate sqrt(N²) --> which 
    ! is already in cn
    ! --> cstar = max(1d-2,intNz/(cvmix_PI*idemix2_const_in%jstar) )
    cstar = max(1d-2,cn/(idemix2_const_in%jstar) )
    
    !___________________________________________________________________________
    ! calculate vertical and horizontal representative group velocities c0 and v0
    ! c0: OE13 Eq. (13)   
    ! alpha_c iwe**2: dissipation of internal wave energy (OE13 Eq. (15))
    do di=1,nlev
        fxa         = sqrt(max(0d0,Nsqr(di))) / ( 1d-22 + abs(coriolis) )
        c0(di)      = max(0d0, idemix2_const_in%gamma*cstar*gofx2(fxa))
        v0(di)      = max(0d0, idemix2_const_in%gamma*cstar*hofx2(fxa))
        alpha_c(di) = max( 1d-4, idemix2_const_in%mu0*acosh(max(1d0,fxa))*abs(coriolis)/cstar**2 )
        
        ! --> ATTENTION: CHECK EFFECT OF THESE LINES !!!
        !! set v0 to zero to prevent horizontal iwe propagation in mixed layer
        !if ( fxa<1d0 ) then
        !   v0(di) = 0d0
        !endif
    enddo
    ! --> ATTENTION: REMEMBER TO CHECK HORIZONTAL STABILTY CRITERIA WHEN 
    !     HORIZOTNAL OPERATION IS DONE !!!
    ! --> pyOM2:    
    !     if (enable_idemix_hor_diffusion .or. enable_idemix_hor_diffusion_iter) then
    !     ! check for stability criterium, lateral diffusion is explicit
    !     !  tau_h v0^2 *dt/dx^2 <= 0.5  ->   v0  <  sqrt( 0.5*dx^2/(dt tau_h)  )
    !     do j=js_pe-onx,je_pe+onx
    !       do i=is_pe-onx,ie_pe+onx
    !          if (enable_idemix_hor_diffusion_iter) then
    !             fxa = 0.2*min( dxt(i)*cost(j), dyt(j) )**2/ max(1D0,dt_tracer/idemix_hor_diffusion_iter *tau_h )
    !          else
    !             fxa = 0.2*min( dxt(i)*cost(j), dyt(j) )**2/ max(1D0,dt_tracer *tau_h )
    !          endif
    !          v0(i,j,:) = min( sqrt(fxa), v0(i,j,:) )
    !       enddo
    !     enddo
    !   endif

    !___________________________________________________________________________
    ! compute group velocity of M2 internal tidal waves
    if (idemix2_const_in%enable_M2) then
        cg_M2=sqrt( max(0d0, omega_M2**2 - coriolis**2 )  ) * cn/omega_M2
    endif

    ! compute group velocity of of near-inertial waves (NIW)
    if (idemix2_const_in%enable_niw) then
        cg_niw=sqrt(  max(0d0, omega_niw**2 - coriolis**2 )  ) * cn/omega_niw
    endif
    
    !___________________________________________________________________________
    ! eq4 in EdenOlbers2014:
    ! 
    ! phi_dot =  [c_n/omega/sqrt(omega²-f²)) * f * grad_f 
    !             + sqrt(omega²-f²)/omega * grad_cn ] * (sin(phi), -cos(phi))
    ! 
    ! phi_dot =  (k_dot_x, kdot_y) * (sin(phi), cos(phi))
    ! 
    
    ! compute wave vector component for the M2 internal tides: 
    ! kdot_M2 = (kdot_x_M2, kdot_y_M2)
    if (idemix2_const_in%enable_M2) then
        ! compute: sqrt(omega²-f²)
        fxa = max(1d-10,omega_M2**2 - coriolis**2 )
        
        ! 1st part coriolis contribution
        ! compute: c_n/omega/sqrt(omega²-f²)) * f * grad_f 
        kdot_y_M2 = -cn/sqrt(fxa)*coriolis/omega_M2*grady_coriol
        !           |
        !           +-> this minus sign is from the -cos(phi)
        !
        ! 2nd part topographic/buoyancy driven contribution
        ! compute: sqrt(omega²-f²)/omega * grad_cn
        kdot_y_M2 = kdot_y_M2 - sqrt(fxa)/omega_M2*cn_grady
        !                     |
        !                     +-> this minus sign is from the -cos(phi)
        kdot_x_M2 = sqrt(fxa)/omega_M2*cn_gradx
    endif
    
    ! compute wave vector component for the near-inertial-waves: 
    ! kdot_niw = (kdot_x_niw, kdot_y_niw)
    if (idemix2_const_in%enable_niw) then
        ! compute: sqrt(omega²-f²)
        fxa = max(1d-10,omega_niw**2 - coriolis**2 )
        
        ! 1st part coriolis contribution
        ! compute: c_n/omega/sqrt(omega²-f²)) * f * grad_f 
        kdot_y_niw = -cn/sqrt(fxa)*coriolis/omega_niw*grady_coriol
        !            |
        !            +-> this minus sign is from the -cos(phi)
        !
        ! 2nd part topographic/buoyancy driven contribution
        ! compute: sqrt(omega²-f²)/omega * grad_cn
        kdot_y_niw = kdot_y_niw - sqrt(fxa)/omega_niw*cn_grady
        !                       |
        !                       +-> this minus sign is from the -cos(phi)
        kdot_x_niw = sqrt(fxa)/omega_niw*cn_gradx
    endif
    
    
    !___________________________________________________________________________
    !zonal, meridional and vertial component of M2 internal tide group velocity
    if (idemix2_const_in%enable_M2) then
        do fbin_i=2,nfbin-1
            u_M2(fbin_i) = cg_M2*cos( phit(fbin_i) )
            v_M2(fbin_i) = cg_M2*sin( phit(fbin_i) ) * coslat 
            
            !!! ATTENTION !!! --> this might be needed to move to nodes
            w_M2(fbin_i) = (kdot_y_M2*cos(phiu(fbin_i)) + kdot_x_M2*sin(phiu(fbin_i)) )
        enddo
    end if 
    
    !zonal, meridional and vertial component of NIW internal wave group velocity
    if (idemix2_const_in%enable_niw) then
        do fbin_i=2,nfbin-1
            u_niw(fbin_i) = cg_niw*cos( phit(fbin_i) )
            v_niw(fbin_i) = cg_niw*sin( phit(fbin_i) ) * coslat 
            
            !!! ATTENTION !!! --> this might be needed to move to nodes
            w_niw(fbin_i) = (kdot_y_niw*cos(phiu(fbin_i)) + kdot_x_niw*sin(phiu(fbin_i)) )
        enddo
    end if
end subroutine compute_groupvel_idemix2    
    
    
! !
! !
! !_______________________________________________________________________________
! subroutine set_compute_groupvel_idemix2(dzt        , &
!                                      nlev       , &
!                                      dtime      , &
!                                      coriolis   , & 
!                                      Nsqr       , &
!                                      omega_M2   , &
!                                      omega_niw  , &
!                                      alpha_c    , & 
!                                      c0         , & 
!                                      v0         , &
!                                      idemix2_const_userdef &
!                                     )    
                                    
                                    
!
!
!_______________________________________________________________________________
subroutine vmix_tke_put_idemix2_real(varname, val, idemix2_userdef_constants)
    ! This subroutine puts real values to IDEMIX variables
    character(len=*),          intent(in) :: varname
    real(cvmix_r8),            intent(in) :: val
    type(idemix2_type), intent(inout), target, optional:: idemix2_userdef_constants
    type(idemix2_type), pointer            :: idemix2_constants_out
    ! do pointer into save variable 
    idemix2_constants_out=>idemix2_constants_saved
    ! if input idemix_userdef_constants present do pointer into this 
    if (present(idemix2_userdef_constants)) idemix2_constants_out=> idemix2_userdef_constants
    select case(trim(varname))
        case('tau_v') ; idemix2_constants_out%tau_v = val
        case('tau_h') ; idemix2_constants_out%tau_h = val
        case('jstar') ; idemix2_constants_out%jstar = val
        case('gamma') ; idemix2_constants_out%gamma = val
        case('mu0'  ) ; idemix2_constants_out%mu0   = val
        case DEFAULT
            print*, "ERROR:", trim(varname), " not a valid choice"
            stop 1
    end select
end subroutine vmix_tke_put_idemix2_real

subroutine vmix_tke_put_idemix2_int(varname,val,idemix2_userdef_constants)
    ! This subroutine puts integer values to IDEMIX variables
    character(len=*),           intent(in) :: varname
    integer,                    intent(in) :: val
    type(idemix2_type), intent(inout), target, optional:: idemix2_userdef_constants
    type(idemix2_type), pointer :: idemix2_constants_out
    idemix2_constants_out=>idemix2_constants_saved
    if (present(idemix2_userdef_constants)) idemix2_constants_out=> idemix2_userdef_constants
    select case(trim(varname))
        case('handle_old_vals'    ) ; idemix2_constants_out%handle_old_vals    = val
        case('nfbin'              ) ; idemix2_constants_out%nfbin              = val
        case('hor_diffusion_niter') ; idemix2_constants_out%hor_diffusion_niter= val
        case DEFAULT
            print*, "ERROR:", trim(varname), " not a valid choice"
            stop 1
    end select
end subroutine vmix_tke_put_idemix2_int

subroutine vmix_tke_put_idemix2_bool(varname,val,idemix2_userdef_constants)
    ! This subroutine puts integer values to IDEMIX variables
    character(len=*),           intent(in) :: varname
    logical,                    intent(in) :: val
    type(idemix2_type), intent(inout), target, optional:: idemix2_userdef_constants
    type(idemix2_type), pointer :: idemix2_constants_out
    idemix2_constants_out=>idemix2_constants_saved
    if (present(idemix2_userdef_constants)) idemix2_constants_out=> idemix2_userdef_constants
    select case(trim(varname))
        case('enable_M2'           ) ; idemix2_constants_out%enable_M2           = val
        case('enable_niw'          ) ; idemix2_constants_out%enable_niw          = val
        case('enable_superbee_adv' ) ; idemix2_constants_out%enable_superbee_adv = val
        case('enable_AB_timestep'  ) ; idemix2_constants_out%enable_AB_timestep  = val
        case('enable_hor_diffusion') ; idemix2_constants_out%enable_hor_diffusion= val
        case DEFAULT
            print*, "ERROR:", trim(varname), " not a valid choice"
            stop 1
    end select
end subroutine vmix_tke_put_idemix2_bool


!
!
!_______________________________________________________________________________
! function g(x) --> adapted from pyOM 
function gofx2(x1)
    implicit none
    real(cvmix_r8) :: gofx2,x1,x2,c
    real(cvmix_r8), parameter :: pi = 3.14159265358979323846264338327950588
    x2=max(3d0,x1)
    c= 1.-(2./pi)*asin(1./x2)
    gofx2 = 2/pi/c*0.9*x2**(-2./3.)*(1-exp(-x2/4.3))
end function gofx2

! function h(x) --> adapted from pyOM
function hofx2(x1)
    implicit none
    real(cvmix_r8) :: hofx2,x1,x2
    real(cvmix_r8), parameter :: pi = 3.14159265358979323846264338327950588
    x2 = max(1d1, x1) ! by_nils: it has to be x2>1
    hofx2 = (2./pi)/(1.-(2./pi)*asin(1./x2)) * (x2-1.)/(x2+1.)
end function hofx2

! function h(x) --> from pyOM
function hofx1(x)
    implicit none
    real(cvmix_r8) :: hofx1,x
    real(cvmix_r8), parameter :: pi = 3.14159265358979323846264338327950588
    hofx1 = (2./pi)/(1.-(2./pi)*asin(1./x)) * (x-1.)/(x+1.)
end function hofx1

end module cvmix_idemix2