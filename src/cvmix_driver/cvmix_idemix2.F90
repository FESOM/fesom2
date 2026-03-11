module cvmix_idemix2
! This module contains the main computations of the IDEMIX 2 parameterization (described in "A Global Model for the Diapycnal
! Diffusivity Induced by Internal Gravity Waves", Olbers&Eden 2013) of Internal wave energy and its dissipation
!
! EdenOlbers2014: 
! --------------
! Eden, C. and Olbers, D. (2014): An Energy Compartment Model for Propagation, 
! Nonlinear Interaction, and Dissipation of Internal Gravity Waves , Journal of 
! Physical Oceanography, 44 (8), pp. 2093-2106 . doi: 10.1175/JPO-D-13-0224.1 

use cvmix_kinds_and_types,       only : cvmix_r8                   , &
                                        CVMIX_OVERWRITE_OLD_VAL    , &
                                        CVMIX_SUM_OLD_AND_NEW_VALS , &
                                        CVMIX_MAX_OLD_AND_NEW_VALS , &
                                        cvmix_data_type            , &
                                        cvmix_PI                   , & 
                                        cvmix_global_params_type
use cvmix_kinds_and_types_addon, only : cvmix_data_type_addon
use cvmix_utils_addon,           only : cvmix_update_tke, solve_tridiag
implicit none

!_______________________________________________________________________________
private 
save
!public member functions
public :: cvmix_idemix2_init
public :: cvmix_idemix2_compute_param
public :: cvmix_idemix2_compute_compart_groupvel
public :: cvmix_idemix2_compute_compart_interact_tscale
public :: cvmix_idemix2_compute_M2_dissipation
public :: cvmix_idemix2_compute_vert_struct_fct
public :: cvmix_idemix2_compute_vdiff_vdiss_Eiw
public :: cvmix_idemix2_compute_Eiw_waveinteract
public :: cvmix_idemix2_compute_Eiw_diss2KvAv
public :: gofx2, hofx1, hofx2


!_______________________________________________________________________________
! Interface to call the IDEMIX parameterization
! interface cvmix_coeffs_idemix2
!     module procedure integrate_idemix2  ! calculation ! FIXME: rename in cvmix_coeffs_low..
!     !module procedure idemix_wrap       ! necessary to handle old/new values and to hand over user_defined constants
! end interface cvmix_coeffs_idemix2

interface cvmix_idemix2_init
    module procedure compute_init
end interface cvmix_idemix2_init

interface cvmix_idemix2_compute_param
    module procedure compute_param
end interface cvmix_idemix2_compute_param

interface cvmix_idemix2_compute_compart_groupvel
    module procedure compute_compart_groupvel
end interface cvmix_idemix2_compute_compart_groupvel

interface cvmix_idemix2_compute_compart_interact_tscale
    module procedure compute_compart_interact_tscale
end interface cvmix_idemix2_compute_compart_interact_tscale

interface cvmix_idemix2_compute_M2_dissipation
    module procedure compute_M2_dissipation
end interface cvmix_idemix2_compute_M2_dissipation

interface cvmix_idemix2_compute_vert_struct_fct
    module procedure compute_vert_struct_fct
end interface cvmix_idemix2_compute_vert_struct_fct   

interface cvmix_idemix2_compute_vdiff_vdiss_Eiw
    module procedure compute_vdiff_vdiss_Eiw
end interface cvmix_idemix2_compute_vdiff_vdiss_Eiw

interface cvmix_idemix2_compute_Eiw_waveinteract
    module procedure compute_Eiw_waveinteract
end interface cvmix_idemix2_compute_Eiw_waveinteract
    
interface cvmix_idemix2_compute_Eiw_diss2KvAv
    module procedure compute_Eiw_diss2KvAv
end interface cvmix_idemix2_compute_Eiw_diss2KvAv    
    
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
    mu0                 , & ! dissipation parameter (dimensionless)
    shelf_dist              ! shelf definition in meters from coast
    
    logical        ::     &
    enable_M2           , &
    enable_niw          , &
    enable_superbee_adv , &
    enable_AB_timestep  , &
    enable_hor_diffusion, &
    enable_hor_diff_iter
    
    integer        ::     &
    handle_old_vals     , &
    nfbin               , &
    hor_diff_niter
    
end type idemix2_type
type(idemix2_type), target :: idemix2_constants_saved 

!_______________________________________________________________________________
contains



!
!
!_______________________________________________________________________________
subroutine compute_init(tau_v               , & ! time scale for vertical symmetrisation
                        tau_h               , & ! time scale for horizontal symmetrisation
                        gamma               , & ! const. derived from the shape of the spectrum in m spac
                        jstar               , & ! mode number scale
                        mu0                 , & ! dissipation parameter (dimensionless)
                        nfbin               , & ! number of spectral bins
                        shelf_dist          , & ! shelf definition in meters from coast
                        enable_M2           , &
                        enable_niw          , & 
                        enable_superbee_adv , &
                        enable_AB_timestep  , &
                        enable_hor_diffusion, &
                        enable_hor_diff_iter, &
                        hor_diff_niter      , & 
                        handle_old_vals, idemix2_userdef_constants)

    ! This subroutine sets user or default values for IDEMIX parameters
    real(cvmix_r8), optional, intent(in) :: tau_v, tau_h, gamma, jstar, mu0, shelf_dist
    logical       , optional, intent(in) :: enable_M2, enable_niw, &
                                            enable_superbee_adv, enable_AB_timestep, &
                                            enable_hor_diffusion, enable_hor_diff_iter
    integer       , optional, intent(in) :: nfbin, hor_diff_niter
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
    
    !___ shelf_dist ____________________________________________________________
    if (present(shelf_dist)) then
        call idemix2_put('shelf_dist', shelf_dist, idemix2_userdef_constants)
    else
        call idemix2_put('shelf_dist', 300000._cvmix_r8, idemix2_userdef_constants)
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
    
    if (present(enable_hor_diff_iter)) then
        call idemix2_put('enable_hor_diff_iter', enable_hor_diff_iter, idemix2_userdef_constants)
    else
        call idemix2_put('enable_hor_diff_iter', .false. , idemix2_userdef_constants)
    end if
    
    if (present(hor_diff_niter)) then
        call idemix2_put('hor_diff_niter', hor_diff_niter, idemix2_userdef_constants)
    else
        call idemix2_put('hor_diff_niter', 5 , idemix2_userdef_constants)
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
end subroutine compute_init


!
!
!_______________________________________________________________________________
! Compute idemix2 representative vertical (v0) and horizontal (c0) group velocites 
! as well as the enery dissipation coefficient alpha_c
subroutine compute_param(             & 
                nlev                , &
                coriolis            , & 
                Nsqr                , & 
                cn                  , &
                alpha_c             , & 
                c0                  , & 
                v0                  , &
                idemix2_const_userdef &
                )
    !___Input___________________________________________________________________
    type(idemix2_type), intent(in) , optional, target :: idemix2_const_userdef
    integer           , intent(in)                    :: nlev           
    real(cvmix_r8)    , intent(in)                    :: coriolis       , &
                                                         cn             
                                                         
    real(cvmix_r8)    , intent(in) , dimension(nlev+1):: Nsqr 
    
    !___Output__________________________________________________________________
    real(cvmix_r8)    , intent(out), dimension(nlev)  :: alpha_c        , &
                                                         c0             , &
                                                         v0
    
    !___Local___________________________________________________________________
    integer                                           :: di
    real(cvmix_r8)                                    :: fxa, cstar !, intNz
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
end subroutine compute_param
    
    
!
!
!_______________________________________________________________________________
! Compute idemix representative vertical (v0) and horizontal (c0) group velocites 
! as well as the enery dissipation coefficient alpha_c
subroutine compute_compart_groupvel(          &
                nfbin               , &   
                coriolis            , & 
                grady_coriol        , &
                coslat              , &
                cn                  , &
                cn_gradx            , &
                cn_grady            , &
                phit                , &
                phiu                , &
                omega_compart       , &
                u_compart           , &
                v_compart           , &
                w_compart           , &
                idemix2_const_userdef &
                )
                            
    !___Input___________________________________________________________________
    type(idemix2_type), intent(in) , optional, target :: idemix2_const_userdef
    integer           , intent(in)                    :: nfbin
    real(cvmix_r8)    , intent(in)                    :: coriolis       , &
                                                         grady_coriol   , &
                                                         coslat         , &
                                                         cn             , &
                                                         cn_gradx       , &
                                                         cn_grady       , &
                                                         omega_compart    
    real(cvmix_r8)    , intent(in) , dimension(nfbin) :: phit           , &
                                                         phiu 
    !___Output__________________________________________________________________
    real(cvmix_r8)    , intent(out), dimension(nfbin) :: u_compart      , &
                                                         v_compart      , &
                                                         w_compart
    !___Local___________________________________________________________________
    integer                                           :: di, fbin_i
    real(cvmix_r8)                                    :: fxa, intNz, cstar, &
                                                         kdot_x, kdot_y, cg
    type(idemix2_type), pointer                       :: idemix2_const_in
    
    ! do pointer into save variable or into user defined input variable 
    idemix2_const_in => idemix2_constants_saved
    if (present(idemix2_const_userdef)) then
        idemix2_const_in => idemix2_const_userdef
    end if
    
    !___________________________________________________________________________
    ! compute group velocity of M2 internal tidal waves
    cg=sqrt( max(0d0, omega_compart**2 - coriolis**2 )  ) * cn/omega_compart
    
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
    
    ! compute: sqrt(omega²-f²)
    fxa = max(1d-10,omega_compart**2 - coriolis**2 )
        
    ! 1st part coriolis contribution
    ! compute: c_n/omega/sqrt(omega²-f²)) * f * grad_f 
    kdot_y = -cn/sqrt(fxa)*coriolis/omega_compart*grady_coriol
    !           |
    !           +-> this minus sign is from the -cos(phi)
    !
    ! 2nd part topographic/buoyancy driven contribution
    ! compute: sqrt(omega²-f²)/omega * grad_cn
    kdot_y = kdot_y - sqrt(fxa)/omega_compart*cn_grady
    !                     |
    !                     +-> this minus sign is from the -cos(phi)
    kdot_x = sqrt(fxa)/omega_compart*cn_gradx
    
    
    !___________________________________________________________________________
    !zonal, meridional and cross-spectral component of M2 internal tide group velocity
    do fbin_i=2,nfbin-1
        u_compart(fbin_i) = cg*cos( phit(fbin_i) )
        v_compart(fbin_i) = cg*sin( phit(fbin_i) ) * coslat 
        w_compart(fbin_i) = (kdot_y*cos(phiu(fbin_i)) + kdot_x*sin(phiu(fbin_i)) )
    end do
    
end subroutine compute_compart_groupvel
    
    

!
!
!_______________________________________________________________________________
subroutine compute_compart_interact_tscale(   &
                dtime               , &
                lat                 , &
                coriolis            , &
                cn                  , &
                zbottom             , &
                topo_hrms           , &
                topo_hlam           , &
                topo_shelf          , &
                omega_compart       , &
                tau_compart         , &
                idemix2_const_userdef &
                )    
    !___Input___________________________________________________________________
    type(idemix2_type), intent(in) , optional, target :: idemix2_const_userdef
    real(cvmix_r8)    , intent(in)                    :: dtime          , &
                                                         lat            , &
                                                         coriolis       , &
                                                         omega_compart  , &
                                                         cn             , &
                                                         zbottom        , &
                                                         topo_hrms      , &
                                                         topo_hlam        
    logical           , intent(in)                    :: topo_shelf
    !___Output__________________________________________________________________
    real(cvmix_r8)    , intent(out)                   :: tau_compart
    !___Local___________________________________________________________________
    real(cvmix_r8)                                    :: fxb, fxc, N0
    type(idemix2_type), pointer                       :: idemix2_const_in
    
    ! do pointer into save variable or into user defined input variable 
    idemix2_const_in => idemix2_constants_saved
    if (present(idemix2_const_userdef)) then
        idemix2_const_in => idemix2_const_userdef
    end if
    
    !___________________________________________________________________________
    ! compute NIW Dissipation Timescale
    if (zbottom>0) then
        N0=cn*cvmix_PI/zbottom
        if (  N0> abs(coriolis) .and. omega_compart> abs(coriolis)  ) then
            fxc = topo_hrms**2.0 * 2.0*cvmix_PI/(1d-12+topo_hlam) 
            fxb = 0.5* N0*( (omega_compart**2.0 + coriolis**2.0)/omega_compart**2.0 )**2.0  &
                            *(omega_compart**2.0 - coriolis**2.0)**0.5/omega_compart
            tau_compart = min(0.5/dtime, fxc*fxb/zbottom) 
        endif
    endif
        
    ! on shelf limitation 
    if (topo_shelf) tau_compart  = 1./(3.*86400.0)
    tau_compart = max(1d0/(50.*86400.0), tau_compart )
    
end subroutine compute_compart_interact_tscale                                     


!
!
!_______________________________________________________________________________
subroutine compute_M2_dissipation(    &
                lat                 , &
                cn                  , &
                zbottom             , &
                alpha_M2_c          , &
                idemix2_const_userdef &
                )    
    !___Input___________________________________________________________________
    type(idemix2_type), intent(in) , optional, target :: idemix2_const_userdef
    real(cvmix_r8)    , intent(in)                    :: lat            , &
                                                         cn             , &
                                                         zbottom        
                                                         
    !___Output__________________________________________________________________
    real(cvmix_r8)    , intent(out)                   :: alpha_M2_c
    !___Local___________________________________________________________________
    real(cvmix_r8)                                    :: N0
    real(cvmix_r8)                                    :: mstar=0.01 ,&
                                                         M2_f=2*cvmix_PI/(12.42*60*60)
    type(idemix2_type), pointer                       :: idemix2_const_in
    
    ! do pointer into save variable or into user defined input variable 
    idemix2_const_in => idemix2_constants_saved
    if (present(idemix2_const_userdef)) then
        idemix2_const_in => idemix2_const_userdef
    end if
    
    ! compute alpha_M2_c M2 Continuous Dissipation Rate (Background M2 tidal 
    ! dissipation rate (not topographic))
    alpha_M2_c = 0.0
    if (zbottom > 0.) then 
        N0 = cn*cvmix_PI/zbottom+1D-20
        if (abs(lat) < 28.5 ) then
            ! lambda+/M2=15*E*mstar/N * (sin(phi-28.5)/sin(28.5))^1/2
            alpha_M2_c=alpha_M2_c+ M2_f*15*mstar/N0* (sin( abs(abs(lat) -28.5)/180.*cvmix_PI )/sin(28.5/180.*cvmix_PI) )**0.5
        endif
        if (abs(lat) < 74.5 ) then
            ! lambda-/M2 =  0.7*E*mstar/N *sin^2(phi)
            alpha_M2_c=alpha_M2_c+ M2_f*0.7*mstar/N0* sin( abs(lat)/180.*cvmix_PI )**2
        endif
        alpha_M2_c=alpha_M2_c/zbottom
    endif
        
    ! ensure background value
    alpha_M2_c = max(0D0, min( 1d-5,alpha_M2_c ))

end subroutine compute_M2_dissipation                     


!
!
!_______________________________________________________________________________
! Calculates the vertical structure functions for internal wave modes in the 
! IDEMIX model. These functions describe how wave energy is distributed vertically 
! in the ocean.
! 
! Key Variables:
!   phin: Vertical structure function φₙ(z) - the modal function itself
!   phinz: Vertical derivative dφₙ/dz - the gradient of the modal function
!   E_struct_M2: Energy structure function for M2 tides
!   E_struct_niw: Energy structure function for near-inertial waves
subroutine compute_vert_struct_fct(   &
                nlev                , &
                dzw                 , &
                coriolis            , &
                omega_M2            , &
                omega_niw           , &
                cn                  , &
                Nsqr                , &
                E_struct_M2         , &
                E_struct_niw        , &
                idemix2_const_userdef &
                )
    !___Input___________________________________________________________________
    type(idemix2_type), intent(in) , optional, target :: idemix2_const_userdef
    integer           , intent(in)                    :: nlev           
    real(cvmix_r8)    , intent(in)                    :: dzw(:)         
    real(cvmix_r8)    , intent(in)                    :: coriolis       , &
                                                         omega_M2       , &
                                                         omega_niw      , &
                                                         cn             
    real(cvmix_r8)    , intent(in)                    :: Nsqr(:) 
    
    !___Output__________________________________________________________________
    real(cvmix_r8)    , intent(inout), optional       :: E_struct_M2(:)    , &
                                                         E_struct_niw(:)   
                                                         
    !___Local___________________________________________________________________
    integer                                           :: di
    real(cvmix_r8)                                    :: fxa, norm, phin_dim1
    real(cvmix_r8), dimension(:)                      :: Nzw(nlev), phin(nlev), dphindz(nlev)
    type(idemix2_type), pointer                       :: idemix2_const_in            
                
    ! do pointer into save variable or into user defined input variable 
    idemix2_const_in => idemix2_constants_saved
    if (present(idemix2_const_userdef)) then
        idemix2_const_in => idemix2_const_userdef
    end if
    
    !___________________________________________________________________________
    ! calculate int_(-h)^z N dz
    phin(:)    = 0.0_cvmix_r8
    Nzw(1)     = max(0.0_cvmix_r8, (Nsqr(1)+Nsqr(2)) * 0.5_cvmix_r8)
    phin(1)    = sqrt(Nzw(1))*dzw(1)
    phin_dim1  = phin(1)
    fxa        = phin(1)/cn
    dphindz(1) = sin(fxa)/Nzw(1)**0.25_cvmix_r8 
    phin(1)    = cos(fxa)*Nzw(1)**0.25_cvmix_r8
    
    do di = 2, nlev
        ! interpolate N² onto tracer levels and cumulative integrate
        Nzw(di)     = max(0.0_cvmix_r8, (Nsqr(di)+Nsqr(di+1)) * 0.5_cvmix_r8)
        
        ! calculate int_(-h)^z N dz
        phin(di)    = phin_dim1 + sqrt(Nzw(di))*dzw(di)
        phin_dim1   = phin(di)
        
        ! calculate phi_n    =    cos( int_(-h)^z N/c_n dz )*N^0.5
        !    and   dphi_n/dz =    sin( int_(-h)^z N/c_n dz )/N^0.5
        fxa         = phin(di)/cn
        dphindz(di) = sin(fxa)/Nzw(di)**0.25_cvmix_r8 
        phin(di)    = cos(fxa)*Nzw(di)**0.25_cvmix_r8
    enddo
    
    
    ! normalisation with int_(-h)^0 dz (dphi_n/dz )^2 /N^2  = 1
    norm = 0.0_cvmix_r8
    do di = 1,nlev
        norm = norm + dphindz(di)**2/Nzw(di)*dzw(di)
    enddo
    do di = 1,nlev   
        if (norm>0.0_cvmix_r8) dphindz(di) = dphindz(di)/norm**0.5
    enddo

    ! normalisation with int_(-h)^0 dz phi_n^2 /c_n^2  = 1
    norm=0.0_cvmix_r8
    do di = 1,nlev
        norm = norm + phin(di)**2/(cn**2)*dzw(di)
    enddo
    do di = 1,nlev   
        if (norm>0.0_cvmix_r8) phin(di) = phin(di)/norm**0.5
    enddo
    
    if (present(E_struct_M2)) then
        ! calculate structure function for energy: 
        ! E(z) = E_0 0.5( (1+f^2/om^2) phi_n^2/c_n^2 + (1-f^2/om^2) (dphi_n/dz)^2/N^2)
        do di = 1,nlev   
            E_struct_M2(di) = 0.5_cvmix_r8*( (1.0_cvmix_r8+coriolis**2/omega_M2**2)*phin(   di)**2 /(cn**2) &
                                            +(1.0_cvmix_r8-coriolis**2/omega_M2**2)*dphindz(di)**2 /Nzw(di) )
        enddo
    endif

    if (present(E_struct_niw)) then
        do di = 1,nlev   
            E_struct_niw(di) = 0.5_cvmix_r8*( (1.0_cvmix_r8+coriolis**2/omega_niw**2)*phin(   di)**2/(cn**2) &
                                             +(1.0_cvmix_r8-coriolis**2/omega_niw**2)*dphindz(di)**2/Nzw(di) )
        enddo
    endif
    
end subroutine compute_vert_struct_fct


!
!
!_______________________________________________________________________________
! Integrate vertical part of internal wave energy idemix equation, solve vertical 
! diffusion an dissipation implicitly over the water column. Keep in mind Eiw lives 
! on the full depth levels
subroutine compute_vdiff_vdiss_Eiw( &
                nlev             , & 
                dzw              , &
                dt               , &
                c0               , &
                alpha_c          , &
                fsrf             , &
                fbot             , &
                Eiw_old          , &
                Eiw_new          , &
                Eiw_dt           , &
                Eiw_vdif         , &
                Eiw_diss         , &
                Eiw_srf          , &
                Eiw_bot          , &
                idemix2_const_userdef &
                )
    !___Input___________________________________________________________________
    type(idemix2_type), intent(in) , optional, target :: idemix2_const_userdef
    integer           , intent(in)                    :: nlev ! --> number of full depth levels          
    real(cvmix_r8)    , intent(in), dimension(nlev-1) :: dzw ! --> layer thickness        
    real(cvmix_r8)    , intent(in)                    :: dt
    real(cvmix_r8)    , intent(in), dimension(nlev)   :: c0
    real(cvmix_r8)    , intent(in), dimension(nlev)   :: alpha_c
    real(cvmix_r8)    , intent(in)                    :: fsrf
    real(cvmix_r8)    , intent(in), dimension(nlev)   :: fbot
!     real(cvmix_r8)    , intent(in), dimension(nlev)   :: fadd --> in moment no additional forcings
    real(cvmix_r8)    , intent(in), dimension(nlev)   :: Eiw_old
    !___Output_________________________________________________________________
    real(cvmix_r8)    , intent(out), dimension(nlev)  :: Eiw_new
    real(cvmix_r8)    , intent(out), dimension(nlev)  :: Eiw_dt
    real(cvmix_r8)    , intent(out), dimension(nlev)  :: Eiw_vdif
    real(cvmix_r8)    , intent(out), dimension(nlev)  :: Eiw_diss
    real(cvmix_r8)    , intent(out), dimension(nlev)  :: Eiw_bot
    real(cvmix_r8)    , intent(out)                   :: Eiw_srf
    !___Local___________________________________________________________________
    integer                                           :: nz 
    real(cvmix_r8), dimension(nlev)                   :: Eiw_max
    real(cvmix_r8), dimension(nlev)                   :: fadd
    real(cvmix_r8), dimension(nlev-1)                 :: c0_zmid
    real(cvmix_r8), dimension(nlev)                   :: dzmid
    real(cvmix_r8), dimension(nlev)                   :: a_dif, b_dif, c_dif
    real(cvmix_r8), dimension(nlev)                   :: a_tri, b_tri, c_tri, d_tri
    type(idemix2_type), pointer                       :: idemix2_const_in     
    
    ! do pointer into save variable or into user defined input variable 
    idemix2_const_in => idemix2_constants_saved
    if (present(idemix2_const_userdef)) then
        idemix2_const_in => idemix2_const_userdef
    end if
    !___________________________________________________________________________
    a_dif       = 0.0_cvmix_r8
    b_dif       = 0.0_cvmix_r8
    c_dif       = 0.0_cvmix_r8
    a_tri       = 0.0_cvmix_r8
    b_tri       = 0.0_cvmix_r8
    c_tri       = 0.0_cvmix_r8
    d_tri       = 0.0_cvmix_r8
    Eiw_new     = 0.0_cvmix_r8
    Eiw_diss    = 0.0_cvmix_r8
    Eiw_bot     = 0.0_cvmix_r8
    
    !___________________________________________________________________________
    ! initialise additional forcings with zeros, there can be additional forcing terms
    ! from EKE, GM disspation, lateral friction, cabbeling heat, bottom friction
    ! This is not yet implemented
    fadd(:) = 0.0_cvmix_r8
    
    !___________________________________________________________________________
    ! prevent negative dissipation of IW energy, might not be neccessary but 
    ! safety first
    Eiw_max = max(0.0_cvmix_r8, Eiw_old)
    
    !___________________________________________________________________________
    ! Solve vertical part implicitly, having Eiw on mid-depth levels like tracer 
    ! for easier advection
    ! 
    ! dEiw/dt  = d/dz(c0*tau_v * d/dz(c0*Eiw)) - alpha_c*Eiw + forcing
    !            |___________________________|         |           
    !             vertical diffusion of E_iw       dissipation           
    !                
    !   ~~~~~~~o zlev=1~~~~~~~
    !          |
    !   zmid=1 +   --> k-0.5
    !          |
    !          o 2 --> k  --> Eiw_k, c0_k
    !          |
    !        2 +   --> k+0.5
    !          |
    !   -------o 3-----------
    !   /////////////////////       
    !  c0_k-0.5 = (c0_k+c0_k-1) / 2
    !  c0_k+0.5 = (c0_k+c0_k+1) / 2
    ! 
    !  d/dz(c0*Eiw)_k = (c0_k*Eiw_k - c0_k-1*Eiw_k-1)/h(k-1) = dcE_k-0.5--> defined on middepth levels
    !  d/dz(c0* dcE_k-0.5) = (c0_k+0.5*dcE_k+0.5 - c0_k-0.5*dcE_k-0.5) * 2 / (h_k+ h_k-1)
    ! 
    ! Discrite Form: 
    ! --------------
    ! Eiw^(t+1) - Eiw^(t) = dz* [d/dz(c0*tau_v * d/dz(c0*Eiw)) - alpha_c*Eiw]^(t+1) + forcing^(t)
    ! Eiw^(t+1)_k = Eiw^(t)_k
    !               + dt * forcing^(t)_k
    !               - dt * alpha_c * Eiw^(t+1)_k
    !               - tau_v * [ c0_k+0.5 * (c0_k+1*Eiw_k+1 - c0_k  *Eiw_k  ) * 2 / h_k   / (h_k+h_k-1)]
    !                          -c0_k-0.5 * (c0_k  *Eiw_k   - c0_k-1*Eiw_k-1) * 2 / h_k-1 / (h_k+h_k-1)]
    ! Rearange for Tridiagonal:
    ! -------------------------
    ! a_k * Eiw^(t+1)_k-1 + b_k * Eiw^(t+1)_k + c_k * Eiw^(t+1) = d_k
    !
    ! Estimate coefficents a,b,c:
    ! ---------------------------
    ! 1) a_k --> (k-1, t+1): 
    !   - dt * 2 / (h_k-1*(h_k+h_k-1)) * c0_(k-0.5)*c0_k-1*Eiw^(t+1)_k-1  = a_k * Eiw^(t+1)_k-1
    !
    ! 2) b_k --> (k  , t+1): 
    !   Eiw^(t+1)_k 
    !   + dt*2/(h_k  *(h_k+h_k-1))*c0_(k+0.5)*c0_k*Eiw^(t+1)_k
    !   + dt*2/(h_k-1*(h_k+h_k-1))*c0_(k-0.5)*c0_k*Eiw^(t+1)_k 
    !   + dt*alpha_c*Eiw^(t+1)_k
    !   = b_k * Eiw^(t+1)_k
    !
    ! 3) c_k --> (k+1  , t+1): 
    !   - dt * 2 / (h_k*(h_k+h_k-1)) * c0_(k+0.5)*c0_k+1*Eiw^(t+1)_k+1  = c_k * Eiw^(t+1)_k+1
    !
    ! 4) d_k --> (k    , t  ): 
    !     Eiw^t_k + dt* Forc^t_k  = d_k
    !
    ! Solve tridiagonal matrix: 
    ! -------------------------
    ! | b_1 c_1   0   0   0 | (Eiw_1) = (d_1)
    ! | a_2 b_2 c_2   0   0 | (Eiw_2) = (d_2)
    ! |   0 a_3 b_3 c_3   0 | (Eiw_3) = (d_3)
    ! |   0   0 a_4 b_4 c_4 | (Eiw_4) = (d_4)
    ! |   0   0   0 a_k b_k | (Eiw_k) = (d_k)
    ! 
    !___________________________________________________________________________ 
    do nz = 2, nlev-1
        c0_zmid(nz-1) = dt * idemix2_const_in%tau_v * (c0(nz)+c0(nz-1))*0.5_cvmix_r8 / dzw(nz-1)
        dzmid(nz)     = (dzw(nz) + dzw(nz-1))*0.5_cvmix_r8
    end do
    dzmid(1)      = dzw(1)*0.5_cvmix_r8
    nz=nlev
    dzmid(nz)     = dzw(nz-1)*0.5_cvmix_r8
    c0_zmid(nz-1) = idemix2_const_in%tau_v * (c0(nz)+c0(nz-1))*0.5_cvmix_r8 / dzw(nz-1)
    
    !___________________________________________________________________________
    ! -- a -- 
    a_dif(1) = 0.0_cvmix_r8 ! not part of the diffusion matrix, thus value is arbitrary
    do nz=2,nlev
        a_dif(nz) = c0_zmid(nz-1)*c0(nz-1) / dzmid(nz) 
    end do
    
    ! -- b -- 
    nz = 1
    b_dif(nz) = c0_zmid(nz  )*c0(nz) / dzmid(nz)
    do nz=2,nlev-1
        b_dif(nz) = (c0_zmid(nz-1)*c0(nz)+c0_zmid(nz)*c0(nz)) / dzmid(nz)
    end do
    nz = nlev
    b_dif(nz) = c0_zmid(nz-1)*c0(nz) / dzmid(nz)
    
    ! -- c-- 
    do nz=1,nlev-1
        c_dif(nz) = c0_zmid(nz)*c0(nz+1) / dzmid(nz) 
    end do
    c_dif(nlev) = 0.0_cvmix_r8
    
    !___________________________________________________________________________
    ! prepare tridiagonal matrix elements: a,b,c,d
    nz = 1
    a_tri(nz)     =  -dt * a_dif(nz)
    b_tri(nz)     = 1.0_cvmix_r8                    &
                    + dt * b_dif(nz)
    c_tri(nz)     =  -dt * c_dif(nz)
    do nz=2,nlev-1
        a_tri(nz) =  -dt * a_dif(nz)
        b_tri(nz) = 1.0_cvmix_r8                    &
                    + dt * b_dif(nz)                &
                    + dt * alpha_c(nz)*Eiw_max(nz)
        c_tri(nz) =  -dt * c_dif(nz)
    end do
    nz = nlev
    a_tri(nz)     =  -dt * a_dif(nz)
    b_tri(nz)     = 1.0_cvmix_r8                    &
                    + dt * b_dif(nz)
    c_tri(nz)     =  -dt * c_dif(nz)
    
    ! -- d ---
    do nz=1,nlev
        ! bottom production term
        Eiw_bot(nz) = fbot(nz)/dzmid(nz)
        ! --> fbot(nz) must be here an array overf depth since we want this 
        !     in FESOM2 to operate on nodes, but our bottom topography is on elements
        !     so bottom forcing is injected over the water column where elements 
        !     that participate to the scalar control volume turn into bottom over the 
        !     vertical depth
        
        d_tri(nz) = Eiw_old(nz) + dt*fadd(nz) + dt*Eiw_bot(nz)
        
    end do
    d_tri(1)   = d_tri(1) + dt*fsrf/dzmid(1)
    
    !___________________________________________________________________________
    ! solve tridiagonal matrix
    call solve_tridiag(a_tri, b_tri, c_tri, d_tri, Eiw_new, nlev)
    
    !___________________________________________________________________________
    ! diagnose implicite tendencies (only for diagnostics)
    ! production of Eiw from vertical diffusion
    nz = 1
    Eiw_vdif(nz) = - b_dif(nz)*Eiw_new(nz) + c_dif(nz)*Eiw_new(nz+1)
    do nz=2,nlev-1
        Eiw_vdif(nz) = a_dif(nz)*Eiw_new(nz-1) - b_dif(nz)*Eiw_new(nz) + c_dif(nz)*Eiw_new(nz+1)
    enddo
    nz = nlev
    Eiw_vdif(nz) = a_dif(nz)*Eiw_new(nz-1) - b_dif(nz)*Eiw_new(nz)
  
    !___________________________________________________________________________
    ! dissipation of E_iw
    Eiw_diss(:) = -alpha_c(: )*Eiw_max(:)*Eiw_new(:)
!     write(*,*) "min/max Eiw_diss: ", minval(Eiw_diss), maxval(Eiw_diss)
!     write(*,*) "min/max alpha_c : ", minval(alpha_c) , maxval(alpha_c)
!     write(*,*) "min/max Eiw_max : ", minval(Eiw_max) , maxval(Eiw_max)
!     write(*,*) "min/max Eiw_new : ", minval(Eiw_new) , maxval(Eiw_new)
!     write(*,*) "min/max Eiw_old : ", minval(Eiw_old) , maxval(Eiw_old)
!     write(*,*)
    
    ! total production of Eiw
    Eiw_dt(:)   = (Eiw_new(:) - Eiw_old(:))/dt
    
    ! surface production term
    Eiw_srf     =  fsrf/dzmid(1)
    
end subroutine compute_vdiff_vdiss_Eiw


!
!
!_______________________________________________________________________________
! Integrate vertical part of internal wave energy idemix equation, solve vertical 
! diffusion an dissipation implicitly over the water column. Keep in mind Eiw lives 
! on the full depth levels
subroutine compute_Eiw_waveinteract( &
                nlev               , & 
                nfbin              , &
                dzw                , &
                dphi               , &
                dt                 , &
                E_iw_old           , &
                E_iw_new           , &
                E_M2_old           , &
                E_M2_new           , &
                E_M2_struct        , &
                alpha_M2_c         , &
                tau_M2             , &
                E_niw_old          , &
                E_niw_new          , &
                E_niw_struct       , &
                tau_niw            , &
                idemix2_const_userdef &
                )
    !___Input/Output____________________________________________________________
    type(idemix2_type), intent(in   ), optional, target :: idemix2_const_userdef
    integer           , intent(in   )                   :: nlev  ! --> number of full depth levels
    integer           , intent(in   )                   :: nfbin ! --> of spectral bins 
    real(cvmix_r8)    , intent(in   )                   :: dzw(nlev-1)   ! --> layer thickness
    real(cvmix_r8)    , intent(in   )                   :: dphi(nfbin) 
    real(cvmix_r8)    , intent(in   )                   :: dt
    real(cvmix_r8)    , intent(in   )                   :: E_iw_old(nlev)
    real(cvmix_r8)    , intent(inout)                   :: E_iw_new(nlev)
    real(cvmix_r8)    , intent(in   ), optional         :: E_M2_old(nfbin)
    real(cvmix_r8)    , intent(inout), optional         :: E_M2_new(nfbin)
    real(cvmix_r8)    , intent(in   ), optional         :: E_M2_struct(nfbin)
    real(cvmix_r8)    , intent(in   ), optional         :: alpha_M2_c
    real(cvmix_r8)    , intent(in   ), optional         :: tau_M2
    real(cvmix_r8)    , intent(in   ), optional         :: E_niw_old(nfbin)
    real(cvmix_r8)    , intent(inout), optional         :: E_niw_new(nfbin)
    real(cvmix_r8)    , intent(in   ), optional         :: E_niw_struct(nfbin)
    real(cvmix_r8)    , optional                        :: tau_niw
    
    !___Local___________________________________________________________________
    integer                                             :: nz, fbini 
    real(cvmix_r8)                                      :: vint, sint, aM2c, M2_diss, fmin 
    type(idemix2_type), pointer                         :: idemix2_const_in     
    
    ! do pointer into save variable or into user defined input variable 
    idemix2_const_in => idemix2_constants_saved
    if (present(idemix2_const_userdef)) then
        idemix2_const_in => idemix2_const_userdef
    end if
    
    !___________________________________________________________________________
    ! vertical integrate E_iw --> compute dzt from dzw
    vint = 0.0_cvmix_r8
    nz   = 1
    vint = vint + E_iw_old(nz)*dzw(nz)*0.5_cvmix_r8
    do nz = 2, nlev-1
        vint = vint + E_iw_old(nz)*(dzw(nz)+dzw(nz-1))*0.5_cvmix_r8
    end do
    nz   = nlev
    vint = vint + E_iw_old(nz)*dzw(nz-1)*0.5_cvmix_r8
    
    !___________________________________________________________________________
    if (present(E_M2_old)) then 
        ! spectrally integrate E_M2
        sint = 0.0_cvmix_r8
        do fbini = 2, nfbin-1
            sint = sint + E_M2_old(fbini)*dphi(fbini)
        end do
        
        ! update M2 energy: interaction of M2 and continuum
        aM2c = alpha_M2_c
        aM2c = max(0.0_cvmix_r8,min(aM2c,1./max(1d-12,dt*vint))) ! FP 2020
        aM2c = max(0.0_cvmix_r8,min(aM2c,1./max(1d-12,dt*sint)))
        do fbini = 2, nfbin-1
            M2_diss = aM2c*vint*E_M2_old(fbini)
            E_M2_new(fbini)= E_M2_new(fbini)-dt*M2_diss 
            
            ! update E_iw from E_M2
            fmin = min( 0.5_cvmix_r8/dt,aM2c*vint ) ! flux limiter
            E_iw_new(fbini) =   E_iw_new(fbini) &
                              + dt * tau_M2* sint * E_M2_struct(fbini) &
                              + dt * fmin  * sint * E_M2_struct(fbini)
            
        end do
    end if
    
    !___________________________________________________________________________
    if (present(E_niw_old)) then 
        ! spectrally integrate E_niw
        sint = 0.0_cvmix_r8
        do fbini = 2, nfbin-1
            sint = sint + E_niw_old(fbini)*dphi(fbini)
        end do
        
        ! update E_iw from E_niw
        do fbini = 2, nfbin-1
            E_iw_new(fbini) =   E_iw_new(fbini) &
                              + dt * tau_niw * sint * E_niw_struct(fbini)
        end do
    end if
    
    
end subroutine compute_Eiw_waveinteract



!
!
!_______________________________________________________________________________
! IDEMIX2 only shortcut: derive diffusivity and viscosity using Osbourne relation
subroutine compute_Eiw_diss2KvAv(nlev, Eiw_diss, Nsqr, KappaH, KappaM)
    !___Input/Output____________________________________________________________
    integer           , intent(in   )                   :: nlev  ! --> number of full depth levels
    real(cvmix_r8)    , intent(in   )                   :: Eiw_diss(nlev)   ! --> layer thickness
    real(cvmix_r8)    , intent(in   )                   :: Nsqr(   nlev)
    real(cvmix_r8)    , intent(inout)                   :: KappaH(  nlev)
    real(cvmix_r8)    , intent(inout)                   :: KappaM(  nlev)

    !___Local___________________________________________________________________
    integer                                             :: nz
    
    KappaH = 0.0_cvmix_r8
    KappaM = 0.0_cvmix_r8
    do nz = 2, nlev-1
        KappaH(nz) =  0.2_cvmix_r8/(1.0_cvmix_r8+0.2_cvmix_r8) * (-1.0_cvmix_r8 * Eiw_diss(nz)) / max(1.0e-12_cvmix_r8, Nsqr(nz))
        KappaH(nz) = max(1.0e-9_cvmix_r8, KappaH(nz))
        KappaH(nz) = min(1.0_cvmix_r8, KappaH(nz))
        KappaM(nz) =  10.0_cvmix_r8 * KappaH(nz)
    end do
end subroutine compute_Eiw_diss2KvAv


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
        case('shelf_dist'  ) ; idemix2_constants_out%shelf_dist   = val
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
        case('hor_diff_niter'     ) ; idemix2_constants_out%hor_diff_niter= val
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
        case('enable_hor_diff_iter') ; idemix2_constants_out%enable_hor_diff_iter= val
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