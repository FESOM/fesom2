MODULE MOD_ICE
USE o_PARAM, only: WP
USE, intrinsic :: ISO_FORTRAN_ENV
USE MOD_WRITE_BINARY_ARRAYS
USE MOD_READ_BINARY_ARRAYS
IMPLICIT NONE
SAVE

!
! 
!_______________________________________________________________________________
! set data array derived type for ice-tracers (area, mice, msnow) more tracer 
! are theretical possible
TYPE T_ICE_DATA
    !___________________________________________________________________________
    real(kind=WP), allocatable, dimension(:)    :: values, values_old, values_rhs, & 
                                                   values_div_rhs, dvalues, valuesl
    integer                                     :: ID
    !___________________________________________________________________________
    contains
#if defined(__PGI)
    private
#endif            
        procedure WRITE_T_ICE_DATA
        procedure READ_T_ICE_DATA
        generic :: write(unformatted) => WRITE_T_ICE_DATA
        generic :: read(unformatted)  => READ_T_ICE_DATA
END TYPE T_ICE_DATA
!
! 
!_______________________________________________________________________________
! set work array derived type for ice
TYPE T_ICE_WORK
    !___________________________________________________________________________
    real(kind=WP), allocatable, dimension(:)    :: fct_tmax, fct_tmin
    real(kind=WP), allocatable, dimension(:)    :: fct_plus, fct_minus
    real(kind=WP), allocatable, dimension(:,:)  :: fct_fluxes
    real(kind=WP), allocatable, dimension(:)    :: fct_massmatrix
    real(kind=WP), allocatable, dimension(:)    :: sigma11, sigma12, sigma22
    real(kind=WP), allocatable, dimension(:)    :: eps11, eps12, eps22
    real(kind=WP), allocatable, dimension(:)    :: ice_strength, inv_areamass, inv_mass
    !___________________________________________________________________________
    contains
#if defined(__PGI)
    private
#endif            
        procedure WRITE_T_ICE_WORK
        procedure READ_T_ICE_WORK
        generic :: write(unformatted) => WRITE_T_ICE_WORK
        generic :: read(unformatted)  => READ_T_ICE_WORK
END TYPE T_ICE_WORK
!
! 
!_______________________________________________________________________________
! set work array derived type for ice
TYPE T_ICE_THERMO
    !___________________________________________________________________________
    real(kind=WP), allocatable, dimension(:)    :: t_skin, thdgr, thdgrsn, thdgr_old, ustar
    !___________________________________________________________________________
    real(kind=WP) :: rhoair=1.3  , inv_rhoair=1./1.3  ! Air density & inverse ,  LY2004 !1.3 AOMIP
    real(kind=WP) :: rhowat=1025., inv_rhowat=1./1025.! Water density & inverse
    real(kind=WP) :: rhofwt=1000., inv_rhofwt=1./1000.! Freshwter density & inverse
    real(kind=WP) :: rhoice=910. , inv_rhoice=1./910. ! Ice density & inverse, AOMIP
    real(kind=WP) :: rhosno=290. , inv_rhosno=1./290. ! Snow density & inverse, AOMIP
    ! Specific heat of air, ice, snow [J/(kg * K)] 
    real(kind=WP) :: cpair=1005., cpice=2106., cpsno=2090. 
!     real(kind=WP) :: cc=rhowat*4190.0  ! Volumetr. heat cap. of water [J/m**3/K](cc = rhowat*cp_water)
!     real(kind=WP) :: cl=rhoice*3.34e5  ! Volumetr. latent heat of ice fusion [J/m**3](cl=rhoice*Lf)
! --> cl and cc are setted in subroutine ice_init(...)
    real(kind=WP) :: cc=1025.*4190.0  ! Volumetr. heat cap. of water [J/m**3/K](cc = rhowat*cp_water)
    real(kind=WP) :: cl=910.*3.34e5  ! Volumetr. latent heat of ice fusion [J/m**3](cl=rhoice*Lf) 
    real(kind=WP) :: clhw=2.501e6      ! Specific latent heat [J/kg]: water	-> water vapor
    real(kind=WP) :: clhi=2.835e6      !                              sea ice-> water vapor
    real(kind=WP) :: tmelt=273.15      ! 0 deg C expressed in K 
    real(kind=WP) :: boltzmann=5.67E-8 ! S. Boltzmann const.*longw. emissivity
    integer       :: iclasses=7        ! Number of ice thickness gradations for ice growth calcs.
    real(kind=WP) :: hmin= 0.01        ! Cut-off ice thickness     !!
    real(kind=WP) :: Armin=0.01        ! Minimum ice concentration !!
    
    ! --- namelist parameter /ice_therm/
    real(kind=WP) :: con= 2.1656, consn = 0.31 ! Thermal conductivities: ice & snow; W/m/K
    real(kind=WP) :: Sice = 4.0        ! Ice salinity 3.2--5.0 ppt.
    real(kind=WP) :: h0=1.0	           ! Lead closing parameter [m] ! 0.5
    real(kind=WP) :: emiss_ice=0.97    ! Emissivity of Snow/Ice, 
    real(kind=WP) :: emiss_wat=0.97    ! Emissivity of open water
    real(kind=WP) :: albsn = 0.81      ! Albedo: frozen snow
    real(kind=WP) :: albsnm= 0.77      !         melting snow
    real(kind=WP) :: albi  = 0.70      !         frozen ice
    real(kind=WP) :: albim = 0.68      !         melting ice
    real(kind=WP) :: albw  = 0.066     !         open water, LY2004
    contains
#if defined(__PGI)
    private
#endif            
        procedure WRITE_T_ICE_THERMO
        procedure READ_T_ICE_THERMO
        generic :: write(unformatted) => WRITE_T_ICE_THERMO
        generic :: read(unformatted)  => READ_T_ICE_THERMO
END TYPE T_ICE_THERMO
!
! 
!_______________________________________________________________________________
! set work array derived type for ice
#if defined (__oasis) || defined (__ifsinterface)
TYPE T_ICE_ATMCOUPL

    !___________________________________________________________________________
    real(kind=WP), allocatable, dimension(:)    :: oce_flx_h, ice_flx_h, tmpoce_flx_h, tmpice_flx_h
#if defined (__oifs) || defined (__ifsinterface)
    !___________________________________________________________________________
    real(kind=WP), allocatable, dimension(:)    :: ice_alb, enthalpyoffuse
    ! !!! DONT FORGET ice_temp rhs_tempdiv rhs_temp is advected for oifs !!! --> becomes additional ice 
    ! tracer in ice%data(4)%values
#endif /* (__oifs)  */
    !___________________________________________________________________________
    contains
#if defined(__PGI)
    private
#endif            
        procedure WRITE_T_ICE_ATMCOUPL
        procedure READ_T_ICE_ATMCOUPL
        generic :: write(unformatted) => WRITE_T_ICE_ATMCOUPL
        generic :: read(unformatted)  => READ_T_ICE_ATMCOUPL
END TYPE T_ICE_ATMCOUPL
#endif /* (__oasis) */

!
! 
!_______________________________________________________________________________
! set main ice derived type contains parameters, data array, work array, u_ice, vice
TYPE T_ICE

    !___________________________________________________________________________
    ! zonal & merdional ice velocity
    real(kind=WP), allocatable, dimension(:)    :: uice, uice_rhs, uice_old, uice_aux
    real(kind=WP), allocatable, dimension(:)    :: vice, vice_rhs, vice_old, vice_aux
    
    ! surface stess atm<-->ice, oce<-->ice
    real(kind=WP), allocatable, dimension(:)    :: stress_atmice_x, stress_iceoce_x
    real(kind=WP), allocatable, dimension(:)    :: stress_atmice_y, stress_iceoce_y
    
    ! oce temp, salt, ssh, and uv at surface
    real(kind=WP), allocatable, dimension(:)    :: srfoce_temp, srfoce_salt, srfoce_ssh
!     real(kind=WP), allocatable, dimension(:,:)  :: srfoce_uv
    real(kind=WP), allocatable, dimension(:)    :: srfoce_u, srfoce_v
    
    ! freshwater & heatflux
    real(kind=WP), allocatable, dimension(:)    :: flx_fw, flx_h
    
    ! maEVP variables
    real(kind=WP), allocatable, dimension(:)    :: alpha_evp_array, beta_evp_array
    
    !___________________________________________________________________________
    ! total number of ice tracers (default=3, 1=area, 2=mice, 3=msnow, (4=ice_temp)
#if defined (__oifs) || defined (__ifsinterface)
    integer                                     :: num_itracers=4
#else
    integer                                     :: num_itracers=3
#endif 

    ! put ice tracers data arrays
    type(t_ice_data), allocatable, dimension(:) :: data
    
    !___________________________________________________________________________
    ! put ice working arrays
    type(t_ice_work)                            :: work
    
    ! put thermodynamics arrays
    type(t_ice_thermo)                          :: thermo
    
#if defined (__oasis) || defined (__ifsinterface)
    !___________________________________________________________________________
    ! put ice arrays for coupled model
    type(t_ice_atmcoupl)                        :: atmcoupl
#endif /* (__oasis) */ 

    !___________________________________________________________________________
    ! set ice model parameters:
    ! --- RHEOLOGY ---
    real(kind=WP)             :: pstar      = 30000.0_WP   ![N/m^2]
    real(kind=WP)             :: ellipse    = 2.0_WP       !
    real(kind=WP)             :: c_pressure = 20.0_WP      !
    real(kind=WP)             :: delta_min  = 1.0e-11      ! [s^(-1)]
    real(kind=WP)             :: Clim_evp   = 615          ! kg/m^2
    real(kind=WP)             :: zeta_min   = 4.0e+8       ! kg/s
    integer                   :: evp_rheol_steps=120       ! EVP rheology cybcycling steps
    real(kind=WP)             :: ice_gamma_fct=0.25_WP     ! smoothing parameter in ice fct advection
    real(kind=WP)             :: ice_diff   = 10.0_WP      ! diffusion to stabilize ice advection
    real(kind=WP)             :: theta_io   =0.0_WP        ! rotation angle (ice-ocean), available
    ! --- in EVP ---
    real(kind=WP)             :: alpha_evp=250, beta_evp=250
    real(kind=WP)             :: c_aevp=0.15               ! 0.1--0.2, but should be adjusted experimentally
    ! --- Ice forcing averaging ---
    integer                   :: ice_ave_steps=1           !ice step=ice_ave_steps*oce_step
    real(kind=WP)             :: cd_oce_ice = 5.5e-3       ! drag coef. oce - ice      
    logical                   :: ice_free_slip=.false.
    integer                   :: whichEVP=0                ! 0=standart; 1=mEVP; 2=aEVP
    
    real(kind=WP)             :: ice_dt                    ! ice step=ice_ave_steps*oce_step
    real(kind=WP)             :: Tevp_inv   
    
    integer                   :: ice_steps_since_upd=0
    logical                   :: ice_update = .true.
    !___________________________________________________________________________
    contains
#if defined(__PGI)
    private
#endif            
        procedure WRITE_T_ICE
        procedure READ_T_ICE
        generic :: write(unformatted) => WRITE_T_ICE
        generic :: read(unformatted)  => READ_T_ICE
END TYPE T_ICE

contains 
!
!
!_______________________________________________________________________________
! Unformatted writing for T_ICE_DATA
subroutine WRITE_T_ICE_DATA(tdata, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_ICE_DATA),      intent(in)     :: tdata
    integer,                intent(in)     :: unit
    integer,                intent(out)    :: iostat
    character(*),           intent(inout)  :: iomsg
    call write_bin_array(tdata%values,         unit, iostat, iomsg)
    call write_bin_array(tdata%values_old,     unit, iostat, iomsg)
    call write_bin_array(tdata%values_rhs,     unit, iostat, iomsg)
    call write_bin_array(tdata%values_div_rhs, unit, iostat, iomsg)
    call write_bin_array(tdata%dvalues,        unit, iostat, iomsg)
    call write_bin_array(tdata%valuesl,        unit, iostat, iomsg)
    write(unit, iostat=iostat, iomsg=iomsg) tdata%ID
end subroutine WRITE_T_ICE_DATA  

! Unformatted reading for T_ICE_DATA
subroutine READ_T_ICE_DATA(tdata, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_ICE_DATA),      intent(inout)  :: tdata
    integer,                intent(in)     :: unit
    integer,                intent(out)    :: iostat
    character(*),           intent(inout)  :: iomsg
    call read_bin_array(tdata%values,         unit, iostat, iomsg)
    call read_bin_array(tdata%values_old,     unit, iostat, iomsg)
    call read_bin_array(tdata%values_rhs,     unit, iostat, iomsg)
    call read_bin_array(tdata%values_div_rhs, unit, iostat, iomsg)
    call read_bin_array(tdata%dvalues,        unit, iostat, iomsg)
    call read_bin_array(tdata%valuesl,        unit, iostat, iomsg)
    read(unit, iostat=iostat, iomsg=iomsg) tdata%ID
end subroutine READ_T_ICE_DATA                                   
!
!
!_______________________________________________________________________________
! Unformatted writing for T_ICE_WORK
subroutine WRITE_T_ICE_WORK(twork, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_ICE_WORK),      intent(in)     :: twork
    integer,                intent(in)     :: unit
    integer,                intent(out)    :: iostat
    character(*),           intent(inout)  :: iomsg
    call write_bin_array(twork%fct_tmax,     unit, iostat, iomsg)
    call write_bin_array(twork%fct_tmin,     unit, iostat, iomsg)
    call write_bin_array(twork%fct_plus,     unit, iostat, iomsg)
    call write_bin_array(twork%fct_minus,    unit, iostat, iomsg)
    call write_bin_array(twork%fct_fluxes,   unit, iostat, iomsg)
    call write_bin_array(twork%fct_massmatrix,unit, iostat, iomsg)
    call write_bin_array(twork%sigma11,      unit, iostat, iomsg)
    call write_bin_array(twork%sigma12,      unit, iostat, iomsg)
    call write_bin_array(twork%sigma22,      unit, iostat, iomsg)
    call write_bin_array(twork%eps11,        unit, iostat, iomsg)
    call write_bin_array(twork%eps12,        unit, iostat, iomsg)
    call write_bin_array(twork%eps22,        unit, iostat, iomsg)
    call write_bin_array(twork%ice_strength, unit, iostat, iomsg)
    call write_bin_array(twork%inv_areamass, unit, iostat, iomsg)
    call write_bin_array(twork%inv_mass,     unit, iostat, iomsg)
end subroutine WRITE_T_ICE_WORK    

! Unformatted reading for T_ICE_WORK
subroutine READ_T_ICE_WORK(twork, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_ICE_WORK),      intent(inout)  :: twork
    integer,                intent(in)     :: unit
    integer,                intent(out)    :: iostat
    character(*),           intent(inout)  :: iomsg
    call read_bin_array(twork%fct_tmax,     unit, iostat, iomsg)
    call read_bin_array(twork%fct_tmin,     unit, iostat, iomsg)
    call read_bin_array(twork%fct_plus,     unit, iostat, iomsg)
    call read_bin_array(twork%fct_minus,    unit, iostat, iomsg)
    call read_bin_array(twork%fct_fluxes,   unit, iostat, iomsg)
    call read_bin_array(twork%fct_massmatrix,unit, iostat, iomsg)
    call read_bin_array(twork%sigma11,      unit, iostat, iomsg)
    call read_bin_array(twork%sigma12,      unit, iostat, iomsg)
    call read_bin_array(twork%sigma22,      unit, iostat, iomsg)
    call read_bin_array(twork%eps11,        unit, iostat, iomsg)
    call read_bin_array(twork%eps12,        unit, iostat, iomsg)
    call read_bin_array(twork%eps22,        unit, iostat, iomsg)
    call read_bin_array(twork%ice_strength, unit, iostat, iomsg)
    call read_bin_array(twork%inv_areamass, unit, iostat, iomsg)
    call read_bin_array(twork%inv_mass,     unit, iostat, iomsg)
end subroutine READ_T_ICE_WORK
!
!
!_______________________________________________________________________________
! Unformatted writing for T_ICE_WORK
subroutine WRITE_T_ICE_THERMO(ttherm, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_ICE_THERMO),    intent(in)     :: ttherm
    integer,                intent(in)     :: unit
    integer,                intent(out)    :: iostat
    character(*),           intent(inout)  :: iomsg
    call write_bin_array(ttherm%t_skin,       unit, iostat, iomsg)
    call write_bin_array(ttherm%thdgr,        unit, iostat, iomsg)
    call write_bin_array(ttherm%thdgrsn,      unit, iostat, iomsg)
    call write_bin_array(ttherm%thdgr_old,    unit, iostat, iomsg)
    call write_bin_array(ttherm%ustar,        unit, iostat, iomsg)
end subroutine WRITE_T_ICE_THERMO    

! Unformatted reading for T_ICE_WORK
subroutine READ_T_ICE_THERMO(ttherm, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_ICE_THERMO),    intent(inout)  :: ttherm
    integer,                intent(in)     :: unit
    integer,                intent(out)    :: iostat
    character(*),           intent(inout)  :: iomsg
    call read_bin_array(ttherm%t_skin,       unit, iostat, iomsg)
    call read_bin_array(ttherm%thdgr,        unit, iostat, iomsg)
    call read_bin_array(ttherm%thdgrsn,      unit, iostat, iomsg)
    call read_bin_array(ttherm%thdgr_old,    unit, iostat, iomsg)
    call read_bin_array(ttherm%ustar,        unit, iostat, iomsg)
end subroutine READ_T_ICE_THERMO
!
!
!_______________________________________________________________________________
! Unformatted writing for T_ICE_ATMCOUPL
#if defined (__oasis) || defined (__ifsinterface)    
subroutine WRITE_T_ICE_ATMCOUPL(tcoupl, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_ICE_ATMCOUPL),  intent(in)     :: tcoupl
    integer,                intent(in)     :: unit
    integer,                intent(out)    :: iostat
    character(*),           intent(inout)  :: iomsg
    call write_bin_array(tcoupl%oce_flx_h,      unit, iostat, iomsg)
    call write_bin_array(tcoupl%ice_flx_h,      unit, iostat, iomsg)
    call write_bin_array(tcoupl%tmpoce_flx_h,   unit, iostat, iomsg)
    call write_bin_array(tcoupl%tmpice_flx_h,   unit, iostat, iomsg)
#if defined (__oifs) || defined (__ifsinterface)  
    call write_bin_array(tcoupl%ice_alb,        unit, iostat, iomsg)
    call write_bin_array(tcoupl%enthalpyoffuse, unit, iostat, iomsg)
#endif /* (__oifs) */
end subroutine WRITE_T_ICE_ATMCOUPL  
#endif /* (__oasis) */    

! Unformatted reading for T_ICE_ATMCOUPL
#if defined (__oasis) || defined (__ifsinterface)
subroutine READ_T_ICE_ATMCOUPL(tcoupl, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_ICE_ATMCOUPL),  intent(inout)  :: tcoupl
    integer,                intent(in)     :: unit
    integer,                intent(out)    :: iostat
    character(*),           intent(inout)  :: iomsg
    call read_bin_array(tcoupl%oce_flx_h, unit, iostat, iomsg)
    call read_bin_array(tcoupl%ice_flx_h, unit, iostat, iomsg)
    call read_bin_array(tcoupl%tmpoce_flx_h, unit, iostat, iomsg)
    call read_bin_array(tcoupl%tmpice_flx_h, unit, iostat, iomsg)
#if defined (__oifs) || defined (__ifsinterface)  
    call read_bin_array(tcoupl%ice_alb, unit, iostat, iomsg)
    call read_bin_array(tcoupl%enthalpyoffuse, unit, iostat, iomsg)
#endif /* (__oifs) */
end subroutine READ_T_ICE_ATMCOUPL
#endif /* (__oasis) */   
!
!
!_______________________________________________________________________________
! Unformatted writing for T_ICE_ATMCOUPL
subroutine WRITE_T_ICE(ice, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_ICE),           intent(in)     :: ice
    integer,                intent(in)     :: unit
    integer,                intent(out)    :: iostat
    character(*),           intent(inout)  :: iomsg
    integer                                :: i
    !___________________________________________________________________________
    write(unit, iostat=iostat, iomsg=iomsg) ice%num_itracers
    do i=1, ice%num_itracers
       write(unit, iostat=iostat, iomsg=iomsg) ice%data(i)
    end do
    !___________________________________________________________________________
    write(unit, iostat=iostat, iomsg=iomsg) ice%thermo
    write(unit, iostat=iostat, iomsg=iomsg) ice%work
#if defined (__oasis) || defined (__ifsinterface)
    write(unit, iostat=iostat, iomsg=iomsg) ice%atmcoupl
#endif /* (__oasis) */       

    !___________________________________________________________________________
    write(unit, iostat=iostat, iomsg=iomsg) ice%pstar
    write(unit, iostat=iostat, iomsg=iomsg) ice%ellipse
    write(unit, iostat=iostat, iomsg=iomsg) ice%c_pressure
    write(unit, iostat=iostat, iomsg=iomsg) ice%delta_min
    write(unit, iostat=iostat, iomsg=iomsg) ice%Clim_evp
    write(unit, iostat=iostat, iomsg=iomsg) ice%zeta_min
    write(unit, iostat=iostat, iomsg=iomsg) ice%evp_rheol_steps
    write(unit, iostat=iostat, iomsg=iomsg) ice%ice_gamma_fct
    write(unit, iostat=iostat, iomsg=iomsg) ice%ice_diff
    write(unit, iostat=iostat, iomsg=iomsg) ice%Tevp_inv
    write(unit, iostat=iostat, iomsg=iomsg) ice%theta_io
    write(unit, iostat=iostat, iomsg=iomsg) ice%alpha_evp
    write(unit, iostat=iostat, iomsg=iomsg) ice%beta_evp
    write(unit, iostat=iostat, iomsg=iomsg) ice%c_aevp
    write(unit, iostat=iostat, iomsg=iomsg) ice%ice_ave_steps
    write(unit, iostat=iostat, iomsg=iomsg) ice%cd_oce_ice
    write(unit, iostat=iostat, iomsg=iomsg) ice%ice_free_slip
    write(unit, iostat=iostat, iomsg=iomsg) ice%whichEVP
    write(unit, iostat=iostat, iomsg=iomsg) ice%ice_dt
    write(unit, iostat=iostat, iomsg=iomsg) ice%Tevp_inv
    write(unit, iostat=iostat, iomsg=iomsg) ice%ice_steps_since_upd
    write(unit, iostat=iostat, iomsg=iomsg) ice%ice_update
    
    !___________________________________________________________________________
    call write_bin_array(ice%uice           , unit, iostat, iomsg)
    call write_bin_array(ice%uice_rhs       , unit, iostat, iomsg)
    call write_bin_array(ice%uice_old       , unit, iostat, iomsg)
    if (ice%whichEVP /= 0) call write_bin_array(ice%uice_aux        , unit, iostat, iomsg)
    call write_bin_array(ice%vice           , unit, iostat, iomsg)
    call write_bin_array(ice%vice_rhs       , unit, iostat, iomsg)
    call write_bin_array(ice%vice_old       , unit, iostat, iomsg)
    if (ice%whichEVP /= 0) call write_bin_array(ice%vice_aux        , unit, iostat, iomsg)
    call write_bin_array(ice%stress_atmice_x, unit, iostat, iomsg)
    call write_bin_array(ice%stress_iceoce_x, unit, iostat, iomsg)
    call write_bin_array(ice%stress_atmice_y, unit, iostat, iomsg)
    call write_bin_array(ice%stress_iceoce_y, unit, iostat, iomsg)
    call write_bin_array(ice%srfoce_u       , unit, iostat, iomsg)
    call write_bin_array(ice%srfoce_v       , unit, iostat, iomsg)
    call write_bin_array(ice%srfoce_temp    , unit, iostat, iomsg)
    call write_bin_array(ice%srfoce_salt    , unit, iostat, iomsg)
    call write_bin_array(ice%srfoce_ssh     , unit, iostat, iomsg)
    call write_bin_array(ice%flx_fw         , unit, iostat, iomsg)
    call write_bin_array(ice%flx_h          , unit, iostat, iomsg)
    if (ice%whichEVP > 0) then
        call write_bin_array(ice%alpha_evp_array        , unit, iostat, iomsg)
        call write_bin_array(ice%beta_evp_array         , unit, iostat, iomsg)
    end if     
    
end subroutine WRITE_T_ICE

! Unformatted reading for T_ICE
subroutine READ_T_ICE(ice, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_ICE),           intent(inout)  :: ice
    integer,                intent(in)     :: unit
    integer,                intent(out)    :: iostat
    character(*),           intent(inout)  :: iomsg
    integer                                :: i
    
    !___________________________________________________________________________
    read(unit, iostat=iostat, iomsg=iomsg) ice%num_itracers
    if (.not. allocated(ice%data)) allocate(ice%data(ice%num_itracers))
    do i=1, ice%num_itracers
       read(unit, iostat=iostat, iomsg=iomsg) ice%data(i)
    end do
    !___________________________________________________________________________
    read(unit, iostat=iostat, iomsg=iomsg) ice%thermo
    read(unit, iostat=iostat, iomsg=iomsg) ice%work
#if defined (__oasis) || defined (__ifsinterface)
    read(unit, iostat=iostat, iomsg=iomsg) ice%atmcoupl
#endif /* (__oasis) */       

    !___________________________________________________________________________
    read(unit, iostat=iostat, iomsg=iomsg) ice%pstar
    read(unit, iostat=iostat, iomsg=iomsg) ice%ellipse
    read(unit, iostat=iostat, iomsg=iomsg) ice%c_pressure
    read(unit, iostat=iostat, iomsg=iomsg) ice%delta_min
    read(unit, iostat=iostat, iomsg=iomsg) ice%Clim_evp
    read(unit, iostat=iostat, iomsg=iomsg) ice%zeta_min
    read(unit, iostat=iostat, iomsg=iomsg) ice%evp_rheol_steps
    read(unit, iostat=iostat, iomsg=iomsg) ice%ice_gamma_fct
    read(unit, iostat=iostat, iomsg=iomsg) ice%ice_diff
    read(unit, iostat=iostat, iomsg=iomsg) ice%Tevp_inv
    read(unit, iostat=iostat, iomsg=iomsg) ice%theta_io
    read(unit, iostat=iostat, iomsg=iomsg) ice%alpha_evp
    read(unit, iostat=iostat, iomsg=iomsg) ice%beta_evp
    read(unit, iostat=iostat, iomsg=iomsg) ice%c_aevp
    read(unit, iostat=iostat, iomsg=iomsg) ice%ice_ave_steps
    read(unit, iostat=iostat, iomsg=iomsg) ice%cd_oce_ice
    read(unit, iostat=iostat, iomsg=iomsg) ice%ice_free_slip
    read(unit, iostat=iostat, iomsg=iomsg) ice%whichEVP
    read(unit, iostat=iostat, iomsg=iomsg) ice%ice_dt
    read(unit, iostat=iostat, iomsg=iomsg) ice%Tevp_inv
    read(unit, iostat=iostat, iomsg=iomsg) ice%ice_steps_since_upd
    read(unit, iostat=iostat, iomsg=iomsg) ice%ice_update
    
    !___________________________________________________________________________
    call read_bin_array(ice%uice            , unit, iostat, iomsg)
    call read_bin_array(ice%uice_rhs        , unit, iostat, iomsg)
    call read_bin_array(ice%uice_old        , unit, iostat, iomsg)
    if (ice%whichEVP /= 0) call read_bin_array(ice%uice_aux     , unit, iostat, iomsg)
    call read_bin_array(ice%vice            , unit, iostat, iomsg)
    call read_bin_array(ice%vice_rhs        , unit, iostat, iomsg)
    call read_bin_array(ice%vice_old        , unit, iostat, iomsg)
    if (ice%whichEVP /= 0) call read_bin_array(ice%vice_aux     , unit, iostat, iomsg)
    call read_bin_array(ice%stress_atmice_x , unit, iostat, iomsg)
    call read_bin_array(ice%stress_iceoce_x , unit, iostat, iomsg)
    call read_bin_array(ice%stress_atmice_y , unit, iostat, iomsg)
    call read_bin_array(ice%stress_iceoce_y , unit, iostat, iomsg)
    call read_bin_array(ice%srfoce_u        , unit, iostat, iomsg)
    call read_bin_array(ice%srfoce_v        , unit, iostat, iomsg)
    call read_bin_array(ice%srfoce_temp     , unit, iostat, iomsg)
    call read_bin_array(ice%srfoce_salt     , unit, iostat, iomsg)
    call read_bin_array(ice%srfoce_ssh      , unit, iostat, iomsg)
    call read_bin_array(ice%flx_fw          , unit, iostat, iomsg)
    call read_bin_array(ice%flx_h           , unit, iostat, iomsg)
    if (ice%whichEVP > 0) then
        call read_bin_array(ice%alpha_evp_array     , unit, iostat, iomsg)
        call read_bin_array(ice%beta_evp_array      , unit, iostat, iomsg)
    end if     
    
end subroutine READ_T_ICE
END MODULE MOD_ICE
!
!
!_______________________________________________________________________________
! interface to initialise derived type for sea ice 
module ice_init_interface
    interface
        subroutine ice_init(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARSUP
        USE MOD_PARTIT
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(inout), target :: mesh
        end subroutine
    end interface
end module
!
!
!_______________________________________________________________________________
! initialise derived type for sea ice 
subroutine ice_init(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    USE o_param, only: WP
    IMPLICIT NONE
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout), target :: mesh
    !___________________________________________________________________________
    integer        :: elem_size, node_size, n, ed(2)
    integer, save  :: nm_unit  = 105       ! unit to open namelist file, skip 100-102 for cray
    integer        :: iost
    !___________________________________________________________________________
    ! define ice namelist parameter
    integer        :: whichEVP, evp_rheol_steps, ice_ave_steps
    real(kind=WP)  :: Pstar, ellipse, c_pressure, delta_min, ice_gamma_fct, &
                      ice_diff, theta_io, alpha_evp, beta_evp, c_aevp, Cd_oce_ice
    namelist /ice_dyn/ whichEVP, Pstar, ellipse, c_pressure, delta_min, evp_rheol_steps, &
                       Cd_oce_ice, ice_gamma_fct, ice_diff, theta_io, ice_ave_steps, &
                       alpha_evp, beta_evp, c_aevp
                       
    real(kind=WP)  :: Sice, h0, emiss_ice, emiss_wat, albsn, albsnm, albi, &
                      albim, albw, con, consn                   
    namelist /ice_therm/ Sice, h0, emiss_ice, emiss_wat, albsn, albsnm, albi, &
                         albim, albw, con, consn
    !___________________________________________________________________________
    ! pointer on necessary derived types
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"   

    !___________________________________________________________________________
    ! open and read namelist.ice for I/O
        open(unit=nm_unit, file='namelist.ice', form='formatted', access='sequential', status='old', iostat=iost )
    if (iost == 0) then
        if (mype==0) write(*,*) '     file   : ', 'namelist.ice',' open ok'
    else
        if (mype==0) write(*,*) 'ERROR: --> bad opening file   : ', 'namelist.ice',' ; iostat=',iost
        call par_ex(partit%MPI_COMM_FESOM, partit%mype)
        stop
    end if
    read(nm_unit, nml=ice_dyn  , iostat=iost)
    read(nm_unit, nml=ice_therm, iostat=iost)
    close(nm_unit)
    
    !___________________________________________________________________________
    ! set parameters in ice derived type from namelist.ice --> namelist /ice_dyn/ 
    ice%whichEVP        = whichEVP 
    ice%pstar           = Pstar
    ice%ellipse         = ellipse
    ice%c_pressure      = c_pressure
    ice%delta_min       = delta_min
    ice%evp_rheol_steps = evp_rheol_steps
    ice%cd_oce_ice      = Cd_oce_ice
    ice%ice_gamma_fct   = ice_gamma_fct
    ice%ice_diff        = ice_diff
    ice%theta_io        = theta_io
    ice%ice_ave_steps   = ice_ave_steps
    ice%alpha_evp       = alpha_evp       
    ice%beta_evp        = beta_evp
    ice%c_aevp          = c_aevp
    
    ! set parameters in ice derived type from namelist.ice --> namelist /ice_therm/ 
    ice%thermo%con      = con
    ice%thermo%consn    = consn
    ice%thermo%Sice     = Sice
    ice%thermo%h0       = h0
    ice%thermo%emiss_ice= emiss_ice    
    ice%thermo%emiss_wat= emiss_wat
    ice%thermo%albsn    = albsn 
    ice%thermo%albsnm   = albsnm
    ice%thermo%albi     = albi    
    ice%thermo%albim    = albim
    ice%thermo%albw     = albw
    
    ice%thermo%cc=ice%thermo%rhowat*4190.0  ! Volumetr. heat cap. of water [J/m**3/K](cc = rhowat*cp_water)
    ice%thermo%cl=ice%thermo%rhoice*3.34e5  ! Volumetr. latent heat of ice fusion [J/m**3](cl=rhoice*Lf)
    
    !___________________________________________________________________________
    ! define local vertice & elem array size
    elem_size=myDim_elem2D+eDim_elem2D
    node_size=myDim_nod2D +eDim_nod2D
    
    !___________________________________________________________________________
    ! allocate/initialise arrays in ice derived type
    ! initialise velocity and stress related arrays in ice derived type 
    allocate(ice%uice(                 node_size))
    allocate(ice%uice_rhs(             node_size))
    allocate(ice%uice_old(             node_size))
    allocate(ice%vice(                 node_size))
    allocate(ice%vice_rhs(             node_size))
    allocate(ice%vice_old(             node_size))
    allocate(ice%stress_atmice_x(      node_size))
    allocate(ice%stress_iceoce_x(      node_size))
    allocate(ice%stress_atmice_y(      node_size))
    allocate(ice%stress_iceoce_y(      node_size))
    ice%uice            = 0.0_WP
    ice%uice_rhs        = 0.0_WP
    ice%uice_old        = 0.0_WP
    ice%stress_atmice_x = 0.0_WP
    ice%stress_iceoce_x = 0.0_WP
    ice%vice            = 0.0_WP
    ice%vice_rhs        = 0.0_WP
    ice%vice_old        = 0.0_WP
    ice%stress_atmice_y = 0.0_WP
    ice%stress_iceoce_y = 0.0_WP
    if (ice%whichEVP /= 0) then
        allocate(ice%uice_aux(         node_size))
        allocate(ice%vice_aux(         node_size))
        ice%uice_aux    = 0.0_WP
        ice%vice_aux    = 0.0_WP
    end if
    if (ice%whichEVP == 2) then
        allocate(ice%alpha_evp_array(  node_size))
        allocate(ice%beta_evp_array(   node_size))
        ice%alpha_evp_array = ice%alpha_evp
        ice%beta_evp_array  = ice%alpha_evp
    end if
    
    !___________________________________________________________________________
    ! initialise surface ocean arrays in ice derived type 
    allocate(ice%srfoce_u(             node_size))
    allocate(ice%srfoce_v(             node_size))
    allocate(ice%srfoce_temp(          node_size))
    allocate(ice%srfoce_salt(          node_size))
    allocate(ice%srfoce_ssh(           node_size))
    ice%srfoce_u         = 0.0_WP
    ice%srfoce_v         = 0.0_WP
    ice%srfoce_temp      = 0.0_WP
    ice%srfoce_salt      = 0.0_WP
    ice%srfoce_ssh       = 0.0_WP
    
    allocate(ice%flx_fw(node_size))
    allocate(ice%flx_h( node_size))
    ice%flx_fw           = 0.0_WP
    ice%flx_h            = 0.0_WP
    
    !___________________________________________________________________________
    ! initialse data array of ice derived type containing "ice tracer" that have
    ! to be advected: a_ice (index=1), m_ice (index=2), m_snow (index=3), 
    ! ice_temp (index=4, only when coupled)
    allocate(ice%data(ice%num_itracers))
    do n = 1, ice%num_itracers
        allocate(ice%data(n)%values(    node_size))
        allocate(ice%data(n)%values_old(node_size))
        allocate(ice%data(n)%values_rhs(node_size))
        allocate(ice%data(n)%values_div_rhs(node_size))
        allocate(ice%data(n)%dvalues(   node_size))
        allocate(ice%data(n)%valuesl(   node_size))
        ice%data(n)%ID             = n
        ice%data(n)%values         = 0.0_WP
        ice%data(n)%values_old     = 0.0_WP
        ice%data(n)%values_rhs     = 0.0_WP
        ice%data(n)%values_div_rhs = 0.0_WP
        ice%data(n)%dvalues        = 0.0_WP
        ice%data(n)%valuesl        = 0.0_WP
        if (n==4) ice%data(n)%values = 265.15_WP
    end do
    
    !___________________________________________________________________________
    ! initialse work array of ice derived type 
    allocate(ice%work%fct_tmax(        node_size))
    allocate(ice%work%fct_tmin(        node_size))
    allocate(ice%work%fct_plus(        node_size))
    allocate(ice%work%fct_minus(       node_size))
    allocate(ice%work%fct_fluxes(      elem_size, 3))
    ice%work%fct_tmax    = 0.0_WP
    ice%work%fct_tmin    = 0.0_WP
    ice%work%fct_plus    = 0.0_WP
    ice%work%fct_minus   = 0.0_WP
    ice%work%fct_fluxes  = 0.0_WP
    
    allocate(ice%work%fct_massmatrix(sum(nn_num(1:myDim_nod2D))))
    ice%work%fct_massmatrix = 0.0_WP
    
    allocate(ice%work%sigma11(         elem_size))
    allocate(ice%work%sigma12(         elem_size))
    allocate(ice%work%sigma22(         elem_size))
    allocate(ice%work%eps11(           elem_size))
    allocate(ice%work%eps12(           elem_size))
    allocate(ice%work%eps22(           elem_size))
    ice%work%sigma11     = 0.0_WP
    ice%work%sigma12     = 0.0_WP
    ice%work%sigma22     = 0.0_WP
    ice%work%eps11       = 0.0_WP
    ice%work%eps12       = 0.0_WP
    ice%work%eps22       = 0.0_WP
    
    allocate(ice%work%ice_strength(    elem_size))
    allocate(ice%work%inv_areamass(    node_size))
    allocate(ice%work%inv_mass(        node_size))
    ice%work%ice_strength= 0.0_WP
    ice%work%inv_areamass= 0.0_WP
    ice%work%inv_mass    = 0.0_WP
    
    !___________________________________________________________________________
    ! initialse thermo array of ice derived type 
    allocate(ice%thermo%ustar(         node_size))
    allocate(ice%thermo%t_skin(        node_size))
    allocate(ice%thermo%thdgr(         node_size))
    allocate(ice%thermo%thdgrsn(       node_size))
    allocate(ice%thermo%thdgr_old(     node_size))
    ice%thermo%ustar     = 0.0_WP
    ice%thermo%t_skin    = 0.0_WP
    ice%thermo%thdgr     = 0.0_WP
    ice%thermo%thdgrsn   = 0.0_WP
    ice%thermo%thdgr_old = 0.0_WP
    
    !___________________________________________________________________________
    ! initialse coupling array of ice derived type 
#if defined (__oasis) || defined (__ifsinterface)    
    allocate(ice%atmcoupl%oce_flx_h(     node_size))
    allocate(ice%atmcoupl%ice_flx_h(     node_size))
    allocate(ice%atmcoupl%tmpoce_flx_h(  node_size))
    allocate(ice%atmcoupl%tmpice_flx_h(  node_size))
    ice%atmcoupl%oce_flx_h     = 0.0_WP
    ice%atmcoupl%ice_flx_h     = 0.0_WP
    ice%atmcoupl%tmpoce_flx_h  = 0.0_WP
    ice%atmcoupl%tmpice_flx_h  = 0.0_WP
#if defined (__oifs) || defined (__ifsinterface)  
    allocate(ice%atmcoupl%ice_alb(       node_size))
    allocate(ice%atmcoupl%enthalpyoffuse(node_size))
    ice%atmcoupl%ice_alb       = 0.6_WP
    ice%atmcoupl%enthalpyoffuse= 0.0_WP
#endif /* (__oifs) */
#endif /* (__oasis) */       

    !___________________________________________________________________________
    ! --> took from oce_mesh.F90 --> subroutine mesh_auxiliary_arrays(partit, mesh)
    ! to here since namelist.ice is now read in ice_init where whichEVP is not available
    ! when  mesh_auxiliary_arrays is called
    !array of 2D boundary conditions is used in ice_maEVP
    if (ice%whichEVP > 0) then
        allocate(mesh%bc_index_nod2D(myDim_nod2D+eDim_nod2D))
        mesh%bc_index_nod2D=1._WP
        do n=1, myDim_edge2D
            ed=mesh%edges(:, n)
            if (myList_edge2D(n) <= mesh%edge2D_in) cycle
            mesh%bc_index_nod2D(ed)=0._WP
        end do
    end if
    
end subroutine ice_init  
!
!
!
!
!
!_______________________________________________________________________________
! initialise derived type for sea ice 
subroutine ice_init_toyocean_dummy(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    USE o_param, only: WP
    IMPLICIT NONE
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer        :: node_size, n
    !___________________________________________________________________________
    ! pointer on necessary derived types
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"   

    !___________________________________________________________________________
    ! define local vertice & elem array size
    node_size=myDim_nod2D+eDim_nod2D
    
    !___________________________________________________________________________
    ! allocate/initialise arrays in ice derived type
    ! initialise velocity and stress related arrays in ice derived type 
    allocate(ice%uice(                  node_size))
    allocate(ice%vice(                  node_size))
    ice%uice               = 0.0_WP
    ice%vice               = 0.0_WP
    allocate(ice%data(ice%num_itracers))
    do n = 1, ice%num_itracers
        allocate(ice%data(n)%values(    node_size))
        ice%data(n)%ID      = n
        ice%data(n)%values  = 0.0_WP
    end do
end subroutine ice_init_toyocean_dummy    
