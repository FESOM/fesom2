MODULE MOD_ICE
USE I_PARAM
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
                                                   values_div_rhs, dvalues, values1
    integer                                     :: ID
    !___________________________________________________________________________
    contains
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
    real(kind=WP), allocatable, dimension(:)    :: tmax, tmin
    !___________________________________________________________________________
    contains
        procedure WRITE_T_ICE_WORK
        procedure READ_T_ICE_WORK
        generic :: write(unformatted) => WRITE_T_ICE_WORK
        generic :: read(unformatted)  => READ_T_ICE_WORK
END TYPE T_ICE_WORK
!
! 
!_______________________________________________________________________________
! set work array derived type for ice
TYPE T_ICE_ATMCOUPL
#if defined (__oasis)
    !___________________________________________________________________________
    real(kind=WP), allocatable, dimension(:)    :: oce_flx_h, ice_flx_h, tmpoce_flx_h, tmpice_flx_h
#if defined (__oifs)
    !___________________________________________________________________________
    real(kind=WP), allocatable, dimension(:)    :: ice_alb, enthalpyoffuse
    ! !!! DONT FORGET ice_temp rhs_tempdiv rhs_temp is advected for oifs !!! --> becomes additional ice 
    ! tracer in ice%data(4)%values
#endif /* (__oifs) */
#endif /* (__oasis) */
    
    !___________________________________________________________________________
    contains
        procedure WRITE_T_ICE_ATMCOUPL
        procedure READ_T_ICE_ATMCOUPL
        generic :: write(unformatted) => WRITE_T_ICE_ATMCOUPL
        generic :: read(unformatted)  => READ_T_ICE_ATMCOUPL
END TYPE T_ICE_ATMCOUPL

!
! 
!_______________________________________________________________________________
! set main ice derived type contains parameters, data array, work array, u_ice, vice
TYPE T_ICE

    !___________________________________________________________________________
    ! zonal & merdional ice velocity
    real(kind=WP), allocatable, dimension(:,:)  :: uvice, uvice_rhs, uvice_old
    
    ! surface stess atm<-->ice, oce<-->ice
    real(kind=WP), allocatable, dimension(:,:)  :: stress_atmice_xy, stress_iceoce_xy
    
    ! oce temp, salt, ssh, and uv at surface
    real(kind=WP), allocatable, dimension(:)    :: srfoce_temp, srfoce_salt, srfoce_ssh
    real(kind=WP), allocatable, dimension(:,:)  :: srfoce_uv
    
    ! freshwater & heatflux
    real(kind=WP), allocatable, dimension(:)    :: flx_fw, flx_h
    
    !___________________________________________________________________________
    ! total number of ice tracers (default=3, 1=area, 2=mice, 3=msnow, (4=ice_temp)
#if defined (__oifs)
    integer                                     :: num_itracers=4
#else
    integer                                     :: num_itracers=3
#endif 

    ! put ice tracers data arrays
    type(t_ice_data), allocatable, dimension(:) :: data
    
    !___________________________________________________________________________
    ! put ice working arrays
    type(t_ice_work)                            :: work
    
#if defined (__oasis)
    !___________________________________________________________________________
    ! put ice arrays for coupled model
    type(t_ice_atmcoupl)                        :: atmcoupl
#endif /* (__oasis) */    
    !___________________________________________________________________________
    ! set ice model parameters:
    ! --- RHEOLOGY ---
    real(kind=WP)             :: Pstar      = 30000.0_WP   ![N/m^2]
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
    
    !___________________________________________________________________________
    contains
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
    call write_bin_array(tdata%values1,        unit, iostat, iomsg)
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
    call read_bin_array(tdata%values1,        unit, iostat, iomsg)
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
    call write_bin_array(twork%tmin, unit, iostat, iomsg)
    call write_bin_array(twork%tmax, unit, iostat, iomsg)
end subroutine WRITE_T_ICE_WORK    

! Unformatted reading for T_ICE_WORK
subroutine READ_T_ICE_WORK(twork, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_ICE_WORK),      intent(inout)  :: twork
    integer,                intent(in)     :: unit
    integer,                intent(out)    :: iostat
    character(*),           intent(inout)  :: iomsg
    call read_bin_array(twork%tmin, unit, iostat, iomsg)
    call read_bin_array(twork%tmax, unit, iostat, iomsg)
end subroutine READ_T_ICE_WORK
!
!
!_______________________________________________________________________________
! Unformatted writing for T_ICE_ATMCOUPL
subroutine WRITE_T_ICE_ATMCOUPL(tcoupl, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_ICE_ATMCOUPL),  intent(in)     :: tcoupl
    integer,                intent(in)     :: unit
    integer,                intent(out)    :: iostat
    character(*),           intent(inout)  :: iomsg
#if defined (__oasis)    
    call write_bin_array(tcoupl%oce_flx_h, unit, iostat, iomsg)
    call write_bin_array(tcoupl%ice_flx_h, unit, iostat, iomsg)
    call write_bin_array(tcoupl%tmpoce_flx_h, unit, iostat, iomsg)
    call write_bin_array(tcoupl%tmpice_flx_h, unit, iostat, iomsg)
#if defined (__oifs)    
    call write_bin_array(tcoupl%ice_alb, unit, iostat, iomsg)
    call write_bin_array(tcoupl%enthalpyoffuse, unit, iostat, iomsg)
#endif /* (__oifs) */
#endif /* (__oasis) */    
end subroutine WRITE_T_ICE_ATMCOUPL  

! Unformatted reading for T_ICE_ATMCOUPL
subroutine READ_T_ICE_ATMCOUPL(tcoupl, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_ICE_ATMCOUPL),  intent(inout)  :: tcoupl
    integer,                intent(in)     :: unit
    integer,                intent(out)    :: iostat
    character(*),           intent(inout)  :: iomsg
#if defined (__oasis)    
    call read_bin_array(tcoupl%oce_flx_h, unit, iostat, iomsg)
    call read_bin_array(tcoupl%ice_flx_h, unit, iostat, iomsg)
    call read_bin_array(tcoupl%tmpoce_flx_h, unit, iostat, iomsg)
    call read_bin_array(tcoupl%tmpice_flx_h, unit, iostat, iomsg)
#if defined (__oifs)    
    call read_bin_array(tcoupl%ice_alb, unit, iostat, iomsg)
    call read_bin_array(tcoupl%enthalpyoffuse, unit, iostat, iomsg)
#endif /* (__oifs) */
#endif /* (__oasis) */   
end subroutine READ_T_ICE_ATMCOUPL
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
    call write_bin_array(ice%uvice, unit, iostat, iomsg)
    call write_bin_array(ice%uvice_rhs, unit, iostat, iomsg)
    call write_bin_array(ice%uvice_old, unit, iostat, iomsg)
    call write_bin_array(ice%stress_atmice_xy, unit, iostat, iomsg)
    call write_bin_array(ice%stress_iceoce_xy, unit, iostat, iomsg)
    call write_bin_array(ice%srfoce_temp, unit, iostat, iomsg)
    call write_bin_array(ice%srfoce_salt, unit, iostat, iomsg)
    call write_bin_array(ice%srfoce_ssh, unit, iostat, iomsg)
    call write_bin_array(ice%flx_fw, unit, iostat, iomsg)
    call write_bin_array(ice%flx_h, unit, iostat, iomsg)
    
    !___________________________________________________________________________
    write(unit, iostat=iostat, iomsg=iomsg) ice%num_itracers
    do i=1, ice%num_itracers
       write(unit, iostat=iostat, iomsg=iomsg) ice%data(i)
    end do
    !___________________________________________________________________________
    write(unit, iostat=iostat, iomsg=iomsg) ice%work
#if defined (__oasis)
    write(unit, iostat=iostat, iomsg=iomsg) ice%atmcoupl
#endif /* (__oasis) */       

    !___________________________________________________________________________
    write(unit, iostat=iostat, iomsg=iomsg) ice%Pstar
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
    call read_bin_array(ice%uvice, unit, iostat, iomsg)
    call read_bin_array(ice%uvice_rhs, unit, iostat, iomsg)
    call read_bin_array(ice%uvice_old, unit, iostat, iomsg)
    call read_bin_array(ice%stress_atmice_xy, unit, iostat, iomsg)
    call read_bin_array(ice%stress_iceoce_xy, unit, iostat, iomsg)
    call read_bin_array(ice%srfoce_temp, unit, iostat, iomsg)
    call read_bin_array(ice%srfoce_salt, unit, iostat, iomsg)
    call read_bin_array(ice%srfoce_ssh, unit, iostat, iomsg)
    call read_bin_array(ice%flx_fw, unit, iostat, iomsg)
    call read_bin_array(ice%flx_h, unit, iostat, iomsg)
    
    !___________________________________________________________________________
    read(unit, iostat=iostat, iomsg=iomsg) ice%num_itracers
    do i=1, ice%num_itracers
       read(unit, iostat=iostat, iomsg=iomsg) ice%data(i)
    end do
    !___________________________________________________________________________
    read(unit, iostat=iostat, iomsg=iomsg) ice%work
#if defined (__oasis)
    read(unit, iostat=iostat, iomsg=iomsg) ice%atmcoupl
#endif /* (__oasis) */       

    !___________________________________________________________________________
    read(unit, iostat=iostat, iomsg=iomsg) ice%Pstar
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
    
end subroutine READ_T_ICE






END MODULE MOD_ICE