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
    real(kind=WP), allocatable, dimension(:)    :: t_skin, thdgr, thdgrsn, thdgr_old
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
TYPE T_ICE_THERMO
    !___________________________________________________________________________
    real(kind=WP), allocatable, dimension(:)    :: t_skin, thdgr, thdgrsn, thdgr_old, ustar
    !___________________________________________________________________________
    contains
        procedure WRITE_T_ICE_THERMO
        procedure READ_T_ICE_THERMO
        generic :: write(unformatted) => WRITE_T_ICE_THERMO
        generic :: read(unformatted)  => READ_T_ICE_THERMO
END TYPE T_ICE_THERMO
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
    real(kind=WP), allocatable, dimension(:,:)  :: uvice, uvice_rhs, uvice_old, uvice_aux
    
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
    
    ! put thermodynamics arrays
    type(t_ice_thermo)                          :: thermo
    
#if defined (__oasis)
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
subroutine WRITE_T_ICE_ATMCOUPL(tcoupl, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_ICE_ATMCOUPL),  intent(in)     :: tcoupl
    integer,                intent(in)     :: unit
    integer,                intent(out)    :: iostat
    character(*),           intent(inout)  :: iomsg
#if defined (__oasis)    
    call write_bin_array(tcoupl%oce_flx_h,      unit, iostat, iomsg)
    call write_bin_array(tcoupl%ice_flx_h,      unit, iostat, iomsg)
    call write_bin_array(tcoupl%tmpoce_flx_h,   unit, iostat, iomsg)
    call write_bin_array(tcoupl%tmpice_flx_h,   unit, iostat, iomsg)
#if defined (__oifs)    
    call write_bin_array(tcoupl%ice_alb,        unit, iostat, iomsg)
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
    if (ice%whichEVP /= 0) call write_bin_array(ice%uvice_aux, unit, iostat, iomsg)
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
    write(unit, iostat=iostat, iomsg=iomsg) ice%thermo
    write(unit, iostat=iostat, iomsg=iomsg) ice%work
#if defined (__oasis)
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
    if (ice%whichEVP /= 0) call read_bin_array(ice%uvice_aux, unit, iostat, iomsg)
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
    read(unit, iostat=iostat, iomsg=iomsg) ice%thermo
    read(unit, iostat=iostat, iomsg=iomsg) ice%work
#if defined (__oasis)
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
        type(t_mesh)  , intent(in)   , target :: mesh
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
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer        :: elem_size, node_size, n
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
    read(nm_unit, nml=ice_dyn,    iostat=iost)
    close(nm_unit)
    
    !___________________________________________________________________________
    ! set parameters in ice derived type from namelist.ice
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
    
    !!PS no namelist paramter  in moment 
    !!PS ice%zeta_min        = zeta_min
    !!PS ice%Tevp_inv        = Tevp_inv
    !!PS ice%ice_free_slip   = ice_free_slip
    !!PS ice%ice_dt          = ice_dt
    !!PS ice%Tevp_inv        = Tevp_inv
    
    !___________________________________________________________________________
    ! define local vertice & elem array size
    elem_size=myDim_elem2D+eDim_elem2D
    node_size=myDim_nod2D+eDim_nod2D
    
    !___________________________________________________________________________
    ! allocate/initialise arrays in ice derived type
    ! initialise velocity and stress related arrays in ice derived type 
    allocate(ice%uvice(             2, node_size))
    allocate(ice%uvice_rhs(         2, node_size))
    allocate(ice%uvice_old(         2, node_size))
    allocate(ice%stress_atmice_xy(  2, node_size))
    allocate(ice%stress_iceoce_xy(  2, node_size))
    ice%uvice            = 0.0_WP
    ice%uvice_rhs        = 0.0_WP
    ice%uvice_old        = 0.0_WP
    ice%stress_atmice_xy = 0.0_WP
    ice%stress_iceoce_xy = 0.0_WP
    if (ice%whichEVP /= 0) then
        allocate(ice%uvice_aux(     2, node_size))
        ice%uvice_aux    = 0.0_WP
    end if 
    
    !___________________________________________________________________________
    ! initialise surface ocean arrays in ice derived type 
    allocate(ice%srfoce_uv(         2, node_size))
    allocate(ice%srfoce_temp(          node_size))
    allocate(ice%srfoce_salt(          node_size))
    allocate(ice%srfoce_ssh(           node_size))
    ice%srfoce_uv        = 0.0_WP
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
    
    allocate(ice%work%t_skin(          node_size))
    allocate(ice%work%thdgr(           node_size))
    allocate(ice%work%thdgrsn(         node_size))
    allocate(ice%work%thdgr_old(       node_size))
    ice%work%t_skin      = 0.0_WP
    ice%work%thdgr       = 0.0_WP
    ice%work%thdgrsn     = 0.0_WP
    ice%work%thdgr_old   = 0.0_WP
    
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
#if defined (__oasis)    
    allocate(ice%tcoupl%oce_flx_h(     node_size))
    allocate(ice%tcoupl%ice_flx_h(     node_size))
    allocate(ice%tcoupl%tmpoce_flx_h(  node_size))
    allocate(ice%tcoupl%tmpice_flx_h(  node_size))
    ice%tcoupl%oce_flx_h     = 0.0_WP
    ice%tcoupl%ice_flx_h     = 0.0_WP
    ice%tcoupl%tmpoce_flx_h  = 0.0_WP
    ice%tcoupl%tmpice_flx_h  = 0.0_WP
#if defined (__oifs)  
    allocate(ice%tcoupl%ice_alb(       node_size))
    allocate(ice%tcoupl%enthalpyoffuse(node_size))
    ice%tcoupl%ice_alb       = 0.0_WP
    ice%tcoupl%enthalpyoffuse= 0.0_WP
#endif /* (__oifs) */
#endif /* (__oasis) */        
end subroutine ice_init    