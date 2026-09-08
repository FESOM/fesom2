MODULE MOD_ICE
USE o_PARAM, only: WP
USE, intrinsic :: ISO_FORTRAN_ENV, only : int32
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
    real(kind=WP), allocatable, dimension(:)    :: inv_areamass, inv_mass
    ! ice_strength: canonical Hibler (1979) ice strength P = P*·h·exp(-C(1-A))
    ! evaluated at nodes from m_ice and a_ice. Exposed to def_stream as the
    ! 'strength_ice' output. EVP/mEVP/aEVP all evaluate the per-element strength
    ! inline from this nodal field; no per-element cache is kept.
    real(kind=WP), allocatable, dimension(:)    :: ice_strength
    !___________________________________________________________________________
    contains
        procedure WRITE_T_ICE_WORK
        procedure READ_T_ICE_WORK
END TYPE T_ICE_WORK
!
!
!_______________________________________________________________________________
! set work array derived type for ice
TYPE T_ICE_THERMO
    !___________________________________________________________________________
    real(kind=WP), allocatable, dimension(:)    :: t_skin, thdgr, thdgrsn, thdgra, thdgr_old, ustar
    real(kind=WP), allocatable, dimension(:)    :: dyngr, dyngrsn, dyngra ! dynamic growth of: ice, snow, area (cmip6 variable for letti!)
    ! melt pond variables
    real(kind=WP), allocatable, dimension(:)    :: apnd, hpnd, ipnd  ! pond area fraction, depth, ice thickness
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
    real(kind=WP) :: cc=1025.*4190.0   ! Volumetr. heat cap. of water [J/m**3/K](cc = rhowat*cp_water)
    real(kind=WP) :: cl=910.*3.34e5    ! Volumetr. latent heat of ice fusion [J/m**3](cl=rhoice*Lf)
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
    real(kind=WP) :: h0=0.5	           ! Lead closing parameter [m] for Nothern Hemisphere! 0.5
    real(kind=WP) :: h0_s=0.5	       ! Lead closing parameter [m] for Southern Hemisphere! 0.5
    real(kind=WP) :: emiss_ice=0.97    ! Emissivity of Snow/Ice,
    real(kind=WP) :: emiss_wat=0.97    ! Emissivity of open water
    real(kind=WP) :: albsn = 0.81      ! Albedo: frozen snow
    real(kind=WP) :: albsnm= 0.77      !         melting snow
    real(kind=WP) :: albi  = 0.70      !         frozen ice
    real(kind=WP) :: albim = 0.68      !         melting ice
    real(kind=WP) :: albw  = 0.066     !         open water, LY2004
    ! Smooth snow→bare-ice albedo blend; alpha_snow = tanh(hsn/h_snowscale).
    ! 0 (default) keeps the legacy step function at hsn > 1 mm.
    ! Typical CICE-style values are 0.02-0.05 m.
    real(kind=WP) :: h_snowscale = 0.0_WP
    ! Coupled snow melt on sea ice requires a surface temperature above 273 K
    ! (.true., legacy). .false. melts snow on the energy budget alone, as the
    ! uncoupled branch does.
    logical       :: snowmelt_tgate = .true.
    ! Width [K] of a linear transition from frozen to melting snow and ice
    ! albedo below 273.15 K. 0 (default) keeps the step at 273.15 K.
    real(kind=WP) :: alb_tramp = 0.0_WP
    real(kind=WP) :: h_ml  = 2.5_WP    ! thickness of uppermost layer deacides how much heat is available

    ! --- additional namelist parameters (Frank.Kauker(at)awi.de 2023/04/04)
    logical       :: snowdist=.true.   ! distribution of snow depth according to ice distribution
    logical       :: new_iclasses=.false. ! ice thickness distribution based on EM observations (Castro-Morales et al., JGR, 2013)
    integer       :: open_water_albedo=0  ! 0=standard; 1=taylor; 2=briegleb
    REAL(kind=WP) :: c_melt=0.5        ! constant in concentration equation for melting conditions
    ! --- melt pond parameters
    logical       :: use_meltponds=.false. ! enable melt pond parameterization
    REAL(kind=WP) :: h_cutoff=3.0      ! cutoff thickness of thickness pdf
    REAL(kind=WP), DIMENSION(15) :: hpdf = (/ 0.066745491, 0.1462317, 0.17769822, 0.13131106, &
         0.11518432, 0.08514193, 0.06871303, 0.05592151, 0.04428673, 0.03584652, 0.02970195, 0.02469673, &
         0.02001543, 0.01653681, 0.0141026 /)  ! pdf of ice thickness based on EM observations
    contains
        procedure WRITE_T_ICE_THERMO
        procedure READ_T_ICE_THERMO
END TYPE T_ICE_THERMO
!
!
!_______________________________________________________________________________
! set work array derived type for ice
#if defined (__cpl_enabled)
TYPE T_ICE_ATMCOUPL

    !___________________________________________________________________________
    real(kind=WP), allocatable, dimension(:)    :: oce_flx_h, ice_flx_h, tmpoce_flx_h, tmpice_flx_h
    !___________________________________________________________________________
    ! Needed by the IFS-family partners. Declared and allocated for every
    ! coupled build because the partner is chosen at run time.
    real(kind=WP), allocatable, dimension(:)    :: ice_alb, enthalpyoffuse, runoff_liquid, runoff_solid, flx_qres, flx_qcon
    ! ist anchor: ice surface temperature as ACTUALLY TRANSMITTED at the last
    ! OASIS send -- the temperature OIFS evaluates its ice-tile fluxes at for
    ! the coming coupling interval. Consumed by the implicit (dQ/dT-linearized)
    ! surface-temperature solve in ice_thermo_cpl.F90 (ice_surftemp).
    real(kind=WP), allocatable, dimension(:)    :: ist_ref
    ! ice_temp / rhs_temp / rhs_tempdiv are advected as the additional ice
    ! tracer in ice%data(ist_itracer_idx)%values.
    !___________________________________________________________________________
    contains
        procedure WRITE_T_ICE_ATMCOUPL
        procedure READ_T_ICE_ATMCOUPL
END TYPE T_ICE_ATMCOUPL
#endif /* (__cpl_enabled) */

!
!
!_______________________________________________________________________________
! set main ice derived type contains parameters, data array, work array, u_ice, vice
TYPE T_ICE

    !___________________________________________________________________________
    ! zonal & merdional ice velocity
    real(kind=WP), allocatable, dimension(:)    :: uice, uice_rhs, uice_old, uice_aux, uice_ib
    real(kind=WP), allocatable, dimension(:)    :: vice, vice_rhs, vice_old, vice_aux, vice_ib
    
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
    ! ice/snow thicknesses in the ice-covered area
    real(kind=WP), allocatable, dimension(:)    :: h_ice, h_snow    

    !___________________________________________________________________________
    ! Ice tracers: 1=area, 2=mice, 3=msnow, then the optional prognostic ice
    ! surface temperature, then the two asynchronous-iceberg tracers (which
    ! icb_dyn.F90 addresses as the last two slots, so they follow the count).
    !
    ! Both are set in ice_init, from the partner atmosphere.
    ! ist_itracer_idx is the slot holding the ice surface temperature, or 0
    ! when the partner does not expect FESOM to carry one. The values here
    ! only apply until then.
    integer                                     :: num_itracers=3
    integer                                     :: ist_itracer_idx=0

    ! put ice tracers data arrays
    type(t_ice_data), allocatable, dimension(:) :: data

    !___________________________________________________________________________
    ! put ice working arrays
    type(t_ice_work)                            :: work

    ! put thermodynamics arrays
    type(t_ice_thermo)                          :: thermo
    
#if defined (__cpl_enabled)
    !___________________________________________________________________________
    ! put ice arrays for coupled model
    type(t_ice_atmcoupl)                        :: atmcoupl
#endif /* (__cpl_enabled) */

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
        procedure, private WRITE_T_ICE
        procedure, private READ_T_ICE
#else
        procedure WRITE_T_ICE
        procedure READ_T_ICE
#endif
        generic :: write(unformatted) => WRITE_T_ICE
        generic :: read(unformatted)  => READ_T_ICE
END TYPE T_ICE

contains

!
!_______________________________________________________________________________
!> Ice tracer layout: how many tracers, and which slot (if any) carries the
!> prognostic ice surface temperature.
!>
!> The IFS-family atmospheres expect FESOM to carry that temperature and
!> advect it; the others do not. The three counts reproduce what the
!> preprocessor used to pick: 6 with the temperature tracer, otherwise 5 with
!> the asynchronous iceberg tracers and 3 without. icb_dyn.F90 addresses the
!> iceberg tracers as the last two slots, so they follow the count.
subroutine set_ice_tracer_layout(ice)
    use cpl_config, only: is_coupled_to_oifs, is_coupled_to_ifs
    implicit none
    type(t_ice), intent(inout) :: ice

    if (is_coupled_to_oifs .or. is_coupled_to_ifs) then
        ice%ist_itracer_idx = 4
        ice%num_itracers    = 6
    else
        ice%ist_itracer_idx = 0
#if defined(__async_icebergs)
        ice%num_itracers    = 5
#else
        ice%num_itracers    = 3
#endif
    end if
end subroutine set_ice_tracer_layout
!
!
!_______________________________________________________________________________
! Unformatted writing for T_ICE_DATA
subroutine WRITE_T_ICE_DATA(tdata, unit)
    IMPLICIT NONE
    class(T_ICE_DATA),      intent(in)     :: tdata
    integer,                intent(in)     :: unit
    integer                                :: iostat
    character(len=1024)                    :: iomsg
    call write_bin_array(tdata%values,         unit, iostat, iomsg)
    call write_bin_array(tdata%values_old,     unit, iostat, iomsg)
    call write_bin_array(tdata%values_rhs,     unit, iostat, iomsg)
    call write_bin_array(tdata%values_div_rhs, unit, iostat, iomsg)
    call write_bin_array(tdata%dvalues,        unit, iostat, iomsg)
    call write_bin_array(tdata%valuesl,        unit, iostat, iomsg)
    write(unit, iostat=iostat, iomsg=iomsg) tdata%ID
end subroutine WRITE_T_ICE_DATA

! Unformatted reading for T_ICE_DATA
subroutine READ_T_ICE_DATA(tdata, unit)
    IMPLICIT NONE
    class(T_ICE_DATA),      intent(inout)  :: tdata
    integer,                intent(in)     :: unit
    integer                                :: iostat
    character(len=1024)                    :: iomsg
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
subroutine WRITE_T_ICE_WORK(twork, unit)
    IMPLICIT NONE
    class(T_ICE_WORK),      intent(in)     :: twork
    integer,                intent(in)     :: unit
    integer                                :: iostat
    character(len=1024)                    :: iomsg
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
    call write_bin_array(twork%inv_areamass, unit, iostat, iomsg)
    call write_bin_array(twork%inv_mass,     unit, iostat, iomsg)
end subroutine WRITE_T_ICE_WORK

! Unformatted reading for T_ICE_WORK
subroutine READ_T_ICE_WORK(twork, unit)
    IMPLICIT NONE
    class(T_ICE_WORK),      intent(inout)  :: twork
    integer,                intent(in)     :: unit
    integer                                :: iostat
    character(len=1024)                    :: iomsg
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
    call read_bin_array(twork%inv_areamass, unit, iostat, iomsg)
    call read_bin_array(twork%inv_mass,     unit, iostat, iomsg)
end subroutine READ_T_ICE_WORK
!
!
!_______________________________________________________________________________
! Unformatted writing for T_ICE_WORK
subroutine WRITE_T_ICE_THERMO(ttherm, unit)
    IMPLICIT NONE
    class(T_ICE_THERMO),    intent(in)     :: ttherm
    integer,                intent(in)     :: unit
    integer                                :: iostat
    character(len=1024)                    :: iomsg
    call write_bin_array(ttherm%t_skin,       unit, iostat, iomsg)
    call write_bin_array(ttherm%thdgr,        unit, iostat, iomsg)
    call write_bin_array(ttherm%thdgrsn,      unit, iostat, iomsg)
    call write_bin_array(ttherm%thdgr_old,    unit, iostat, iomsg)
    call write_bin_array(ttherm%ustar,        unit, iostat, iomsg)
    
    ! melt pond variables
    call write_bin_array(ttherm%apnd,         unit, iostat, iomsg)
    call write_bin_array(ttherm%hpnd,         unit, iostat, iomsg)
    call write_bin_array(ttherm%ipnd,         unit, iostat, iomsg)
    
    ! dynamic growth of: ice, snow, area (cmip6 variable for letti!)
    call write_bin_array(ttherm%dyngr,        unit, iostat, iomsg)
    call write_bin_array(ttherm%dyngrsn,      unit, iostat, iomsg)
    call write_bin_array(ttherm%dyngra,       unit, iostat, iomsg)
    
end subroutine WRITE_T_ICE_THERMO

! Unformatted reading for T_ICE_WORK
subroutine READ_T_ICE_THERMO(ttherm, unit)
    IMPLICIT NONE
    class(T_ICE_THERMO),    intent(inout)  :: ttherm
    integer,                intent(in)     :: unit
    integer                                :: iostat
    character(len=1024)                    :: iomsg
    call read_bin_array(ttherm%t_skin,       unit, iostat, iomsg)
    call read_bin_array(ttherm%thdgr,        unit, iostat, iomsg)
    call read_bin_array(ttherm%thdgrsn,      unit, iostat, iomsg)
    call read_bin_array(ttherm%thdgr_old,    unit, iostat, iomsg)
    call read_bin_array(ttherm%ustar,        unit, iostat, iomsg)
    ! melt pond variables
    call read_bin_array(ttherm%apnd,         unit, iostat, iomsg)
    call read_bin_array(ttherm%hpnd,         unit, iostat, iomsg)
    call read_bin_array(ttherm%ipnd,         unit, iostat, iomsg)
    
    ! dynamic growth of: ice, snow, area (cmip6 variable for letti!)
    call read_bin_array(ttherm%dyngr,        unit, iostat, iomsg)
    call read_bin_array(ttherm%dyngrsn,      unit, iostat, iomsg)
    call read_bin_array(ttherm%dyngra,       unit, iostat, iomsg)
end subroutine READ_T_ICE_THERMO
!
!
!_______________________________________________________________________________
! Unformatted writing for T_ICE_ATMCOUPL
#if defined (__cpl_enabled)
subroutine WRITE_T_ICE_ATMCOUPL(tcoupl, unit)
    IMPLICIT NONE
    class(T_ICE_ATMCOUPL),  intent(in)     :: tcoupl
    integer,                intent(in)     :: unit
    integer                                :: iostat
    character(len=1024)                    :: iomsg
    call write_bin_array(tcoupl%oce_flx_h,      unit, iostat, iomsg)
    call write_bin_array(tcoupl%ice_flx_h,      unit, iostat, iomsg)
    call write_bin_array(tcoupl%tmpoce_flx_h,   unit, iostat, iomsg)
    call write_bin_array(tcoupl%tmpice_flx_h,   unit, iostat, iomsg)
! TODO: still compile-time, so the bin-restart stream layout is unchanged.
! Making this follow ist_itracer_idx means writing a leading flag/count and
! branching on it when reading -- a restart format change, handled
! separately.
#if defined (__cpl_direct) || defined (__cpl_oasis50)
    call write_bin_array(tcoupl%ice_alb,        unit, iostat, iomsg)
    call write_bin_array(tcoupl%enthalpyoffuse, unit, iostat, iomsg)
    call write_bin_array(tcoupl%runoff_liquid, unit, iostat, iomsg)
    call write_bin_array(tcoupl%runoff_solid, unit, iostat, iomsg)
#endif /* IFS-family partners */

end subroutine WRITE_T_ICE_ATMCOUPL  
#endif /* (__cpl_enabled) */

! Unformatted reading for T_ICE_ATMCOUPL
#if defined (__cpl_enabled)
subroutine READ_T_ICE_ATMCOUPL(tcoupl, unit)
    IMPLICIT NONE
    class(T_ICE_ATMCOUPL),  intent(inout)  :: tcoupl
    integer,                intent(in)     :: unit
    integer                                :: iostat
    character(len=1024)                    :: iomsg
    call read_bin_array(tcoupl%oce_flx_h, unit, iostat, iomsg)
    call read_bin_array(tcoupl%ice_flx_h, unit, iostat, iomsg)
    call read_bin_array(tcoupl%tmpoce_flx_h, unit, iostat, iomsg)
    call read_bin_array(tcoupl%tmpice_flx_h, unit, iostat, iomsg)
! TODO: still compile-time, so the bin-restart stream layout is unchanged.
! Making this follow ist_itracer_idx means writing a leading flag/count and
! branching on it when reading -- a restart format change, handled
! separately.
#if defined (__cpl_direct) || defined (__cpl_oasis50)
    call read_bin_array(tcoupl%ice_alb, unit, iostat, iomsg)
    call read_bin_array(tcoupl%enthalpyoffuse, unit, iostat, iomsg)
    call read_bin_array(tcoupl%runoff_liquid, unit, iostat, iomsg)
    call read_bin_array(tcoupl%runoff_solid, unit, iostat, iomsg)
#endif /* IFS-family partners */
end subroutine READ_T_ICE_ATMCOUPL
#endif /* (__cpl_enabled) */
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
       call ice%data(i)%WRITE_T_ICE_DATA(unit)
    end do
    !___________________________________________________________________________
    call ice%thermo%WRITE_T_ICE_THERMO(unit)
    call ice%work%WRITE_T_ICE_WORK(unit)
#if defined (__cpl_enabled)
    call ice%atmcoupl%WRITE_T_ICE_ATMCOUPL(unit)
#endif /* (__cpl_enabled) */

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
    ! The tracer layout is set by set_ice_tracer_layout before the restart is
    ! read, so a stream written under a different layout does not fit the
    ! allocated ice%data. Reading it would run past the end of the array.
    read(unit, iostat=iostat, iomsg=iomsg) ice%num_itracers
    if (.not. allocated(ice%data)) then
        allocate(ice%data(ice%num_itracers))
    else if (size(ice%data) /= ice%num_itracers) then
        write(iomsg, '(a,i0,a,i0,a)') &
            'ice restart holds ', ice%num_itracers, ' tracers, this run has ',&
            size(ice%data), '. Restart and run must use the same ice '// &
            'tracer layout; see set_ice_tracer_layout.'
        iostat = 1
        return
    end if
    do i=1, ice%num_itracers
       call ice%data(i)%READ_T_ICE_DATA(unit)
    end do
    !___________________________________________________________________________
    call ice%thermo%READ_T_ICE_THERMO(unit)
    call ice%work%READ_T_ICE_WORK(unit)
#if defined (__cpl_enabled)
    call ice%atmcoupl%READ_T_ICE_ATMCOUPL(unit)
#endif /* (__cpl_enabled) */

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
