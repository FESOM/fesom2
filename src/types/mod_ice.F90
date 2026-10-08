module mod_ice
    ! Sea-ice type (decision D5: EVOLVING). M3 (sea ice EVP). Faithful transcription
    ! of FESOM2 v2.7.3 MOD_ICE.F90 (the standard EVP path: whichEVP=0; NO icepack /
    ! meltponds / cavity / oasis-yac / oifs). The advected ice tracers (a_ice / m_ice /
    ! m_snow) live in data(1:3)%values exactly like FESOM2 — the ice FCT advection (M3c)
    ! loops over them generically, so the data(:) layout makes that port byte-identical
    ! (mirrors FESOM3's own t_tracer). The EVP stress tensor (sigma11/12/22, on elements)
    ! is PROGNOSTIC — elastic memory persists across steps.
    !
    ! PRECISION DISCIPLINE (L: M2.10b verbatim literals): the parameter initialisers are
    ! copied CHARACTER-FOR-CHARACTER from MOD_ICE.F90, including the UN-SUFFIXED ones
    ! (e.g. rhoair=1.3, cd_oce_ice=5.5e-3, delta_min=1.0e-11). An un-suffixed literal is
    ! default real (single) converted to WP at declaration — `1.3_WP` would differ ~1 ULP.
    ! Keep them verbatim so the EVP/thermo arithmetic (M3b/M3d) byte-matches FESOM2.
    use mod_precision, only: WP
    implicit none
    save
    private
    public :: t_ice, t_ice_data, t_ice_work, t_ice_thermo

    !___________________________________________________________________________
    ! ice-tracer data: 1=a_ice (concentration), 2=m_ice (ice mass = h*a),
    ! 3=m_snow (snow mass). Mirrors FESOM2 T_ICE_DATA.
    type t_ice_data
        real(kind=WP), allocatable, dimension(:) :: values, values_old, values_rhs
        real(kind=WP), allocatable, dimension(:) :: values_div_rhs, dvalues, valuesl
        integer :: ID = 0
    end type t_ice_data

    !___________________________________________________________________________
    ! ice work arrays: FCT (on nodes/elements) + EVP stress/strain (on elements).
    type t_ice_work
        real(kind=WP), allocatable, dimension(:)   :: fct_tmax, fct_tmin
        real(kind=WP), allocatable, dimension(:)   :: fct_plus, fct_minus
        real(kind=WP), allocatable, dimension(:,:) :: fct_fluxes
        real(kind=WP), allocatable, dimension(:)   :: fct_massmatrix   ! sized ssh_stiff%nza
        real(kind=WP), allocatable, dimension(:)   :: fct_lumpedmass   ! row sums of fct_massmatrix (nodes)
        real(kind=WP), allocatable, dimension(:)   :: sigma11, sigma12, sigma22  ! EVP stress (elem)
        real(kind=WP), allocatable, dimension(:)   :: eps11, eps12, eps22        ! strain rate (elem)
        real(kind=WP), allocatable, dimension(:)   :: ice_strength, inv_areamass, inv_mass
    end type t_ice_work

    !___________________________________________________________________________
    ! ice thermodynamics: work fields + the (mostly literal) thermo parameters.
    ! Initialisers copied VERBATIM from MOD_ICE.F90 (un-suffixed = single->WP).
    type t_ice_thermo
        real(kind=WP), allocatable, dimension(:) :: t_skin, thdgr, thdgrsn, thdgra, thdgr_old, ustar
        real(kind=WP), allocatable, dimension(:) :: dyngr, dyngrsn, dyngra
        !_______________________________________________________________________
        real(kind=WP) :: rhoair=1.3  , inv_rhoair=1./1.3    ! Air density & inverse,  LY2004
        real(kind=WP) :: rhowat=1025., inv_rhowat=1./1025.  ! Water density & inverse
        real(kind=WP) :: rhofwt=1000., inv_rhofwt=1./1000.  ! Freshwater density & inverse
        real(kind=WP) :: rhoice=910. , inv_rhoice=1./910.   ! Ice density & inverse, AOMIP
        real(kind=WP) :: rhosno=290. , inv_rhosno=1./290.   ! Snow density & inverse, AOMIP
        real(kind=WP) :: cpair=1005., cpice=2106., cpsno=2090.  ! Specific heats [J/(kg K)]
        real(kind=WP) :: cc=1025.*4190.0   ! Volumetr. heat cap. of water [J/m**3/K]
        real(kind=WP) :: cl=910.*3.34e5    ! Volumetr. latent heat of ice fusion [J/m**3]
        real(kind=WP) :: clhw=2.501e6      ! Specific latent heat [J/kg]: water  -> vapor
        real(kind=WP) :: clhi=2.835e6      !                              sea ice -> vapor
        real(kind=WP) :: tmelt=273.15      ! 0 deg C expressed in K
        real(kind=WP) :: boltzmann=5.67E-8 ! S. Boltzmann const.*longw. emissivity
        integer       :: iclasses=7        ! ice thickness gradations for growth calcs
        real(kind=WP) :: hmin= 0.01        ! Cut-off ice thickness
        real(kind=WP) :: Armin=0.01        ! Minimum ice concentration
        ! --- namelist /ice_therm/ ---
        real(kind=WP) :: con= 2.1656, consn = 0.31 ! Thermal conductivities: ice & snow; W/m/K
        real(kind=WP) :: Sice = 4.0        ! Ice salinity 3.2--5.0 ppt
        real(kind=WP) :: h0=0.5            ! Lead closing parameter [m] (NH)
        real(kind=WP) :: h0_s=0.5          ! Lead closing parameter [m] (SH)
        real(kind=WP) :: emiss_ice=0.97    ! Emissivity of snow/ice
        real(kind=WP) :: emiss_wat=0.97    ! Emissivity of open water
        real(kind=WP) :: albsn = 0.81      ! Albedo: frozen snow
        real(kind=WP) :: albsnm= 0.77      !         melting snow
        real(kind=WP) :: albi  = 0.70      !         frozen ice
        real(kind=WP) :: albim = 0.68      !         melting ice
        real(kind=WP) :: albw  = 0.066     !         open water, LY2004
        real(kind=WP) :: h_ml  = 2.5_WP    ! thickness of uppermost layer (heat available)
        logical       :: snowdist=.true.
        logical       :: new_iclasses=.false.
        integer       :: open_water_albedo=0
        real(kind=WP) :: c_melt=0.5
        logical       :: use_meltponds=.false.
        ! new_iclasses (EM thickness distribution, Castro-Morales 2013) params. therm_ice
        ! pointer-associates these UNCONDITIONALLY (MOD_ICE.F90:94-97), but only reads them
        ! when new_iclasses=.true. — always .false. in the gated path, so they never affect
        ! the byte-result; kept verbatim for a faithful transcription. h_cutoff is in the
        ! &ice_therm namelist (=3.0, exact); hpdf is hardcoded (single->WP both sides).
        real(kind=WP) :: h_cutoff=3.0
        real(kind=WP), dimension(15) :: hpdf = (/ 0.066745491, 0.1462317, 0.17769822, 0.13131106, &
             0.11518432, 0.08514193, 0.06871303, 0.05592151, 0.04428673, 0.03584652, 0.02970195, 0.02469673, &
             0.02001543, 0.01653681, 0.0141026 /)
    end type t_ice_thermo

    !___________________________________________________________________________
    ! main ice type. Mirrors FESOM2 T_ICE (standard EVP path only).
    type t_ice
        ! zonal & meridional ice velocity (nodes)
        real(kind=WP), allocatable, dimension(:) :: uice, uice_rhs, uice_old, uice_aux
        real(kind=WP), allocatable, dimension(:) :: vice, vice_rhs, vice_old, vice_aux
        ! surface stress atm<->ice, oce<->ice (nodes)
        real(kind=WP), allocatable, dimension(:) :: stress_atmice_x, stress_iceoce_x
        real(kind=WP), allocatable, dimension(:) :: stress_atmice_y, stress_iceoce_y
        ! oce temp/salt/ssh and uv at surface (nodes)
        real(kind=WP), allocatable, dimension(:) :: srfoce_temp, srfoce_salt, srfoce_ssh
        real(kind=WP), allocatable, dimension(:) :: srfoce_u, srfoce_v
        ! freshwater & heat flux (nodes)
        real(kind=WP), allocatable, dimension(:) :: flx_fw, flx_h
        ! ice/snow thickness in the ice-covered area (nodes)
        real(kind=WP), allocatable, dimension(:) :: h_ice, h_snow
        ! node boundary mask for the mEVP velocity solve (whichEVP/=0): 1 interior,
        ! 0 at nodes touching a boundary edge. FESOM2 stores this in mesh%bc_index_nod2D
        ! (MOD_ICE.F90:889); kept in t_ice here (the FESOM3 t_mesh has type-bound I/O —
        ! adding a component there trips an ifort generic-serialization cascade). Built in
        ! ice_setup; 0/1 so WP-exact. Value byte-identical to FESOM2's mesh field.
        real(kind=WP), allocatable, dimension(:) :: bc_index_nod2D
        !_______________________________________________________________________
        ! ice tracers (1=a_ice, 2=m_ice, 3=m_snow)
        integer :: num_itracers = 3
        type(t_ice_data), allocatable, dimension(:) :: data
        type(t_ice_work)   :: work
        type(t_ice_thermo) :: thermo
        !_______________________________________________________________________
        ! ice model parameters (VERBATIM defaults from MOD_ICE.F90).
        ! --- RHEOLOGY ---
        real(kind=WP) :: pstar      = 30000.0_WP   ! [N/m^2]
        real(kind=WP) :: ellipse    = 2.0_WP
        real(kind=WP) :: c_pressure = 20.0_WP
        real(kind=WP) :: delta_min  = 1.0e-11      ! [s^(-1)]  (un-suffixed: single->WP)
        real(kind=WP) :: Clim_evp   = 615          ! kg/m^2
        real(kind=WP) :: zeta_min   = 4.0e+8       ! kg/s
        integer       :: evp_rheol_steps=120       ! EVP subcycling steps
        real(kind=WP) :: ice_gamma_fct=0.25_WP     ! smoothing parameter in ice fct advection
        real(kind=WP) :: ice_diff   = 10.0_WP      ! diffusion to stabilise ice advection
        real(kind=WP) :: theta_io   =0.0_WP        ! rotation angle (ice-ocean)
        ! --- in EVP ---
        real(kind=WP) :: alpha_evp=250, beta_evp=250
        real(kind=WP) :: c_aevp=0.15               ! (un-suffixed: single->WP)
        ! --- Ice forcing averaging ---
        integer       :: ice_ave_steps=1           ! ice step = ice_ave_steps*oce_step
        real(kind=WP) :: cd_oce_ice = 5.5e-3       ! drag coef. oce-ice (un-suffixed)
        logical       :: ice_free_slip=.false.
        integer       :: whichEVP=0                ! 0=standard; 1=mEVP; 2=aEVP
        real(kind=WP) :: ice_dt                    ! ice step = ice_ave_steps*oce_step
        real(kind=WP) :: Tevp_inv
        integer       :: ice_steps_since_upd=0
        logical       :: ice_update = .true.
    end type t_ice

end module mod_ice
