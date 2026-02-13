! ============================================================================
! ============== Namelist file for FESOM2 sea ice model =====================
! ============================================================================
! This file contains configuration for sea ice dynamics and thermodynamics:
! - EVP (Elastic-Viscous-Plastic) rheology options
! - Ice strength and deformation parameters
! - Ocean-ice drag
! - Ice thermodynamics and thickness distribution
! - Albedo parameterizations
! ============================================================================

! ============================================================================
! SEA ICE DYNAMICS
! ============================================================================
&ice_dyn
! --- EVP Rheology Options ---
whichEVP       = 1              ! EVP solver type:
                                !   0 = standard EVP
                                !   1 = modified EVP (mEVP)
                                !   2 = adaptive EVP (aEVP)

! --- Ice Strength Parameters ---
Pstar          = 30000.0        ! ice strength parameter [N/m²] (typical: 20000-30000)
ellipse        = 2.0            ! aspect ratio of yield curve ellipse (dimensionless)
c_pressure     = 20.0           ! ice concentration parameter for strength computation (dimensionless)

! --- Ice Deformation ---
delta_min      = 1.0e-11        ! minimum strain rate for viscosity regularization [s⁻¹]

! --- EVP Subcycling ---
evp_rheol_steps = 120           ! number of EVP subcycles per ice time step

! --- mEVP Stability Parameters (for whichEVP=1) ---
alpha_evp      = 500            ! mEVP stability constant (adjust with resolution)
beta_evp       = 500            ! mEVP stability constant (adjust with resolution)

! --- aEVP Tuning (for whichEVP=2) ---
c_aevp         = 0.15           ! aEVP tuning constant (adjust with resolution)

! --- Ocean-Ice Coupling ---
Cd_oce_ice     = 0.0055         ! drag coefficient between ocean and ice (dimensionless, typical: 0.0055)

! --- Numerical Stabilization ---
ice_gamma_fct  = 0.5            ! smoothing parameter for ice dynamics (0.0-1.0)
ice_diff       = 0.0            ! artificial diffusion for numerical stability [m²/s]
theta_io       = 0.0            ! rotation angle for ice-ocean stress [degrees]

! --- Time Stepping ---
ice_ave_steps  = 1              ! ice time step = ice_ave_steps × ocean time step
/

! ============================================================================
! SEA ICE THERMODYNAMICS
! ============================================================================
&ice_therm
! --- Ice Properties ---
Sice           = 4.0            ! bulk ice salinity [ppt] (typical range: 3.2-5.0)

! --- Ice Thickness Distribution ---
iclasses       = 7              ! number of ice thickness categories (default: 7)
                                ! set to 15 if using EM distribution (new_iceclasses=.true.)
new_iclasses   = .false.        ! use ice thickness distribution from EM observations
                                ! (Castro-Morales et al., JGR, 2013)
h_cutoff       = 3.0            ! thickness cutoff for new_iclasses [m]

! --- Lead Closing Parameters ---
h0             = 0.5            ! lead closing parameter for Northern Hemisphere [m]
h0_s           = 0.5            ! lead closing parameter for Southern Hemisphere [m]

! --- Minimum Thresholds ---
hmin           = 0.01           ! minimum ice thickness [m]
armin          = 0.01           ! minimum ice concentration (dimensionless)

! --- Emissivity (Longwave Radiation) ---
emiss_ice      = 0.97           ! emissivity of snow/ice surface (dimensionless, 0-1)
emiss_wat      = 0.97           ! emissivity of open water (dimensionless, 0-1)

! --- Albedo (Shortwave Radiation) ---
albsn             = 0.81           ! albedo of frozen snow (dimensionless, 0-1)
albsnm            = 0.77           ! albedo of melting snow (dimensionless, 0-1)
albi              = 0.7            ! albedo of frozen ice (dimensionless, 0-1)
albim             = 0.68           ! albedo of melting ice (dimensionless, 0-1)
albw              = 0.1            ! albedo of open water (dimensionless, 0-1)
open_water_albedo = 0           ! open water albedo scheme:
                                !   0 = default (constant albw)
                                !   1 = Taylor et al.
                                !   2 = Briegleb et al.

! --- Thermal Conductivity ---
con            = 2.1656         ! thermal conductivity of ice [W/m/K]
consn          = 0.31           ! thermal conductivity of snow [W/m/K]

! --- Snow Distribution ---
snowdist       = .true.         ! distribute snow depth according to ice thickness distribution

! --- Melting Parameters ---
c_melt         = 0.5            ! constant in concentration equation for melting conditions (0-1)
/
