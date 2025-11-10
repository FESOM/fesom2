! Ice namelist
&ice_dyn
whichEVP=0             ! 0=standart; 1=mEVP; 2=aEVP
Pstar=30000.0          ! [N/m^2]
ellipse=2.0
c_pressure=20.0        ! ice concentration parameter used in ice strength computation
delta_min=1.0e-11      ! [s^(-1)]
evp_rheol_steps=120    ! number of EVP subcycles
alpha_evp=250          ! constant that control numerical stability of mEVP. Adjust with resolution. 
beta_evp=250           ! constant that control numerical stability of mEVP. Adjust with resolution.
c_aevp=0.15            ! a tuning constant in aEVP. Adjust with resolution.
Cd_oce_ice=0.0055      ! drag coef. oce - ice 
ice_gamma_fct=0.5      ! smoothing parameter
ice_diff=0.0           ! diffusion to stabilize
theta_io=0.0           ! rotation angle
ice_ave_steps=1        ! ice step=ice_ave_steps*oce_step
/

&ice_therm
Sice=4.0               ! Ice salinity 3.2--5.0 ppt.
iclasses=7            ! default = 7; in case of EM distribution ('new_iceclasses=.true.') must be set to 15
h0=0.5                 ! Lead closing parameter for Nothern Hemisphere [m], default 0.5
h0_s=0.5               ! Lead closing parameter [m] for Southern Hemisphere, default 0.5
hmin=0.01              ! default=0.01
armin=0.01             ! default=0.01
emiss_ice=0.97         ! Emissivity of Snow/Ice,
emiss_wat=0.97         ! Emissivity of open water
albsn=0.81             ! Albedo: frozen snow
albsnm=0.77            !         melting snow
albi=0.7               !         frozen ice
albim=0.68             !         melting ice
albw=0.1               !         open water
con=2.1656             ! Thermal conductivities: ice; W/m/K
consn=0.31             !                         snow
snowdist=.true.       ! distribution of snow depth according to ice distribution - default: .true.
new_iclasses=.false.    ! default=.false.; ice thickness distribution based on EM observations (Castro-Morales et al., JGR, 2013)
open_water_albedo=0    ! 0=default; 1=taylor; 2=briegleb
c_melt=0.5             ! constant in concentration equation for melting conditions - default=0.5
h_cutoff=3.0           ! only used for new_iclasses=.true.
use_meltponds=.false.  ! enable melt pond parameterization - default=.false.
/
