! Ice namelist
&ice_dyn
whichEVP=0
Pstar=30000.0            
delta_min=1.0e-11
evp_rheol_steps=150
Cd_oce_ice=0.0055
ice_gamma_fct=0.5
ice_diff=0.0
theta_io=0.0                  !0.436
ice_ave_steps=1 !ice step=ice_ave_steps*oce_step
/

&ice_therm
Sice=4.0
h0=.5
emiss_ice=0.97
emiss_wat=0.97
albsn=0.81
albsnm=0.77
albi=0.7
albim=0.68
albw=0.1
/
