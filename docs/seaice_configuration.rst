.. _chap_seaice_configuration:

Standard sea ice configuration (namelist.ice)
*********************************************

Sections of the namelist
========================

Section &ice_dyn
""""""""""""""""

- **whichEVP=0** Version of EVP. ``0`` - standard EVP, ``1`` - modified EVP (as in Kimmritz et al. (2015)). For this option don't forget to adjust alpha and beta in the code, ``2`` - adaptive EVP (as in Kimmritz et al., 2016). See Koldunov et al., 2019 for hints on tuning.
- **Pstar=30000.0**
- **delta_min=1.0e-11**
- **evp_rheol_steps=150** Number of EVP substeps. See Koldunov et al., 2019 for hints on tuning.
- **Cd_oce_ice=0.0055** Drag coefficient between ocean and ice.
- **ice_gamma_fct=0.5**
- **ice_diff=0.0**
- **theta_io=0.0**
- **ice_ave_steps=1** Ice step=ice_ave_steps*oce_step


Section &ice_therm
""""""""""""""""""

- **Sice=4.0**
- **h0=.5**
- **emiss_ice=0.97**
- **emiss_wat=0.97**
- **albsn=0.81**
- **albsnm=0.77**
- **albi=0.7**
- **albim=0.68**
- **albw=0.1**

