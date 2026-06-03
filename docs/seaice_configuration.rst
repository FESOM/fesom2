.. _chap_seaice_configuration:

Standard sea ice configuration (namelist.ice)
*********************************************

Sea-ice options are split between the dynamics (&ice_dyn) and thermodynamics (&ice_therm) sections. The same namelist is used whether the standard FESOM ice thermodynamics or Icepack is active (Icepack overrides the thermodynamics when enabled).

Sections of the namelist
========================

Section &ice_dyn
""""""""""""""""

- **whichEVP=0** EVP rheology flavour: ``0`` standard EVP, ``1`` modified EVP (mEVP, Kimmritz et al. 2015; tune ``alpha_evp``/``beta_evp``), ``2`` adaptive EVP (aEVP, Kimmritz et al. 2016; tune ``c_aevp``).
- **Pstar=30000.0**, **ellipse=2.0**, **c_pressure=20.0** yield-curve parameters controlling ice strength and ellipticity.
- **delta_min=1.0e-11** minimum strain rate used to regularise viscosities.
- **evp_rheol_steps=120** number of EVP sub-cycles per ice time step.
- **Cd_oce_ice=0.0055** ocean–ice drag coefficient.
- **ice_gamma_fct=0.5**, **ice_diff=0.0**, **theta_io=0.0** numerical stabilisation: damping factor for EVP rhs, small Laplacian diffusion and optional stress rotation angle.
- **ice_ave_steps=1** number of ocean time steps per ice dynamics step (sets sub-cycling).
- **alpha_evp=250**, **beta_evp=250**, **c_aevp=0.15** stability/tuning constants for mEVP/aEVP; use higher values on fine meshes to improve convergence.

Section &ice_therm
"""""""""""""""""

- **Sice=4.0** bulk ice salinity (ppt) used in freezing-point and enthalpy calculations.
- **iclasses=7**, **new_iclasses=.false.**, **h_cutoff=3.0** number of ice thickness categories. When ``new_iclasses`` is true, an observationally derived distribution (Castro-Morales et al., 2013) with ``iclasses=15`` and cutoff ``h_cutoff`` is used instead of the default even spacing.
- **h0=0.5**, **h0_s=0.5** lead-closing parameters (Northern/Southern Hemisphere) controlling how quickly open water is refrozen.
- **hmin=0.01**, **armin=0.01** minimum ice thickness and concentration to avoid numerical underflow.
- **emiss_ice=0.97**, **emiss_wat=0.97** longwave emissivities for snow/ice and open water.
- **albsn=0.81**, **albsnm=0.77**, **albi=0.7**, **albim=0.68**, **albw=0.1**, **open_water_albedo=0** shortwave albedo parameters. ``open_water_albedo`` chooses the scheme: ``0`` fixed ``albw``, ``1`` Taylor et al., ``2`` Briegleb et al.
- **con=2.1656**, **consn=0.31** thermal conductivities of ice and snow (W/m/K).
- **snowdist=.true.** distribute snowfall across thickness categories instead of a single class.
- **use_meltponds=.false.** enable the built-in melt-pond parameterization; when active, melt-pond diagnostics can be written from ``namelist.io``.
- **c_melt=0.5** concentration decay rate during melting in the standard FESOM ice thermodynamics.
- **h_ml=0.0** optional mixed-layer depth used in some simplified melt formulations (kept at 0 unless explicitly needed).
