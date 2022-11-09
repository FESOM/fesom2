.. _chap_ocean_configuration:

Ocean configuration (namelist.oce)
**********************************

Sections of the namelist
========================

Section &oce_dyn
""""""""""""""""

- **C_d=0.0025**, Bottom drag, nondimensional.
- **A_ver= 1.e-4**, Vertical viscosity, m^2/s
- **laplacian=.false.**,  Use Laplacian viscosity
- **A_hor=0.**, Background horizontal viscosity
- **A_hor_max=0.**,  Maximum viscosity allowed (to limit Smag and Leith contributions when they are too large
- **Div_c=.5**,  Modified Leith viscosity, nondimensional, 0.3 -- 1.0
- **Leith_c=.05**, The strength of the Leith viscosity.
- **tau_c= 0.**, Controls the strength of filters (1.5 and 0.2 for dt=1min and 15min, respectively)
- **Smag_c=0.**, Smagorinsky viscosity, nondimensional, 0.1 --0.2
- **biharmonic=.false.**, Use biharmonic viscosity.
- **visc_option=5**, Option 2 is to use Laplacian+Leith+biharmonic background, Option 5 is to use easy backscatter.
- **easy_bs_scale = 35.**,  Area scaling, to be used with visc_option=5 (easy backscatter)
- **easy_bs_return= 1.5**,  Coefficient for returned sub-gridscale energy, to be used with visc_option=5 (easy backscatter)
- **Abh0=0.**, Biharmonic viscosity, m^4/s
- **scale_area=5.8e9**,  Viscosity and diffusion are for an element with ``scale_area``.
- **mom_adv=2**,  1=vector CV, p1 vel, 2=sca. CV, 3=vector inv.
- **free_slip=.false.**, Switch for the free slip.
- **i_vert_visc=.true.**
- **w_split=.false.**
- **w_exp_max=1.e-3**
- **SPP=.false.**,  Salt Plume Parameterization.
- **Fer_GM=.true.**, Switch for the GM after Ferrari et al. 2010
- **K_GM_max=3000.0**, Maximum GM thickness diffusivity (m2/s)
- **K_GM_min=2.0**, Maximum GM thickness diffusivity (m2/s)
- **K_GM_bvref=2**,  def of bvref in ferreira scaling 0=srf,1=bot mld,2=mean over mld,3=weighted mean over mld
- **K_GM_rampmax=40.0**,  Resol >K_GM_rampmax[km] GM on
- **K_GM_rampmin=30.0**, Resol <K_GM_rampmin[km] GM off, in between linear scaled down
- **scaling_Ferreira=.true.**,  GM vertical scaling after Ferreira et al.(2005) (as also implemented by Qiang in FESOM 1.4)
- **scaling_Rossby=.false.**, GM is smoothly switched off according to Rossby radius (from 1. in coarse areas to 0. where resolution reaches 2 points/Rossby radius).
- **scaling_resolution =.true.**, GM is spatially scaled with resolution; A value of K_GM corresponds then to a resolution of 100km.
- **scaling_FESOM14=.true.**, Special treatment of GM in the NH (as also implemented by Qiang in FESOM 1.4; it is zero within the boundary layer)
- **Redi=.true.**
- **visc_sh_limit=5.0e-3**, For KPP, max visc due to shear instability.
- **mix_scheme='KPP'**,  Vertical mixing scheme: KPP, PP.
- **Ricr=0.3**, Critical bulk Richardson Number.
- **concv  = 1.6**, Constant for pure convection (eqn. 23) (Large 1.5-1.6; MOM default 1.8).



Section &oce_tra
""""""""""""""""

- **diff_sh_limit=0.05**, For KPP, max diff due to shear instability.
- **Kv0_const=.false.**
- **double_diffusion=.false.**, For KPP,dd switch
- **K_ver=1.0e-5**
- **K_hor=3000.**
- **surf_relax_T=0.0**
- **surf_relax_S=1.929e-06**, Its 50m/300days 6.43e-07! m/s 10./(180.*86400.)
- **balance_salt_water =.true.**, Balance virtual-salt or freshwater flux or not
- **clim_relax=0.0** 1/s, geometrical information has to be supplied
- **ref_sss_local=.true.**
- **ref_sss=34.**
- **i_vert_diff=.true.**
- **tracer_adv =2**, 1=MUSCL, 2=MUSCL+FCT
- **num_tracers=2**, Number of all tracers
- **tracer_ID=0,1**, Their IDs (0 and 1 are reserved for temperature and salinity)

Implemented trassers (3d restoring):

- 301 - Fram Strait.
- 302 - Bering Strait
- 303 - Barents Sea Opening

To turn on, for example, the Fram Strait, one has to change the namelist in the following way::

   num_tracers=3
   tracer_ID=0,1,301


Initial conditions &oce_init3d
""""""""""""""""""""""""""""""

Initial conditions are loaded from netCDF file (or files) that contain values on regular grid. On the format of netCDF files see :ref:`chap_data_processing`.

One have to specify several parameters in the ``namelist.oce``. The section of the ``namelist.oce`` might look like this:

::

    &oce_init3d                               ! initial conditions for tracers
    n_ic3d   = 2                              ! number of tracers to initialize
    idlist   = 1, 0                           ! their IDs (0 is temperature, 1 is salinity, etc.). The reading order is defined here!
    filelist = 'phc3.0_winter.nc', 'phc3.0_winter.nc' ! list of files in ClimateDataPath to read (one file per tracer), same order as idlist
    varlist  = 'salt', 'temp'                 ! variables to read from specified files
    t_insitu = .true.                         ! if T is insitu it will be converted to potential after reading it



- **n_ic3d** - how many tracers you want to initialise. As a minimum you have to initialise temperature and salinity. They have reserved id's 0 and 1 respectivelly (see **idlist** below). In this example only two tracers (temperature and salinity) are initialised.
- **idlist** - IDs of tracers. In this variable you define the reading order. Here we on purpose change the order of temperature and salinity to demonstrate that the order can be arbitrary. Once again remember that ID 0 and 1 are reserved for temperature and salinity respectively.
- **filelist** - coma separated list of files (each in qoutation marks) that contain initial conditions (see the section below about requirements to the file format). The path to the folder with this files is defined in ``namelist.config`` (**ClimateDataPath** variable). In this case the file ``phc3.0_winter.nc`` is the same for temperature and salinity since it contains both variables.
- **varlist** - names of the variables in the netCDF files specified above. Note again the order of the variables in the example, it can be arbitrary and in this case temperature comes after salinity.
- **t_insitu** - most of climatologies are distributed with in situ temperature, while model needs potential temperature. This flag allows to do the conversion (UNESCO equation) on the fly.

