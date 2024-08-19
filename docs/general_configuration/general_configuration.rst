.. _chap_general_configuration:

General configuration (namelist.config)
***************************************

General configuation is defined in the ``namelist.conf``. Here you define time stepping and restart frequency, details of the ALE and mesh geometry.

Sections of the namelist
========================

Section &modelname
""""""""""""""""""

- **runid='fesom'** define name of the run. It will be used as part of the file name in model output and restarts. Don't change it if you don't have a good reason, since many post processing scripts assume it to be ``fesom``.

Section &timestep
"""""""""""""""""

- **step_per_day=32** define how many steps per day the model will have. The variable ``step_per_day`` must be an integer multiple of 86400 ``(mod(86400,step_per_day)==0)``. Valid values are, for example: 32(45min), 36(40min), 48(30min), 60(24min), 72(20min), 144(10min), 288(5min), 1440(1min).
- **run_length= 62** length of the run in ``run_length_unit``.
- **run_length_unit='y'** units of the ``run_length``. Valid values are year (``y``), month (``m``), day (``d``), and model time step (``s``).

Section &clockinit
""""""""""""""""""

- **timenew=0.0**, **daynew=1**, **yearnew=1948** gives the seconds, day and year of the initialisation time point, respectively. If the initialisation time is identical with the first line of the clock file runid.clock the model performs a cold start. If the initialisation time and the first line of the clock file are not identical the model assumes that a restart file must exist and tries to do a warm start.


Section &paths
""""""""""""""

- **Meshpath=''**, path to your mesh directory
- **ClimateDataPath=''**, path to the location of your 3D initialisation data for temperatur and salinity. 
- **Resultpath=''**, directory where your results should be stored


Section &restart_log
""""""""""""""""""""

- **restart_length=1**, how often should restart file be written in units of  ``restart_length_unit``
- **restart_length_unit='y'** units of the ``restart_length``. Valid values are year (``y``), month (``m``), day (``d``), and model time step (``s``).
- **logfile_outfreq=960** the frequency (in timesteps), the model state information should be written into the job monitor .log/.out file.


Section &ale_def
""""""""""""""""

- **which_ALE='linfs'**, which Arbitrary Lagrangian Eulerian (ALE) approach should be used? Options are 1) ``linfs`` - vertical grid is fixed in time, 2) ``zlevel`` - only the surface layer is allowed to move with the change in ssh all other levels are fixed in time 3) ``zstar`` - all layers, except the bottom layer are allowed to move, the change in ssh is equally distributed over all layers. It is recommended to use either ``linfs`` or ``zstar``.
- **use_partial_cell=.false.**, switch if partial bottom cells should be used or not. Partial cell means that the bottom layer thickness can be different from the full depth levels to be closer to the real bottom topography
- **min_hnode=0.5**, for ``zlevel``: layer thickness should not become smaller than min_hnode [in fraction from 0 to 1] of original layer thickness. If it happens switch from ``zlevel`` to local ``zstar``.
- **lzstar_lev=4**, for ``zlevel``  in case min_hnode criteria is reached over how many level should ssh change be distributed for local zstar
- **max_ice_loading=5.0**, for ``use_floatice=.True.`` in case of floating sea ice how much ice loading is allowed [unit m] the excess is discarded

Section &geometry
"""""""""""""""""

- **cartesian     =.false.**, use flat cartesian coordinates (idealized geometry)
- **fplane        =.false.**, use fplane approximation, coriolis force is lat independent coriolis=2*omega*0.71
- **rotated_grid  =.true.**, should the model perform on rotated grid.
- **force_rotation=.false.**, if input mesh is unrotated it must be rotated in FESOM2.0 than ``.true.``, if input mesh is already rotated ``.false.``

- **alphaEuler    =50.0**, rotated Euler angles, alpha [degree]
- **betaEuler     =15.0**, rotated Euler angles, beta [degree]
- **gammaEuler    =-90.0**, rotated Euler angles, gamma [degree]


Section &calendar
"""""""""""""""""

- **include_fleapyear=.false.**, should be ``.true.`` when the forcing contains fleapyears (i.e. NCEP...)


Section &run_config
"""""""""""""""""""

- **use_ice     =.true.**, simulate ocean + sea ice
- **use_floatice=.false.**, allow floating sea ice only possible with ``zlevel`` or ``zstar``
- **use_sw_pene =.true.**, use parameterisation for short wave penetration. Incoming short wave radiation is distributed over several layers


Section &machine
""""""""""""""""

- **n_levels = 2**, number of hierarchy level for mesh partitioning
- **n_part   = 12, 36**, number of partitions on each hierarchy level, the last number should optimal corresponds with the number of cores per computational node




