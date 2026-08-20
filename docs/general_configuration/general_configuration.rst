.. _chap_general_configuration:

General configuration (namelist.config)
***************************************

``namelist.config`` contains the switches that define the overall experiment setup: timing, restart cadence, paths to meshes and input data, vertical coordinate choice, grid geometry and which optional model components are active.

Sections of the namelist
========================

Section &modelname
""""""""""""""""""

- **runid='fesom'** identifier used in output, restart and log file names. Changing it also changes the expected clock file name (``<runid>.clock``) and any downstream scripts that look for ``fesom`` by default.

Section &timestep
"""""""""""""""""

- **step_per_day=36** number of baroclinic steps per day. The code enforces that ``86400 mod step_per_day == 0`` and prints the supported values when the check fails. It sets the fundamental time step ``dt = 86400/step_per_day`` seconds that is shared by ocean and ice components. Valid values are, for example: 32(45min), 36(40min), 48(30min), 60(24min), 72(20min), 144(10min), 288(5min), 1440(1min).
- **run_length=1** length of the submitted integration segment.
- **run_length_unit='y'** unit for ``run_length``; one of ``'y'`` (years), ``'m'`` (months), ``'d'`` (days), or ``'s'`` (time steps). The calendar in :ref:`Section &calendar<chap_general_configuration_calendar>` determines how months/years are counted.

Section &clockinit
""""""""""""""""""

- **timenew=0.0**, **daynew=1**, **yearnew=1958** define the starting time stamp. If this triplet matches the first line of ``<runid>.clock`` the model performs a cold start from initial conditions; otherwise it expects restart files at ``RestartInPath`` and does a warm start. The driver resets ``timenew`` to ``0`` when it encounters ``86400`` to keep times in the ``[0,86400)`` range.

Section &paths
""""""""""""""

- **MeshPath** location of the unstructured grid description (``nod2d.out``, ``elem2d.out``, depth files, etc.).
- **ClimateDataPath** directory containing initial hydrography and optional restoring climatologies.
- **TideForcingPath** directory for external tidal potential files if tides are used.
- **ResultPath** root directory for all model output, including ``<runid>.clock``.
- **RestartInPath**, **RestartOutPath** allow separating where restarts are read from and where they are written. If they are left empty, both default to ``ResultPath``; the setup routine adds trailing slashes automatically.
- **MeshId** optional string tag written into coupler metadata (useful when multiple meshes are supported in a workflow).

Section &restart_log
""""""""""""""""""""

- **restart_length=1**, **restart_length_unit='y'** cadence for portable netCDF restarts (units: ``y``, ``m``, ``d``, ``h``, ``s``, or ``off`` to disable). When a portable restart is written, the raw/binary formats are written as well unless explicitly disabled.
- **raw_restart_length=1**, **raw_restart_length_unit='y'** frequency of raw binary “core dump” restarts. If the unit is ``off`` only the portable restarts are produced.
- **bin_restart_length=1**, **bin_restart_length_unit='y'** frequency of derived-type binary restarts (Fortran binary) used for fast restarts on identical hardware/compiler stacks.
- **logfile_outfreq=960** number of model steps between status lines written to the stdout/stderr log.

Section &ale_def
""""""""""""""""

- **which_ALE='zstar'** choice of vertical coordinate: ``'linfs'`` (linear free surface), ``'zlevel'`` (fixed levels with moving surface layer only), or ``'zstar'`` (all layers except the bottom move with sea surface height). ``zstar`` keeps column thickness proportional to local depth and is the recommended default.
- **use_partial_cell=.false.** enables bottom partial cells. When enabled, bottom layer thickness can deviate from the canonical level spacing to better fit bathymetry.
- **partial_cell_thresh=0.0** minimum full-cell thickness (in metres) allowed before applying a thinner partial bottom cell. Prevents already thin layers from being reduced further.
- **min_hnode=0.5**, **lzstar_lev=4** controls the ``zlevel`` → local ``zstar`` fallback: if a surface layer would shrink below ``min_hnode`` of its nominal thickness, the surface height anomaly is distributed over ``lzstar_lev`` levels instead of one.
- **max_ice_loading=5.0** cap on how much overburden from floating ice is transferred to the ocean when ``use_floatice`` is active; excess is discarded.

Section &geometry
"""""""""""""""""

- **cartesian=.false.** use a planar Cartesian grid instead of spherical geometry (useful for idealised setups).
- **fplane=.false.** use a constant Coriolis parameter instead of latitude-dependent rotation.
- **cyclic_length=360** longitudinal domain length in degrees; change for limited-area or regional cyclic domains.
- **rotated_grid=.true.**, **force_rotation=.true.** control whether grid coordinates are rotated away from the geographic pole. ``force_rotation`` enforces rotation even if the mesh was already rotated (used for coupled runs to avoid a land pole).
- **alphaEuler=50.**, **betaEuler=15.**, **gammaEuler=-90.** Euler angles (degrees) describing the rotation applied when ``rotated_grid`` is used; converted to radians internally.
- **which_depth_n2e='mean'** how nodal depths are mapped to elements when no element depth file is supplied: ``'mean'`` (default), ``'min'`` (use shallowest node), ``'max'`` (use deepest node).
- **use_depthonelem=.false.**, **use_depthfile='aux3d'** control which bathymetry file is read: ``depth@elem.out`` (element-based) when ``use_depthonelem`` is true, otherwise ``depth@node.out``. ``use_depthfile`` switches between reading auxiliary depths from ``aux3d.out`` or the ``depth@`` files.
- **use_cavityonelem=.false.** analogous flag for cavity bathymetry (``cavity_depth@elem.out``).
- **metric_factor_zero=.false.** set the metric factor to zero, effectively switching to Cartesian metric factors even on spherical meshes (expert/debug option).

.. _chap_general_configuration_calendar:

Section &calendar
"""""""""""""""""

- **include_fleapyear=.true.** include leap years in the model clock. Set to ``.false.`` for perpetual-365-day calendars (e.g. CORE forcing).
- **use_flpyrcheck=.true.** aborts early if the leap-year setting is inconsistent with the calendar encoded in forcing files (helpful when swapping between CORE/JRA/ERA forcings).

Section &run_config
"""""""""""""""""""

- **use_ice=.true.** enable the dynamic/thermodynamic sea ice model.
- **use_floatice=.false.** allow floating ice above the ocean (ice shelf/iceberg loading); requires ``zlevel`` or ``zstar`` ALE.
- **use_sw_pene=.true.** apply shortwave penetration below the surface; requires chlorophyll data or a constant value (see ``namelist.forcing``).
- **use_cavity=.false.**, **use_cavity_partial_cell=.false.** activate ice-shelf cavities and optional surface partial cells inside cavities.
- **cavity_partial_cell_thresh=0.0** minimum cavity layer thickness before applying surface partial cells (prevents extremely thin top cavity layers).
- **use_cavity_fw2press=.true.** whether freshwater fluxes in cavities affect the hydrostatic pressure field.
- **toy_ocean=.false.**, **which_toy='soufflet'** enable idealised “toy” forcing/geometry setups.
- **flag_debug=.false.**, **flag_warn_cflz=.true.** runtime verbosity and vertical CFL warnings.
- **lwiso=.false.** enable water isotope tracers (adds isotope tracers internally).
- **use_transit=.false.** enable transient tracer package (CFCs, SF6, etc.; controlled via ``namelist.transit``).
- **compute_oasis_corners=.false.** compute grid cell corners for conservative coupling through OASIS.

Section &machine
""""""""""""""""

- **n_levels=1**, **n_part=...** hierarchical partitioning settings for the mesh partitioner. The product of ``n_part`` entries gives the total number of MPI ranks targeted by the partition.

Section &icebergs
"""""""""""""""""

- **use_icebergs=.false.**, **use_icesheet_coupling=.false.** toggle the iceberg module and coupling to an ice-sheet model.
- **turn_off_hf=.false.**, **turn_off_fw=.false.** disable latent heat or freshwater fluxes from icebergs when needed for debugging.
- **lbalance_fw=.true.**, **cell_saturation=2** controls for preventing excessive freshwater injection into small grid cells.
- **lmin_latent_hf=.true.**, **lverbose_icb=.false.** control numerical safety and verbosity of iceberg thermodynamics.
- **ib_num=0**, **steps_per_ib_step=8** number of iceberg classes and sub-cycling of iceberg dynamics relative to the ocean step.
- **ib_async_mode=0**, **thread_support_level_required=3** OpenMP-assisted asynchronous iceberg computation; the default keeps iceberg calculations synchronous.


