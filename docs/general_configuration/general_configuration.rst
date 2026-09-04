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
- **ib_num=1** number of icebergs read from the seeding files described in :ref:`the iceberg module section<chap_general_configuration_icebergs>`. The shipped configuration pairs ``ib_num=1`` with ``use_icebergs=.false.``, so the single entry stays inactive until the module is switched on.
- **steps_per_ib_step=8** number of ocean time steps per iceberg time step; iceberg dynamics and thermodynamics are evaluated with the longer time step ``dt*steps_per_ib_step``.
- **l_allowgrounding=2** grounding mode: ``0`` free drift (grounding ignored), ``1`` slow drift with reduced velocity, ``2`` stationary (see :ref:`the iceberg module section<chap_general_configuration_icebergs>`).
- **l_cap_ibhf_n=.false.** cap the iceberg heat flux applied to interior ocean cells at a safe temperature floor per cell and time step.
- **use_icb_iron=.false.**, **icb_iron_const=50.0e-6**, **l_icb_iron_file=.false.** passive iron tracer carried by icebergs: melting releases iron in proportion to the meltwater flux, using either the constant concentration (mol per cubic metre) or per-iceberg values from ``icb_iron.dat``. The resulting flux field is diagnostic and does not feed back on the ocean.
- **ib_async_mode=0**, **thread_support_level_required=3** OpenMP-assisted asynchronous iceberg computation: ``0`` keeps iceberg calculations synchronous (reference results), ``1`` overlaps them with the ocean and sea-ice computation, ``2`` keeps the OpenMP code active but serialized for testing. The thread support level requests ``MPI_THREAD_SERIALIZED`` (``2``) or ``MPI_THREAD_MULTIPLE`` (``3``) from the MPI library.


.. _chap_general_configuration_icebergs:

The iceberg module
==================

With ``use_icebergs=.true.`` FESOM carries a set of Lagrangian icebergs that
drift over the mesh and melt, returning freshwater and latent heat to the
ocean. The icebergs are point objects with a prescribed length, width and
height. They do not occupy grid cells and do not block the flow, so the
module is a source of meltwater and heat rather than a geometric obstacle.

Each iceberg is advanced with its own momentum balance. The forces are the
drag exerted by the ocean, the atmosphere and the sea ice, the radiation
force of waves against the iceberg side, the Coriolis force and the sea
surface slope. The wave force follows Bigg et al. (1997) and uses a wave
amplitude estimated from the wind speed; it can be switched off per iceberg.
The slope term is the gradient of the sea surface height and carries its own
scaling factor. The Coriolis term is treated semi-implicitly and the water
drag implicitly, which matters because an iceberg step spans several ocean
steps.

Four melt terms are computed: melting at the base, lateral melting driven by
buoyant convection, wave erosion, and a lateral melt rate that uses the
three-equation formulation with the local temperature, salinity and velocity
of the water column. The melt fluxes enter the ocean as a freshwater flux and
a latent heat flux, both of which can be switched off individually with
``turn_off_fw`` and ``turn_off_hf`` when isolating one effect. The heat flux
is distributed over the depth range the iceberg occupies rather than applied
at the surface alone.

Icebergs are advanced only every ``steps_per_ib_step`` ocean time steps, so
with the shipped value of 8 the iceberg module runs once per eight ocean
steps.

Seeding files
"""""""""""""

The initial iceberg population is read from plain text files in the run
directory. Each file holds one value per line and is read ``ib_num`` times,
so every file must have at least ``ib_num`` lines:

- ``icb_longitude.dat`` and ``icb_latitude.dat``, the release position in degrees.
- ``icb_length.dat`` and ``icb_height.dat``, the iceberg dimensions in metres. The width defaults to the length file, so square icebergs are obtained by supplying only these two files.
- ``icb_scaling.dat``, a scaling factor that lets one model iceberg represent many real ones.
- ``icb_calving_day.dat``, the day of the year on which each iceberg enters the simulation, counted in days since the start of the run. A value of 2.5 releases the iceberg at noon on the third day. On restart the calving day is reduced by the number of days already simulated, so an unmodified file continues to work.
- ``icb_iron.dat``, only read when ``l_icb_iron_file=.true.``, giving the iron concentration carried by each iceberg.

``ib_num`` must match the number of lines in these files. The shipped
``namelist.config`` pairs ``ib_num=1`` with ``use_icebergs=.false.``, so the
single placeholder entry stays inactive until the module is switched on.
Running with an empty iceberg population is supported and simply produces no
iceberg fluxes. A realistic circum-Antarctic
distribution of nearly 7000 icebergs, derived from synthetic aperture radar
imagery, is available from :cite:`wesche2015`.

Grounding
"""""""""

An iceberg whose scaled draft exceeds the local water depth is flagged as
grounded. What happens then is controlled by ``l_allowgrounding``. With
``0`` grounding is ignored and icebergs drift freely over shallow
topography. With ``1`` a grounded iceberg keeps drifting at a reduced
velocity, and with ``2``, the default, it is held stationary until it has
melted enough to float again.

Restarts and coupled runs
"""""""""""""""""""""""""

The per-iceberg state is written to and read from its own restart file,
separate from the ocean restart. It is a text file with one record per
iceberg holding the geometry, the position, the velocity, the drag
coefficients and densities, the calving day and the grounded and melted
flags. When
``use_icesheet_coupling=.true.`` the restart also carries the state needed to
exchange calving with an ice-sheet model. Setting ``lverbose_icb=.true.``
prints per-iceberg diagnostics.
