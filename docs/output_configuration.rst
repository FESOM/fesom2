.. _chap_output_configuration:

Output configuration (namelist.io)
**********************************

``namelist.io`` controls which diagnostics are computed, how often they are written, and how output files are split. The namelist contains two parts: diagnostic flags in ``&diag_list`` and the output stream definitions in ``&nml_general``/``&nml_list``.

Sections of the namelist
========================

Section &diag_list
""""""""""""""""""

- **ldiag_solver** write solver convergence diagnostics (``ssh_rhs``, iteration counts).
- **lcurt_stress_surf** enable curl of the surface stress field (activates ``curl_surf`` output).
- **ldiag_curl_vel3** compute 3-D relative vorticity from velocity (``curl_u``).
- **ldiag_Ri** compute and output bulk Richardson number diagnostics (``shear``, ``Ri``).
- **ldiag_TurbFlux** turbulent flux diagnostics (``KvdTdz``, ``KvdSdz``).
- **ldiag_salt3D** extra 3-D salinity diagnostics.
- **ldiag_dMOC** density-space MOC diagnostics (velocity transports in density bins).
- **ldiag_DVD** Discrete Variance Decay diagnostics for numerical mixing estimates.
- **ldiag_forc** comprehensive forcing diagnostics bundle (``FORC`` stream).
- **ldiag_extflds** extended field diagnostics (additional helper fields used by some analyses).
- **ldiag_destine** vertically integrated heat content diagnostics (``hc300m``, ``hc700m``, ``hc``).
- **ldiag_trflx** tracer flux diagnostics (``utemp``, ``vtemp``, ``usalt``, ``vsalt``).
- **ldiag_ice** extra ice-volume diagnostics (``vol_ice``, ``vol_snow``).
- **ldiag_uvw_sqr** squared velocity diagnostics (``u2``, ``v2``, ``w2``).
- **ldiag_trgrd_xyz** horizontal/vertical tracer gradients (``temp_grd*``, ``salt_grd*``).
- **ldiag_cmor** CMIP-style fields (``tos``, ``sos``, ``pbo``, ``volo``, etc.) that mirror CMOR definitions.

Section &nml_general
""""""""""""""""""""

- **io_listsize=120** allocate space for the number of streams specified in ``io_list``. The code trims the list if it encounters ``'unknown   '`` IDs, so it is safe to set a generous value.
- **vec_autorotate=.false.** automatically rotate vector fields to geographic coordinates on output when a rotated grid is used.
- **lnextGEMS=.false.**, **nlev_upper=0** optionally write a reduced-column output tailored for nextGEMS (temperature/salinity/velocity in the upper ``nlev_upper`` levels at 3-hourly resolution).
- **filesplit_freq='y'** output splitting: ``'y'`` writes one file per year, ``'m'`` splits monthly files.
- **compression_level=1** netCDF compression level (1–9).

Section &nml_list
"""""""""""""""""

- **io_list** enumerates output variables as tuples ``'id', freq, unit, precision``. ``freq`` is an integer, ``unit`` one of ``y/m/d/h/s`` for years/months/days/hours/steps, and ``precision`` is 4 (single) or 8 (double). The ``config/namelist.io`` file contains a full catalog of supported IDs at the bottom; copy the desired lines into this section to enable them.
- Some IDs require specific flags in ``&diag_list`` or other namelists: e.g. ``curl_surf`` needs ``lcurt_stress_surf=.true.``, CMOR variables require ``ldiag_cmor=.true.``, tracer gradient outputs need ``ldiag_trgrd_xyz=.true.``, and melt-pond outputs require ``use_meltponds`` in ``namelist.ice``.
- When ``use_ice`` is active, ice-specific fields (``uice``, ``vice``, ``a_ice``, ``m_ice``, ``m_snow``, etc.) can be requested like any other stream. When Icepack is used, additional outputs are configured via ``&nml_list_icepack`` in ``namelist.icepack``.

Icepack output (namelist.icepack)
""""""""""""""""""""""""""""""""""

- **&nml_general/io_listsize** in ``namelist.icepack`` mirrors ``namelist.io`` but is read by the Icepack driver. Set ``io_list_icepack`` to the desired Icepack diagnostics (see the comments in ``config/namelist.icepack``). Icepack output is written separately from the ocean/sea-ice streams defined above.

.. _sec_parallel_io:

Parallel output and restart I/O
===============================

By default every output or restart field is gathered onto a single rank, which
then writes it with serial netCDF. On large meshes with many MPI ranks this
gather becomes the bottleneck: one rank must hold each global field and all
other ranks wait for it. The ``&io_parallel`` group in ``namelist.config`` (not
``namelist.io``) switches to a collective mode in which a subset of ranks, the
writers, writes its own contiguous slice of a shared netCDF file through
MPI-IO.

Because a rank owns a scattered set of global indices, the data is first
redistributed (one ``MPI_Alltoallv``) into uniform contiguous blocks, one per
writer. Each writer then issues identical collective ``put_var`` calls, and the
netCDF chunk along the horizontal dimension equals the writer block, so
compressed parallel writing needs no chunk redistribution inside HDF5. The same
machinery serves restart writing and restart reading, with separately tunable
rank counts.

The writers are a subset of the compute ranks, spread over the machine; no
extra ranks need to be reserved. A 3-D field moves in a single exchange, since
FESOM stores a node's whole column contiguously, and each writer then writes
its block level by level.

Section &io_parallel
""""""""""""""""""""

- **parallel_write=.false.** enable collective parallel output. Switch it on only at a file boundary (a new year, or a new month with ``filesplit_freq='m'``): the collective path cannot append to a file that the serial gather path created.
- **n_writers=0** number of output writer ranks. ``0`` uses as many as the minimum-block-size guard allows; a larger request is silently clamped.
- **n_writers_restart=-1**, **n_readers_restart=-1** writer and reader counts for restart files, ``-1`` meaning the same as ``n_writers``. The three optima differ in practice, and reading is typically fastest with far fewer ranks than writing, hence the separate knobs.
- **chunk_levels=8** vertical levels per netCDF chunk in 3-D output. Lower values favour reading horizontal maps, higher values favour reading vertical profiles.

.. attention::
   ``parallel_write=.true.`` requires a netCDF library built against parallel
   HDF5 (MPI support). With a serial netCDF build the collective
   ``nf90_create_par`` call fails at runtime. The serial gather default works
   with any netCDF build.

The files produced by both paths are identical in layout (netCDF-4 classic
model), so post-processing does not need to know which path wrote them.

.. planned figure: schematic of compute ranks feeding a writer subset via redistribution, writers holding contiguous blocks of one shared netCDF file.

.. _sec_xios_output:

Output via XIOS
===============

FESOM can hand its output to the `XIOS <https://forge.ipsl.jussieu.fr/ioserver>`_
I/O server instead of writing files itself; consult the XIOS upstream
documentation for anything beyond the FESOM-specific setup described here. This
path is intended for coupled setups (AWI-ESM3) in which OASIS has already split
the world communicator and dedicated ``xios_server.exe`` ranks run alongside
the model components.

Build with

::

    export XIOS_ROOT=/path/to/xios/installation
    cmake -DFESOM_WITH_XIOS=ON ...

``XIOS_ROOT`` must contain ``xios.mod`` under ``inc`` or ``include`` and the
XIOS library under ``lib``. The build turns ``ASYNCHRONOUS_IO_THREADS`` off,
since XIOS needs a deterministic field send order across ranks.

When XIOS is active, the stream definitions in ``namelist.io`` are ignored:
every known stream id is registered internally, and the XIOS XML alone decides
which fields are written, at what cadence, and with what time averaging. The
XML files live in ``docs/xios_xml/`` and are copied into the run directory:

- **iodef.xml** the top-level simulation file: one context per model component plus an ``xios`` context with server and buffer tuning.
- **context_fesom.xml** the FESOM context. It sets the Gregorian calendar, includes the definition files below, and can override ``&diag_list`` flags from XML (a flag omitted there keeps its namelist value).
- **field_def_fesom.xml** declares every field id with units, standard name, averaging operation and missing-value handling.
- **file_def_fesom.xml** maps fields to files and output frequencies. One file per variable, named ``<var>.fesom`` plus the XIOS split suffix, split yearly.
- **domain_def_fesom.xml** the unstructured node and element domains; sizes and coordinates are filled at runtime from the mesh.
- **grid_def_fesom.xml** composes the domains with vertical axes into 2-D and 3-D grids.
- **axis_def_fesom.xml** the vertical axes ``nz`` (layer mid-points) and ``nz1`` (layer interfaces), plus the density axis used by the dMOC diagnostic.

To add an output, pick a field id from ``field_def_fesom.xml`` (or declare a
new one there) and add a ``<file>`` entry referencing it to a file group in
``file_def_fesom.xml`` with the desired ``output_freq``; a whole group can be
switched off with ``enabled="false"``. Sea-ice fields are declared with
``detect_missing_value="true"``: FESOM sends fill values where there is no ice,
and XIOS excludes those samples from the time mean. Dry cells below the bottom
are masked the same way.

Track and transect output
"""""""""""""""""""""""""

An XIOS build can additionally write curtain output along user-supplied
polylines (ship tracks, mooring arrays). Each polyline CSV is intersected with
the mesh once at startup: node-based variables such as ``temp`` are linearly
interpolated at the mesh edges the line crosses, while element-based variables
such as ``u`` are sampled on the triangles along the crossing path.

Activation needs ``ltracks=.true.`` in ``&nml_general`` of ``namelist.io``
(default ``.false.``) together with a non-empty ``track_files``; in XIOS mode
the same variables can instead be set in the FESOM context XML, which then
overrides the namelist. The ``&nml_tracks`` group in ``namelist.io`` holds the
details:

- **track_files** semicolon-separated list of polyline CSV files, one per track.
- **track_vars='temp'** variable to sample per track; a single entry applies to all tracks.
- **track_names** short names used in field and file ids; defaults derived from the file names.
- **track_output_freq='1h'** output frequency for all tracks, in XIOS duration syntax.

Each track and variable pair becomes its own XIOS file named
``<var>_track_<name>.fesom`` plus the split suffix, collected in a file group
that splits yearly.

.. _sec_cmor_diag:

CMOR diagnostics
================

Setting ``ldiag_cmor=.true.`` in ``&diag_list`` activates a set of fields that
follow CMOR conventions for CMIP-style output, ported from FESOM 1.4. Two kinds
of fields are computed:

- 2-D fields on the surface or sea floor: ``tos``, ``sos``, ``pbo`` and the column-integrated tendency terms (``opottemptend`` and related heat and salt budget fields).
- global scalars: ``volo``, ``soga``, ``thetaoga`` and the hemispheric sea-ice scalars (``siarean``/``siareas``, ``siextentn``/``siextents``, ``sivoln``/``sivols``).

The corresponding ``io_list`` entries ship commented out at the bottom of
``config/namelist.io``; uncomment the ones you need, keeping
``ldiag_cmor=.true.``, for example

::

    'tos       ',1, 'm', 8,
    'volo      ',1, 'm', 8,

Scalar diagnostics are written by rank 0 into small per-variable time-series
files named ``<var>.<runid>.<year>.nc``. In an XIOS run the CMOR fields go
through XIOS like any other field, and ``ldiag_cmor`` can be forced from
``context_fesom.xml``.

.. _sec_output_time_axis:

Output time axis
================

Since FESOM 2.8 the time axis of mean output follows CF conventions for time
averages. Each record is stamped at the midpoint of its averaging interval, and
every file carries a ``time_bounds(time, axis_nbounds)`` variable holding the
exact start and end of the interval that was actually accumulated. The bounds
are honest for partial intervals: the first record after a cold start or a
restart into the middle of an output interval covers only the simulated part.
Mean variables carry ``cell_methods = "time: mean"``, and the ``time`` variable
carries a ``calendar`` attribute.

Appending to a file written by an older FESOM version keeps working: new
records get the new midpoint stamps, but no bounds variable is added to the old
file, and a warning is printed.

When a period is re-run, for example after restarting from an earlier restart
file, the existing record for that period is overwritten in place instead of a
duplicate being appended. The matching is done on the time stamp, so it also
works when the time step changed between the runs. XIOS output already behaves
this way, so both output paths deliver the same time axis.
