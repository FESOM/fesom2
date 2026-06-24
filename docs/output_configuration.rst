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
