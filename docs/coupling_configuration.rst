Coupling configuration (namelist.cpl)
=====================================

FESOM talks to external components through exactly **one** interface, and that
choice is a compile-time property: it decides which coupler library is linked
and which API generation the driver code is written against. *Which* external
components FESOM then talks to over that interface is a run-time choice, made
in ``namelist.cpl``.

Choosing the interface (build time)
-----------------------------------

One cmake variable selects it::

    cmake -DFESOM_COUPLING=<standalone|direct|oasis28|oasis50|yac> ...

============== ============================================ ===================
Value          Interface                                    Macro defined
============== ============================================ ===================
``standalone`` no external component                        ``__standalone``
``direct``     IFS drives FESOM in-process, no coupler      ``__cpl_direct``
``oasis28``    OASIS3-MCT 2.8 (``prism_*`` API generation)  ``__cpl_oasis28``
``oasis50``    OASIS3-MCT 5.0 (``oasis_*`` API generation)  ``__cpl_oasis50``
``yac``        YAC                                          ``__cpl_yac``
============== ============================================ ===================

Because the variable is a single scalar, "exactly one interface" is structural
rather than something that has to be validated afterwards; an unknown value is
rejected at configure time. Shared code tests the derived group macros instead
of enumerating interfaces:

- ``__cpl_enabled`` -- an external component is present (anything but
  ``standalone``).
- ``__cpl_coupler`` -- a coupler library owns MPI initialisation, the FESOM
  communicator and teardown (``oasis28``, ``oasis50``, ``yac``).
- ``__cpl_oasis`` -- the OASIS driver is compiled (either generation).

The older boolean switches ``FESOM_COUPLED``, ``OIFS_COUPLED``, ``USE_YAC`` and
``ENABLE_IFS_INTERFACE`` still work for one release: they map onto
``FESOM_COUPLING`` and emit a deprecation warning. Combinations that used to
configure and then fail at compile or link time -- ``FESOM_COUPLED`` together
with ``USE_YAC``, or ``OIFS_COUPLED`` without ``FESOM_COUPLED`` -- are now
rejected at configure time.

Choosing the partner components (run time)
------------------------------------------

``namelist.cpl`` lives in the run directory alongside ``namelist.config``. It
is only read by a coupled build, and every entry has a compiled-in default, so
run directories that predate the file keep working: the partner implied by the
compiled interface is assumed.

Section &coupling_partner
"""""""""""""""""""""""""

Exactly one atmosphere must be selected. Each partner requires the interface it
talks over, and the model aborts at startup if none or several are selected, or
if the selection does not match the compiled interface.

- **is_coupled_to_echam=.false.** ECHAM6 atmosphere; requires ``oasis28``.
- **is_coupled_to_oifs=.false.** OpenIFS atmosphere; requires ``oasis50``.
- **is_coupled_to_icon_a=.false.** ICON-A atmosphere; requires ``yac``.
- **is_coupled_to_ifs=.false.** IFS atmosphere; requires ``direct``.

The selection drives the exchanged field set (the ``cpl_send``/``cpl_recv``
name tables and their counts), whether the ECHAM-only flux-correction
machinery runs (``force_flux_consv`` and the out-of-band net-flux exchange in
``net_rec_from_atm``), and whether the ice model carries a prognostic ice
surface temperature as a fourth advected ice tracer, which the IFS-family
atmospheres expect and the others do not.

Section &coupling_oasis
"""""""""""""""""""""""

Read only for ``FESOM_COUPLING=oasis28`` or ``oasis50``.

- **cpl_comp_name='fesom'** component name registered with OASIS; must match
  the ``namcouple`` / SMIOC entry.
- **cpl_grid_name='feom'** grid name written to the OASIS grid files.
- **compute_oasis_corners=.false.** also write grid corners, which
  first-order conservative remapping (``CONSERV`` in the ``namcouple``)
  needs. This is a property of the configured remapping, not of the partner,
  which is why it stays a user option. It superseded the entry of the same
  name in ``&run_config`` of ``namelist.config``; that location still works
  for one release, and either switches it on.

Section &coupling_yac
"""""""""""""""""""""

Read only for ``FESOM_COUPLING=yac``.

- **cpl_comp_name='fesom2'** component name registered with YAC; must match
  ``coupling.yaml``.
- **cpl_grid_name='fesom_grid'** grid name passed to ``yac_fdef_grid``.
- **cpl_config_file='coupling.yaml'** YAC configuration file, relative to the
  run directory. No such file ships with FESOM.

Current limitation
------------------

One binary serves any partner whose field tables its interface knows, but
``oasis28`` and ``oasis50`` remain separate builds: ECHAM needs OASIS3-MCT
2.8, whose ``prism_*`` API generation differs from 5.0's. So in practice
``oasis28`` runs ECHAM and ``oasis50`` runs OpenIFS until ECHAM is ported to
OASIS 5.0, at which point it becomes a namelist switch with no code change.
