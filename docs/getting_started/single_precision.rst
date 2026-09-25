.. _single_precision:

Running FESOM2 in single precision
**********************************

All real arithmetic in FESOM2 is written against the working-precision kind
``WP``, which is double precision (real64) by default. A build option narrows
``WP`` to single precision (real32) throughout the model. This roughly halves
the memory footprint and memory traffic of the model state, which speeds up
runs that are limited by memory bandwidth, as large-mesh FESOM2 runs usually
are.

Building
========

Precision is a property of the build, not of the run. Configure with::

    cmake -DUSE_SINGLE_PRECISION=ON

together with your usual options. At startup the model prints its working
precision to the log::

    FESOM working precision: WP=4 bytes (SINGLE PRECISION MODE)

A double-precision build prints ``WP=8 bytes (DOUBLE PRECISION MODE)``. Check
this banner if you are unsure which executable a job ran.

What stays in double precision
==============================

Narrowing every number to float32 would let rounding errors grow without bound
in long sums and would destabilise the iterative SSH solver. A second kind,
``WP_full``, therefore remains real64 in every build, and a small set of
accumulators and control quantities uses it deliberately. The bandwidth-heavy
bulk of the model follows ``WP``:

======================================================  ==============================================
Component                                               Precision in a single-precision build
======================================================  ==============================================
Ocean and ice state, tracers, fluxes, mesh arrays       single (``WP``)
SSH solver matrix, vectors and preconditioner (SpMV)    single (``WP``)
Global integrals and total ocean areas                  double accumulation (``WP_full``)
CG scalar products, search direction, stopping test     double (``WP_full``)
SSH stiffness matrix accumulation                       double shadow, rounded to ``WP`` once per step
Time axis and output running means                      double (real64)
OASIS exchange buffers                                  double (``WP_full``)
CVMix mixing schemes                                    double (fixed-kind library, wrappers convert)
======================================================  ==============================================

The double-precision shadow of the SSH stiffness matrix exists because the
matrix is updated incrementally every time step; accumulating in real64 and
rounding once per step keeps the rounding from compounding over the run, while
the solver itself keeps the single-precision memory bandwidth. Likewise the
conjugate-gradient scalars steer the iteration, so they are computed from
real64 reductions even though the vectors they act on stay single.

Salinity as an anomaly
======================

Around a typical open-ocean salinity of 35 psu, float32 has a spacing of about
2e-6 psu. Storing the anomaly relative to 35 instead moves the state close to
zero, where the spacing is 30 to 250 times finer:

- **use_salt_anomaly=.false.** (section ``&oce_dyn`` of ``namelist.oce``)
  store the salinity state as ``S-35`` internally. Output and restart fields
  remain absolute salinity (the offset is added back once, after averaging).
  With the default ``.false.`` the model is bit-identical to earlier versions.
  Some configurations, for example ice cavities or icebergs and several
  salinity diagnostics, do not support the anomaly form yet; the model then
  refuses to start and prints the exact list.

Coupled configurations
======================

Coupling paths that still exchange ``WP`` buffers as hardcoded MPI doubles
would corrupt memory in a single-precision build, so cmake refuses such
combinations at configure time with a fatal error. OpenIFS coupling via OASIS
(``OIFS_COUPLED``, as in AWI-ESM3) has been converted and is allowed; OASIS
coupling without OpenIFS, YAC and the IFS interface are refused. Setting
``FESOM_ALLOW_SINGLE_PRECISION_COUPLED=ON`` overrides the guard, but is meant
for developing those conversions, not for production runs.

What to expect
==============

- Results are not bit-comparable with a double-precision run. Differences
  start at rounding level and grow chaotically, as for any small perturbation;
  they are not a sign of a broken build.
- Restart files written through netCDF keep double-precision variables on
  disk in both builds, but the state of a single-precision run carries only
  float32 information. The raw binary dumps of the mesh derived types store
  arrays at working precision, so keep one precision throughout an experiment
  and its restarts.
- Output precision is unaffected: each variable's precision in the ``io_list``
  of ``namelist.io`` is chosen independently of ``WP``, and running means are
  always accumulated in real64 before being written.
- The standard test suite runs in both precisions: continuous integration
  builds the same full configuration once per precision and runs the same
  tests against each executable.

When to use it
==============

Use single precision for throughput-limited standalone ocean and ocean-ice
production runs, in particular on large meshes or for ensembles, where the
memory-bandwidth saving translates directly into speed. Stay with the default
double precision for reference solutions, for debugging suspected numerical
problems, and for the coupled configurations that the build system still
refuses.
