.. _profiler:

Built-in profiler
*****************

FESOM2 contains its own profiler that reports where the model spends its time
and how evenly that time is distributed over the MPI ranks. It is compiled out
by default and costs nothing when disabled, because both the timing calls and
the diagnostic calls sit behind preprocessor guards.

Enabling
========

::

    cmake -DFESOM_PROFILING=ON ..

- **FESOM_PROFILING=OFF** compile the profiling and scalar-diagnostics call sites. Without it the macros expand to nothing.

Reading the output
==================

At the end of a run the profiler writes ``fesom.stats`` into ``ResultPath``,
next to the model output. The file contains one row per instrumented section
with the mean, minimum and maximum wall time over the MPI ranks, the number
of calls, the share of the total runtime and the share of the parent section,
and two separate measures of how badly the work is distributed:

- **RngImb(%)**, the range-based imbalance, computed as the spread between the slowest and the fastest rank divided by the mean. It reacts strongly to a single slow rank.
- **StdImb(%)**, the distribution-based imbalance, computed as the standard deviation over the ranks divided by the mean. It stays small when one rank is late and grows when many ranks are.

Reading the two together tells you what kind of imbalance you have, and the
file repeats this legend in its own header: a high range imbalance with a low
standard imbalance means one outlier rank while the others are balanced, both
being high means a general imbalance across many ranks, and a low range
imbalance with a high standard imbalance means several performance clusters
without a single outlier.

Alongside the timing table the same file carries the scalar diagnostics
registered through ``fesom_rtdiagnostics``, which record per-call quantities
such as solver iterations, residuals and CFL numbers. These are accumulated
in double precision regardless of the working precision of the build, so the
numbers from a single-precision and a double-precision run can be compared
directly.
