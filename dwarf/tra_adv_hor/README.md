# Horizontal tracer advection dwarf

Standalone timing of FESOM's per-tracer horizontal advection sequence (low-order flux, low-order
tendency gather, MFCT flux, flux-to-tendency gather) on real per-rank data, for development. The
model decides whether a change is kept.

1. Dump the inputs from a FESOM run: `FESOM_DWARF_DUMP=<k> FESOM_DWARF_DIR=<dir>` dumps the k-th
   call of `do_oce_adv_tra` for tracer 1 (`src/dwarf_dump.F90`), one file per rank.
2. `make` (FESOM's gfortran flags, pure MPI), then `mpirun -n <ranks of the dump> ./dwarf.x <dir> [repeats]`.

`kernels_ref.F90` holds verbatim copies of FESOM's kernels; the dwarf checks that they reproduce the
dumped results bit for bit, and compares each variant in `kernels_opt.F90` with them. Each rank
works on its own data without communication, so 1 versus 8 ranks measures how the kernels share the
node's memory bandwidth. `dwarf_ranks.txt` has per-rank work and times.
