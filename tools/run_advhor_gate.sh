#!/usr/bin/env bash
# End-to-end M1.1 horizontal-advection operator byte-gate: dump FESOM2's advection
# fields (1-rank pi, real kernels on a prescribed velocity/tracer), dump FESOM3's
# (1-rank pi, transcribed kernels on the identical prescribed field), and compare
# for max|delta|=0 on adv_flux_hor + del_ttf_advhoriz (upw1 + muscl) and every
# intermediate (helem, edge_up_dn_tri, tr_xy).
#
# Prereqs (one-time, see docs/HANDOFF.md):
#   - FESOM2 rebuilt with src/fesom_advhor_dump.F90 wired into ocean_setup
#     (build/bin/fesom.x -> build/lib64/libfesom.so).
#   - pi mesh has the hand-crafted dist_1/ (shared with the geometry gate).
#   - FESOM3 built: build_intel_dp/bin/fesom_advhordump.
set -euo pipefail
F3=/home/a/a270088/fesom3
RUN="${1:-/scratch/a/a270088/advhordump_pi}"

echo "[1/3] FESOM2 oracle advection dump (1-rank pi)"
bash "$F3/tools/run_advhordump_pi.sh" "$RUN" "$RUN/advhor_f2.bin" >/dev/null

echo "[2/3] FESOM3 advection dump (1-rank pi)"
source "$F3/env.sh" intel >/dev/null 2>&1
export FESOM3_ADVHOR_OUT="$RUN/advhor_f3.bin"
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n 1 \
    "$F3/build_intel_dp/bin/fesom_advhordump" >/dev/null 2>&1

echo "[3/3] compare"
python3 "$F3/tools/advhor_diff.py" "$RUN/advhor_f2.bin" "$RUN/advhor_f3.bin"
