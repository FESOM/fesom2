#!/usr/bin/env bash
# CONSERVATION + NO-LEAKAGE GATE for the bottom-at-vertices change (docs/plans/completed/
# 20260910-fesom3-bottom-at-vertices.md, Task 2). This is the numerical regression net
# that replaces the retired FESOM2 byte-gates: it is FESOM_F-internal, needs no oracle,
# and must stay green through every task of that plan.
#
# What it checks, via src/drivers/fesom_conserve.F90:
#   zstar + Redi, np 1 -> the Redi path. diff_ver_part_redi_expl is the ONLY caller of
#       the tr_xynodes node-average whose denominator the depth-independent area changed,
#       and it never runs unless Redi is on. Note that conservation alone does NOT validate
#       that denominator -- the Redi tendency is a telescoping flux divergence, so it
#       conserves whatever tr_xynodes contains. The denominator is pinned by test_bottom's
#       'Redi: node-averaged gradient is exact on the WET area' check instead.
#   zstar, np 1 and 2 -> HARD conservation gate. Total heat and salt content must not
#       drift beyond TOL. This is the sharp test of the bottom change: if the flux areas
#       and the cell volumes ever disagree -- the central risk when area() becomes
#       depth-independent and nlevels() flips to vertex-defined -- it shows up here as
#       drift long before it shows up anywhere else.
#   linfs, np 1 -> INVARIANTS ONLY (T10 velocity leakage, helem == mean(hnode),
#       zbar_e_bot, finiteness). Drift is reported but NOT gated, because the linear free
#       surface is not tracer-conserving by construction: hnode is frozen, so the surface
#       vertical advective flux -w*T*area at nzmin is a real source/sink with no thickness
#       change to balance it. Measured baseline on pi over 20 steps: heat -2.2e-04,
#       salt -1.2e-06, non-monotone (a free-surface adjustment transient).
#
#   zstar + wsplit, np 1 and 8 (+ momadv-vinv-wsplit and GM+wsplit, np 1) -> the smooth
#       w = w_e + w_i split (use_wsplit, oce_wsplit) with ALL its consumers live: the
#       momentum TDMA (w_i, flux form at momadv_opt=2 / advective form at 1), momentum
#       advection (w_e), the tracer-diffusion TDMA (do_wimpl, w_i) and the FCT
#       adv_tra_vert_impl (w_i, then the LO flux on the full w); GM+wsplit adds the bolus
#       fer_w, which is added to w AND w_e at owned+halo right where the split's halo
#       contract lives and bypasses the cap. The 20-step cold start never reaches
#       CFL_z ~ 1 (measured last-step max 2.76e-2 on pi, split off), so the cap is set
#       to FESOM3_WSPLIT_MAXCFL=0.005 (onset 0.0025 = the driver's 0.5*maxcfl default;
#       the triple is WSPLIT_CFG below, the same as test/CMakeLists.txt's _wsplit_cfg):
#       measured 8632 owned faces with w_i /= 0, 500 with f = |w_i|/|w| >= 0.5, max f
#       0.82 (np 1/2 identical to +-1 face). FESOM3_WSPLIT_EXPECT_SPLIT=1000 makes the
#       driver error-stop unless >= 1000 faces split and one has f >= 0.5, so these
#       configs cannot pass vacuously. The driver prints the split statistics line for
#       EVERY config (the line before the drift), which is why every tail below is one
#       line longer than the drift block alone.
#
#   linfs + sea ice, np 2 and 8 -> SHARED-ELEMENT gate (docs/plans/completed/20261012-shared-element-sync.md).
#       A border element is owned by every rank that owns one of its vertices; update_vel syncs
#       uv to the designated owner's copy. FESOM3_DIVSUM_TOL=1e-12 makes the driver stop unless
#       the divergence-sum ratio |sum ssh_rhs_old| / sum|ssh_rhs_old| stays below 1e-12 in every
#       step (a linfs mass leak raises it) and the copy spreads of uv, uv_rhsAB and sigma are 0.
#       Vector-invariant momentum with rvo_upwind 1 (the CORE2 leak's trigger) at np 2 and 8,
#       the scalar advection (momadv 2, the default) at np 2. The same vinv case with
#       FESOM3_SHARED_SYNC=0 is report-only: it shows whether the copies drift without the sync.
#       Measured 2026-10-12 (pi, 20 steps, cold start): max ratio 1.6e-16 (vinv np 2), 3.1e-16
#       (vinv np 8), 1.8e-16 (scalar np 2), every spread 0. Sync off: uv spread 6.9e-18 -- the
#       roundoff drift -- and uv_rhsAB spread 5.8e-11, i.e. the rvo_upwind blend has amplified it
#       by ~1e7 within 20 steps (sigma 0; the ratio stays 1.6e-16 this early at these small
#       velocities). The other configs report ratios <= 3.8e-16 and zero spreads.
#   Every driver run prints one shared-element line after the drift block, so every tail
#   below is one line longer than before.
#
# Measured zstar baseline on pi, 20 steps, before any bottom change:
#   np=1  heat -4.87e-15  salt -1.66e-14
#   np=2  heat  0.00e+00  salt -7.80e-15
# TOL=1e-12 leaves roughly two orders of headroom over that.
#
# Usage:  bash tools/run_conserve_pi.sh [NSTEPS]
set -euo pipefail

F3="${F3:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BIN="${BIN:-$F3/build_intel_dp/bin/fesom_conserve}"
PIMESH="${PIMESH:-/home/a/a270088/port2/fesom2/tests/data/MESHES/pi}"
ICFILE="${ICFILE:-/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc}"
NSTEPS="${1:-20}"
TOL="${TOL:-1e-12}"
MPIFLAGS="${MPIFLAGS:---mca pml ob1 --mca btl self,vader --mca btl_vader_single_copy_mechanism none --oversubscribe}"  # LESSONS L61
# the tuned wsplit triple (header); re-tune the cap here AND in test/CMakeLists.txt
WSPLIT_CFG="${WSPLIT_CFG:-FESOM3_WSPLIT=1 FESOM3_WSPLIT_MAXCFL=0.005 FESOM3_WSPLIT_EXPECT_SPLIT=1000}"

[ -x "$BIN" ] || { echo "run_conserve_pi: missing $BIN (build first)"; exit 1; }
[ -d "$PIMESH" ] || { echo "run_conserve_pi: missing pi mesh $PIMESH"; exit 1; }

export FESOM3_MESH_DIR="$PIMESH" FESOM3_IC_FILE="$ICFILE" FESOM3_NSTEPS="$NSTEPS"

fail=0

# np=8 matters beyond conservation: at np=1 nNodL == mesh%nod2D, so a routine that
# over-runs a LOCAL-sized array using a dummy declared with the GLOBAL nod2D is exactly
# in-bounds and cannot fail. Only a rank count where nNodL << nod2D exposes it.
for np in 1 2 8; do
    echo "=== zstar, np=$np, $NSTEPS steps (conservation gate, tol=$TOL) ==="
    if FESOM3_WHICH_ALE=zstar FESOM3_CONSERVE_TOL="$TOL" \
         mpirun $MPIFLAGS -n "$np" "$BIN" 2>&1 | tail -8; then
        :
    else
        echo "run_conserve_pi: FAILED (zstar np=$np)"; fail=1
    fi
done

# Physics combinations. The area/areasvol consumers are spread across the vertical
# diffusion TDMA, the Redi isoneutral flux and the KPP non-local / shortwave terms, and
# each is only reachable with its own scheme switched on -- a gate that runs only the
# default configuration proves nothing about the others.
# NOTE read from fd 3: mpirun inherits stdin and CONSUMES the heredoc, so a plain
# `while read ... done <<CFG` silently runs only the FIRST config and drops the rest.
# The heredoc is unquoted so $WSPLIT_CFG expands (no other line contains $ or `).
while read -r -u 3 label vars; do
    [ -z "$label" ] && continue
    echo "=== zstar + $label, np=1, $NSTEPS steps (conservation gate, tol=$TOL) ==="
    if env $vars FESOM3_WHICH_ALE=zstar FESOM3_CONSERVE_TOL="$TOL" \
         mpirun $MPIFLAGS -n 1 "$BIN" 2>&1 | tail -9; then
        :
    else
        echo "run_conserve_pi: FAILED (zstar+$label np=1)"; fail=1
    fi
done 3<<CFG
Redi FESOM3_REDI=1
KPP FESOM3_MIX_KPP=1
TKE FESOM3_MIX_TKE=1
GM+Redi+TKE FESOM3_FER_GM=1 FESOM3_REDI=1 FESOM3_MIX_TKE=1
GM+KPP FESOM3_FER_GM=1 FESOM3_MIX_KPP=1
splines FESOM3_SHEAR_SPLINES=1
splines+KPP FESOM3_SHEAR_SPLINES=1 FESOM3_MIX_KPP=1
splines+TKE FESOM3_SHEAR_SPLINES=1 FESOM3_MIX_TKE=1
N2splines FESOM3_N2_SPLINES=1
momadv-vinv FESOM3_MOMADV_OPT=1
momadv-vinv-upw FESOM3_MOMADV_OPT=1 FESOM3_RVO_UPWIND=0.7
wsplit $WSPLIT_CFG
momadv-vinv-wsplit FESOM3_MOMADV_OPT=1 $WSPLIT_CFG
GM+wsplit FESOM3_FER_GM=1 $WSPLIT_CFG
bothsplines+GM+Redi+TKE FESOM3_SHEAR_SPLINES=1 FESOM3_N2_SPLINES=1 FESOM3_FER_GM=1 FESOM3_REDI=1 FESOM3_MIX_TKE=1
CFG

# Vector-invariant momentum advection at HIGH rank. The new code adds two nodal scatters
# and two halo exchanges, so np=8 (where nNodL << nod2D) is the configuration that matters
# most -- the same reason the zstar loop above runs it.
echo "=== zstar + momadv-vinv, np=8, $NSTEPS steps (conservation gate, tol=$TOL) ==="
if FESOM3_MOMADV_OPT=1 FESOM3_WHICH_ALE=zstar FESOM3_CONSERVE_TOL="$TOL" \
     mpirun $MPIFLAGS -n 8 "$BIN" < /dev/null 2>&1 | tail -8; then
    :
else
    echo "run_conserve_pi: FAILED (zstar+momadv-vinv np=8)"; fail=1
fi

# upwind blend at np=8: its exchange_elem of omega_e is the one new MR communication.
echo "=== zstar + momadv-vinv-upw, np=8, $NSTEPS steps (conservation gate, tol=$TOL) ==="
if FESOM3_MOMADV_OPT=1 FESOM3_RVO_UPWIND=0.7 FESOM3_WHICH_ALE=zstar FESOM3_CONSERVE_TOL="$TOL" \
     mpirun $MPIFLAGS -n 8 "$BIN" < /dev/null 2>&1 | tail -8; then
    :
else
    echo "run_conserve_pi: FAILED (zstar+momadv-vinv-upw np=8)"; fail=1
fi

# w split at np=8: compute_Wvel_split produces w_e/w_i at owned+HALO with no trailing
# exchange (FESOM2's layout) and the NEXT step's momentum TDMA reads them at the halo
# vertices of owned elements, so the rank count where nNodL << nod2D is the one that runs
# that consumer at scale (the halo values themselves are pinned by test_wsplit S5).
echo "=== zstar + wsplit, np=8, $NSTEPS steps (conservation gate, tol=$TOL) ==="
if env $WSPLIT_CFG FESOM3_WHICH_ALE=zstar FESOM3_CONSERVE_TOL="$TOL" \
     mpirun $MPIFLAGS -n 8 "$BIN" < /dev/null 2>&1 | tail -8; then
    :
else
    echo "run_conserve_pi: FAILED (zstar+wsplit np=8)"; fail=1
fi

echo "=== linfs, np=1, $NSTEPS steps (invariants only; drift reported, not gated) ==="
if FESOM3_WHICH_ALE=linfs mpirun $MPIFLAGS -n 1 "$BIN" 2>&1 | tail -7; then
    :
else
    echo "run_conserve_pi: FAILED (linfs invariants np=1)"; fail=1
fi

# SHARED: linfs + sea ice, the shared-element gate (header). Own blocks, not the heredoc loop:
# that loop runs np 1 only. The sync-off case is report-only (no FESOM3_DIVSUM_TOL).
SHARED_ENV="FESOM3_WHICH_ALE=linfs FESOM3_CONSERVE_ICE=1"
while read -r -u 3 label np vars; do
    [ -z "$label" ] && continue
    echo "=== linfs + ice + $label, np=$np, $NSTEPS steps (shared-element gate, divsum tol=$TOL) ==="
    if env $SHARED_ENV $vars FESOM3_DIVSUM_TOL="$TOL" \
         mpirun $MPIFLAGS -n "$np" "$BIN" < /dev/null 2>&1 | tail -8; then
        :
    else
        echo "run_conserve_pi: FAILED (linfs+ice+$label np=$np)"; fail=1
    fi
done 3<<CFG
momadv-vinv-upw1 2 FESOM3_MOMADV_OPT=1 FESOM3_RVO_UPWIND=1.0
momadv-vinv-upw1 8 FESOM3_MOMADV_OPT=1 FESOM3_RVO_UPWIND=1.0
momadv-scalar 2 FESOM3_MOMADV_OPT=2
CFG
echo "=== linfs + ice + momadv-vinv-upw1, np=2, sync OFF, $NSTEPS steps (report only) ==="
if env $SHARED_ENV FESOM3_MOMADV_OPT=1 FESOM3_RVO_UPWIND=1.0 FESOM3_SHARED_SYNC=0 \
     mpirun $MPIFLAGS -n 2 "$BIN" < /dev/null 2>&1 | tail -8; then
    :
else
    echo "run_conserve_pi: FAILED (linfs+ice+momadv-vinv-upw1 np=2 sync off)"; fail=1
fi

if [ "$fail" -ne 0 ]; then
    echo "run_conserve_pi: GATE FAILED"
    exit 1
fi
echo "run_conserve_pi: GATE OK"
