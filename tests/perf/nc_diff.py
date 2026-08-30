#!/usr/bin/env python3
"""Numeric diff of two FESOM NetCDF output files.

Usage: nc_diff.py <baseline>.nc <new>.nc
Prints per-variable max abs/rel difference and a WORST_REL / NAN summary.
Acceptance bound for the precision-hygiene refactor: WORST_REL <= 1.2e-7
(~1 float32 ULP, the floor set by float32 NetCDF output) and NAN == 0.
"""
import sys, numpy as np
from netCDF4 import Dataset

SKIP = ('time', 'time_bnds', 'nz', 'nz1', 'nod2', 'elem')


def load(p):
    """Return {varname: raw float array}. Note we deliberately do NOT md5 the
    file itself: FESOM stamps the git SHA into a global attribute, so two
    numerically identical runs from different commits have different md5s."""
    d = Dataset(p)
    d.set_auto_mask(False)
    out = {}
    for k, v in d.variables.items():
        if k in SKIP:
            continue
        a = np.array(v[:])
        if a.dtype.kind == 'f':
            out[k] = a
    d.close()
    return out


def bitwise_equal(x, y):
    if x.shape != y.shape or x.dtype != y.dtype:
        return False
    return np.array_equal(x.view(np.uint8), y.view(np.uint8))


a, b = load(sys.argv[1]), load(sys.argv[2])
missing = sorted(set(a) ^ set(b))
if missing:
    print(f"  VARIABLE SET DIFFERS: {missing}")
bits_equal = not missing and all(bitwise_equal(a[k], b[k]) for k in a)
worst = 0.0
nan = 0
# A shape mismatch is a HARD failure, not a skip. It used to `continue` without
# recording anything, so a run whose record count changed (say 4 time records
# against 1) left worst=0 and reported VALUES-IDENTICAL with exit 0 -- the gate
# passed a file it had never actually compared.
shape_diff = []
for k in sorted(set(a) & set(b)):
    if a[k].shape != b[k].shape:
        print(f"  SHAPE {k}: {a[k].shape} vs {b[k].shape}")
        shape_diff.append(k)
        continue
    d = np.abs(a[k].astype(np.float64) - b[k].astype(np.float64))
    rel = (d / np.maximum(np.abs(b[k].astype(np.float64)), 1e-30)).max()
    nanc = int(np.isnan(a[k]).sum())
    nan += nanc
    if d.max() > 0 or nanc:
        print(f"  {k}: max|d|={d.max():.3e} maxrel={rel:.3e} nan={nanc}")
    worst = max(worst, rel)
if shape_diff or missing:
    status = 'SHAPE/VARIABLE MISMATCH'
elif bits_equal and nan == 0:
    status = 'BITWISE-IDENTICAL'
elif worst == 0 and nan == 0:
    status = 'VALUES-IDENTICAL'
else:
    status = 'DIFF'
print(f"WORST_REL={worst:.3e}  NAN={nan}  ({status})")
sys.exit(0 if (worst <= 1.2e-7 and nan == 0 and not missing and not shape_diff) else 1)
