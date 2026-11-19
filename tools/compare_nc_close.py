#!/usr/bin/env python3
"""Compare two FESOM result directories within a tolerance.

    compare_nc_close.py REF_DIR TEST_DIR [--rtol R] [--atol A] [--vars v1,v2] [--exclude v3,v4]

Every *.nc file below REF_DIR must exist below TEST_DIR. For each floating
point variable (restricted to --vars, without --exclude),

    max|test - ref| <= rtol * max|ref| + atol

must hold. With the default rtol = atol = 0 this is an exact comparison of
the values. Exit status 0 only if every compared variable passes.
"""
import argparse
import os
import sys

import netCDF4
import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_nc_exact import nc_files, raw  # noqa: E402


def main():
    p = argparse.ArgumentParser()
    p.add_argument("ref")
    p.add_argument("test")
    p.add_argument("--rtol", type=float, default=0.0)
    p.add_argument("--atol", type=float, default=0.0)
    p.add_argument("--vars", default="", help="comma separated variable names")
    p.add_argument("--exclude", default="", help="comma separated variable names")
    a = p.parse_args()
    only = set(v for v in a.vars.split(",") if v)
    skip = set(v for v in a.exclude.split(",") if v)

    fr, ft = nc_files(a.ref), nc_files(a.test)
    bad = 0
    for f in sorted(fr - ft):
        print(f"MISSING in test: {f}")
        bad += 1
    seen = set()
    nvars = 0
    for f in sorted(fr & ft):
        with netCDF4.Dataset(os.path.join(a.ref, f)) as da, netCDF4.Dataset(os.path.join(a.test, f)) as db:
            for name in sorted(da.variables):
                if (only and name not in only) or name in skip:
                    continue
                xa = raw(da.variables[name])
                if xa.dtype.kind != "f":
                    continue
                if name not in db.variables:
                    print(f"MISSING {f}: {name}")
                    bad += 1
                    continue
                xb = raw(db.variables[name])
                seen.add(name)
                nvars += 1
                if xa.shape != xb.shape:
                    print(f"DIFF {f}: {name}: shape {xa.shape} vs {xb.shape}")
                    bad += 1
                    continue
                xa, xb = xa.astype(np.float64), xb.astype(np.float64)
                d = float(np.nanmax(np.abs(xa - xb))) if xa.size else 0.0
                scale = float(np.nanmax(np.abs(xa))) if xa.size else 0.0
                lim = a.rtol * scale + a.atol
                ok = d <= lim and np.array_equal(np.isnan(xa), np.isnan(xb))
                rel = d / scale if scale > 0 else d
                print(f"{'ok  ' if ok else 'FAIL'} {f}: {name}: max|diff| {d:.3e} rel {rel:.3e} (limit {lim:.3e})")
                bad += not ok
    for name in sorted(only - seen):
        print(f"MISSING variable in all files: {name}")
        bad += 1
    if nvars == 0:
        print("no variables compared")
        bad += 1
    print(f"{'PASS' if bad == 0 else 'FAIL'}: {nvars} variables compared, {a.test} vs {a.ref}")
    sys.exit(0 if bad == 0 else 1)


if __name__ == "__main__":
    main()
