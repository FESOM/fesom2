#!/usr/bin/env python3
"""Compare two FESOM result directories bit for bit.

    compare_nc_exact.py REF_DIR TEST_DIR

Every *.nc file below REF_DIR (outputs and restarts, also inside restart
directories) must exist below TEST_DIR with the same variables, shapes,
dtypes and identical bytes. Files only present in one tree are reported too.
Exit status 0 only if everything matches.
"""
import os
import sys

import netCDF4
import numpy as np


def nc_files(root):
    out = set()
    for d, _, files in os.walk(root):
        for f in files:
            if f.endswith(".nc"):
                out.add(os.path.relpath(os.path.join(d, f), root))
    return out


def raw(var):
    var.set_auto_maskandscale(False)
    return np.ascontiguousarray(var[...])


def compare_file(a, b):
    errs = []
    with netCDF4.Dataset(a) as fa, netCDF4.Dataset(b) as fb:
        va, vb = set(fa.variables), set(fb.variables)
        if va != vb:
            errs.append(f"variables differ: only ref {sorted(va - vb)}, only test {sorted(vb - va)}")
        for name in sorted(va & vb):
            xa, xb = raw(fa.variables[name]), raw(fb.variables[name])
            if xa.shape != xb.shape or xa.dtype != xb.dtype:
                errs.append(f"{name}: shape/dtype {xa.shape} {xa.dtype} vs {xb.shape} {xb.dtype}")
            elif xa.tobytes() != xb.tobytes():
                if xa.dtype.kind == "f":
                    d = np.abs(xa.astype(np.float64) - xb.astype(np.float64))
                    n = int(np.count_nonzero(xa.view(f"u{xa.itemsize}") != xb.view(f"u{xb.itemsize}")))
                    errs.append(f"{name}: {n} values differ, max |diff| {np.nanmax(d):.3e}")
                else:
                    errs.append(f"{name}: values differ")
    return len(va & vb), errs


def main():
    ref, test = sys.argv[1], sys.argv[2]
    fr, ft = nc_files(ref), nc_files(test)
    bad = 0
    for f in sorted(fr ^ ft):
        print(f"MISSING in {'test' if f in fr else 'ref'}: {f}")
        bad += 1
    nvars = 0
    for f in sorted(fr & ft):
        n, errs = compare_file(os.path.join(ref, f), os.path.join(test, f))
        nvars += n
        for e in errs:
            print(f"DIFF {f}: {e}")
        bad += bool(errs)
    status = "IDENTICAL" if bad == 0 else f"DIFFERENT ({bad} files)"
    print(f"{status}: {len(fr & ft)} files, {nvars} variables compared, {test} vs {ref}")
    sys.exit(0 if bad == 0 else 1)


if __name__ == "__main__":
    main()
