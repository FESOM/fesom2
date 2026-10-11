#!/usr/bin/env python3
"""Stretch the bundled CORE2 test forcing into a longer series for benchmarking.

The files in tests/data/FORCING/CORE2 hold only 5 six-hourly records, so a pi
run sees four or five record changes per field and the forcing refresh path
can hardly be timed. This writes copies on the same grid with NREC records at
STEP_HOURS spacing, cycling through the original records. The values are real
CORE2 fields, only their time stamps are new, so the run stays physical.

Files without a year in their name (runoff, SSS climatology, ...) are linked.

    tools/make_synthetic_forcing.py tests/data/FORCING/CORE2 OUTDIR --nrec 240 --step-hours 1
"""
import argparse
import os
import re
import sys

import netCDF4
import numpy as np


def stretch(src, dst, nrec, step_hours):
    with netCDF4.Dataset(src) as fin, netCDF4.Dataset(dst, "w", format=fin.data_model) as fout:
        tdim = next(n for n, d in fin.dimensions.items() if d.isunlimited())
        nsrc = len(fin.dimensions[tdim])
        for name, dim in fin.dimensions.items():
            fout.createDimension(name, None if dim.isunlimited() else len(dim))
        fout.setncatts({k: fin.getncattr(k) for k in fin.ncattrs()})
        cycle = np.arange(nrec) % nsrc
        for name, var in fin.variables.items():
            fill = var.getncattr("_FillValue") if "_FillValue" in var.ncattrs() else None
            out = fout.createVariable(name, var.dtype, var.dimensions, fill_value=fill)
            out.setncatts({k: var.getncattr(k) for k in var.ncattrs() if k != "_FillValue"})
            if not var.dimensions or var.dimensions[0] != tdim:
                out[...] = var[...]
            elif name == tdim:
                out[:] = (np.arange(nrec) * step_hours).astype(var.dtype)
            elif name == tdim + "_bnds":
                # bounds are in days in the CORE2 files
                t0 = np.arange(nrec) * step_hours / 24.0
                out[:] = np.stack([t0, t0 + step_hours / 24.0], axis=1)
            else:
                for r in range(nrec):
                    out[r, ...] = var[cycle[r], ...]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("srcdir")
    ap.add_argument("outdir")
    ap.add_argument("--nrec", type=int, default=240)
    ap.add_argument("--step-hours", type=float, default=1.0)
    args = ap.parse_args()

    if not os.path.isdir(args.srcdir):
        sys.exit(f"no such directory: {args.srcdir}")
    os.makedirs(args.outdir, exist_ok=True)
    for fname in sorted(os.listdir(args.srcdir)):
        src = os.path.abspath(os.path.join(args.srcdir, fname))
        dst = os.path.join(args.outdir, fname)
        if os.path.lexists(dst):
            os.remove(dst)
        if re.search(r"\.\d{4}\.nc$", fname):
            stretch(src, dst, args.nrec, args.step_hours)
            print(f"stretched {fname}: {args.nrec} records every {args.step_hours} h")
        else:
            os.symlink(src, dst)


if __name__ == "__main__":
    main()
