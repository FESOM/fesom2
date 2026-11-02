#!/usr/bin/env python3
"""Write a synthetic ocean-fraction mask on the grid of a forcing file.

    make_test_ocean_mask.py FORCING_FILE OUT_FILE

For the ocean-only forcing tests (use_ocean_only_forcing): the bundled CORE2
test data has no land mask, so this makes one on the same LON x LAT grid as
FORCING_FILE. Variable sftof(LAT, LON) in percent, 0 on land, 100 over ocean.
Land is a few fixed blocks: two continents wide enough that filling them takes
several sweeps, a polar cap, and a one-cell island, so the fill is exercised at
several depths and the blocks touch model ocean nodes on any global mesh.
"""
import os
import sys

import netCDF4
import numpy as np


def main():
    src, out = sys.argv[1], sys.argv[2]
    with netCDF4.Dataset(src) as f:
        lon = np.asarray(f.variables["LON"][:], dtype=np.float64)
        lat = np.asarray(f.variables["LAT"][:], dtype=np.float64)
    lo, la = np.meshgrid(lon % 360.0, lat)           # (LAT, LON), like the forcing fields
    land = ((lo >= 60) & (lo < 140) & (la >= -20) & (la <= 55))   # wide continent
    land |= ((lo >= 250) & (lo < 300) & (la >= 0) & (la <= 60))   # second continent
    land |= la <= -70                                             # polar cap
    i, j = len(lon) // 3, len(lat) // 2
    land[j, i] = True                                             # single-cell island
    ensure_dir(out)
    with netCDF4.Dataset(out, "w") as g:
        g.createDimension("LAT", len(lat))
        g.createDimension("LON", len(lon))
        g.createVariable("LAT", "f8", ("LAT",))[:] = lat
        g.createVariable("LON", "f8", ("LON",))[:] = lon
        v = g.createVariable("sftof", "f4", ("LAT", "LON"))
        v.units = "%"
        v[:] = np.where(land, 0.0, 100.0).astype(np.float32)
    print(f"wrote {out}: {int(land.sum())} land cells of {land.size}")


def ensure_dir(path):
    d = os.path.dirname(os.path.abspath(path))
    os.makedirs(d, exist_ok=True)


if __name__ == "__main__":
    main()
