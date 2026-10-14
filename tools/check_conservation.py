#!/usr/bin/env python3
"""Check the conservation budgets a FESOM log reports (CONS lines, conservation_freq > 0).

    check_conservation.py <log> [--rtol-volume R] [--rtol-heat R] [--rtol-salt R] [--json FILE]

A CONS line holds, per quantity (volume, heat, salt): the integral, its change since the
start of the run, the flux the model applied at the surface since the start, their
difference (the residual) and the residual relative to the change. This check judges the
residual against the INTEGRAL (|residual| / |integral| at the last reported step), because
the change itself can be a sum of cancelling terms; the defaults are 1e-12 for volume and
salt and 1e-10 for heat. Prints one line per quantity and the largest relative residual
of the run; --json writes the whole series for plots.
"""
import json
import sys

NAMES = ("volume", "heat", "salt")
DEF = {"volume": 1e-12, "heat": 1e-10, "salt": 1e-12}


def main():
    args = sys.argv[1:]
    log = args[0]
    rtol = dict(DEF)
    out = None
    for i, a in enumerate(args):
        for q in NAMES:
            if a == "--rtol-" + q:
                rtol[q] = float(args[i + 1])
        if a == "--json":
            out = args[i + 1]
    ser = {q: [] for q in NAMES}
    for line in open(log, errors="replace"):
        w = line.split()
        if len(w) == 8 and w[0] == "CONS" and w[2] in ser:
            step = int(w[1])
            integral, change, flux, res, rel = [float(x) for x in w[3:8]]
            ser[w[2]].append({"step": step, "integral": integral, "change": change, "flux": flux,
                              "residual": res, "rel_change": rel,
                              "rel_integral": abs(res) / max(abs(integral), 1e-300)})
    if not any(ser.values()):
        print("check_conservation: no CONS lines in %s (conservation_freq = 0?)" % log)
        return 2
    if out:
        json.dump(ser, open(out, "w"))
    bad = False
    for q in NAMES:
        if not ser[q]:
            continue
        last = ser[q][-1]
        worst = max(s["rel_integral"] for s in ser[q])
        ok = last["rel_integral"] <= rtol[q]
        bad = bad or not ok
        print("check_conservation: %-6s step %6d  integral %.6e  change %.3e  flux %.3e  residual %.3e  |res|/|integral| %.2e (max %.2e, tol %.0e) %s" % (
            q, last["step"], last["integral"], last["change"], last["flux"], last["residual"], last["rel_integral"], worst, rtol[q], "ok" if ok else "FAILED"))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
