#!/usr/bin/env python3
"""Verify that every !$OMP ORDERED region binds to a loop carrying ORDERED.

Guards FESOM/fesom2#444 and FESOM/fesom2#528.

The __openmp_reproducible build replaces the OpenMP lock-based scatter and
reduction loops with !$OMP ORDERED regions so that accumulation happens in
loop-iteration order.  OpenMP requires the enclosing worksharing loop to carry
the ORDERED clause; none of them did, so every affected file failed to compile:

    Error: 'ordered' region must be closely nested inside a loop region with
           an 'ordered' clause

This check drives the real Fortran preprocessor, so it sees exactly the
directive stream the compiler sees and needs no #if reasoning of its own.
Original file/line numbers are recovered from the cpp line markers.

Two modes, both meaningful invariants:

  reproducible  every active !$OMP ORDERED region is closely nested inside a
                worksharing loop whose directive carries the ORDERED clause.
                Violations are exactly the compile errors of #444/#528.

  default       no ORDERED region and no ORDERED clause is active at all.
                Guards against "fixing" #444 by adding the clause
                unconditionally, which would serialize the default build.
"""

import argparse
import os
import re
import subprocess
import sys

LINEMARKER = re.compile(r'^# (\d+) "([^"]+)"')
OMP_PREFIX = re.compile(r'^[ \t]*!\$OMP[ \t&]*', re.IGNORECASE)
CONTINUED = re.compile(r'&[ \t]*$')

END_ORDERED = re.compile(r'^END[ \t]+ORDERED\b')
# Fortran block statements, needed to tell whether an ordered region sits inside
# an inner DO loop of its worksharing loop (see scan()).
FORTRAN_DO  = re.compile(r'^[ \t]*(\w+[ \t]*:[ \t]*)?DO\b(?![ \t]*$)', re.IGNORECASE)
FORTRAN_END_DO = re.compile(r'^[ \t]*(END[ \t]*DO|ENDDO)\b', re.IGNORECASE)
ORDERED_REGION = re.compile(r'^ORDERED([ \t(]|$)')
END_LOOP = re.compile(r'^END[ \t]+(PARALLEL[ \t]+)?(DO|PARALLEL)([ \t]|$)')
BEGIN_LOOP = re.compile(r'^(PARALLEL[ \t]+)?DO([ \t]|$)')
ORDERED_CLAUSE = re.compile(r'[ \t]ORDERED([ \t(]|$)')


def preprocess(compiler, flags, defines, include_dir, path):
    cmd = [compiler] + flags + defines + ["-I" + include_dir, path]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            "failed to preprocess %s:\n  %s\n%s" % (path, " ".join(cmd), proc.stderr)
        )
    return proc.stdout


def directives(text, path):
    """Yield (file, line, uppercased directive text) for each OpenMP directive.

    Fortran continuation lines are joined so that clauses are never split.
    """
    cur_file, cur_line = path, 0
    pending, pending_at = None, None

    for raw in text.split("\n"):
        marker = LINEMARKER.match(raw)
        if marker:
            cur_line = int(marker.group(1)) - 1
            cur_file = marker.group(2)
            continue
        cur_line += 1

        if not OMP_PREFIX.match(raw):
            stripped = raw.strip()
            if stripped and not stripped.startswith('!'):
                yield cur_file, cur_line, None, stripped
            continue
        body = OMP_PREFIX.sub("", raw.upper()).strip()

        if pending is None:
            pending, pending_at = body, (cur_file, cur_line)
        else:
            pending += " " + body
        if CONTINUED.search(pending):
            continue

        yield pending_at[0], pending_at[1], pending.replace("&", " ").strip(), None
        pending, pending_at = None, None


def scan(path, text, mode, root, violations):
    """Walk one preprocessed translation unit, appending any violations."""
    loop = None  # (file, line, directive text) of the innermost open loop
    n_ordered = 0
    do_depth = 0  # Fortran DO nesting inside the current worksharing loop

    for dfile, dline, directive, code in directives(text, path):
        if directive is None:
            # plain Fortran statement - only track DO nesting inside a loop
            if loop is not None:
                if FORTRAN_END_DO.match(code):
                    do_depth = max(0, do_depth - 1)
                elif FORTRAN_DO.match(code):
                    do_depth += 1
            continue
        def rel(p):
            return os.path.relpath(p, root) if os.path.isabs(p) else p

        if END_ORDERED.match(directive):
            continue

        if ORDERED_REGION.match(directive):
            n_ordered += 1
            if mode == "default":
                violations.append(
                    "%s:%d: !$OMP ORDERED region is active in the DEFAULT build; "
                    "ORDERED belongs to the __openmp_reproducible path only"
                    % (rel(dfile), dline)
                )
            elif loop is None:
                violations.append(
                    "%s:%d: !$OMP ORDERED is not closely nested inside a worksharing loop"
                    % (rel(dfile), dline)
                )
            elif do_depth > 1:
                # OpenMP: a thread must not execute more than one ordered region
                # per iteration of the associated loop. An ordered region inside
                # an inner DO runs once per inner iteration -> undefined
                # behaviour, and in practice the ordering guarantee is lost.
                violations.append(
                    "%s:%d: !$OMP ORDERED is nested inside %d inner DO loop(s) of the "
                    "worksharing loop at %s:%d, so it executes more than once per "
                    "iteration (undefined behaviour)"
                    % (rel(dfile), dline, do_depth - 1, rel(loop[0]), loop[1])
                )
            elif not ORDERED_CLAUSE.search(" " + loop[2]):
                violations.append(
                    "%s:%d: !$OMP ORDERED binds to the loop at %s:%d "
                    "('!$OMP %s') which lacks the ORDERED clause"
                    % (rel(dfile), dline, rel(loop[0]), loop[1], loop[2])
                )

        elif END_LOOP.match(directive):
            loop = None
            do_depth = 0

        elif BEGIN_LOOP.match(directive):
            loop = (dfile, dline, directive)
            do_depth = 0
            if mode == "default" and ORDERED_CLAUSE.search(" " + directive):
                violations.append(
                    "%s:%d: ORDERED clause on '!$OMP %s' is active in the DEFAULT "
                    "build and would serialize it"
                    % (rel(dfile), dline, directive)
                )

    return n_ordered


HINT = """
A worksharing loop containing an !$OMP ORDERED region must itself carry the
ORDERED clause, and only on the __openmp_reproducible path:

    #if defined(__openmp_reproducible)
    !$OMP DO ORDERED
    #else
    !$OMP DO
    #endif

See https://github.com/FESOM/fesom2/issues/444 and
    https://github.com/FESOM/fesom2/issues/528
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-dir", required=True, help="repository root")
    ap.add_argument("--compiler", required=True, help="Fortran compiler")
    ap.add_argument("--preprocess-flags", default="-cpp -E",
                    help="flags that make the compiler preprocess only")
    ap.add_argument("--mode", required=True, choices=["reproducible", "default"])
    ap.add_argument("--extra-defines", default="",
                    help="additional -D flags, e.g. to reach code behind __oifs")
    args = ap.parse_args()

    src_dir = os.path.join(args.source_dir, "src")
    flags = args.preprocess_flags.replace(";", " ").split()
    defines = ["-D_OPENMP=201511"]
    if args.mode == "reproducible":
        defines.append("-D__openmp_reproducible")
    defines += args.extra_defines.replace(";", " ").split()

    sources = sorted(
        os.path.join(src_dir, f) for f in os.listdir(src_dir) if f.endswith(".F90")
    )

    violations, scanned, n_ordered = [], 0, 0
    for src in sources:
        # Only files that mention ORDERED at all can violate either invariant.
        with open(src, "r", errors="replace") as fh:
            if "ORDERED" not in fh.read().upper():
                continue
        scanned += 1
        text = preprocess(args.compiler, flags, defines, src_dir, src)
        n_ordered += scan(src, text, args.mode, args.source_dir, violations)

    if violations:
        print("OpenMP ORDERED check FAILED for mode=%s: %d violation(s)"
              % (args.mode, len(violations)), file=sys.stderr)
        for v in violations:
            print("  " + v, file=sys.stderr)
        print(HINT, file=sys.stderr)
        return 1

    print("OpenMP ORDERED check passed (mode=%s): %d file(s) scanned, "
          "%d ORDERED region(s) active" % (args.mode, scanned, n_ordered))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except RuntimeError as exc:
        print("OpenMP ORDERED check ERROR: %s" % exc, file=sys.stderr)
        sys.exit(2)
