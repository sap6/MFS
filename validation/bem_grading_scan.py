"""Bempp on the cube at its best: scan of the grading exponent g together
with the quadrature orders and the mesh size, using the mirror symmetry
(bem_cube_sym.py). NOTES.md items 58, 59.

    python3 bem_grading_scan.py k1          # k = 2 pi, standard wave
    python3 bem_grading_scan.py k1o         # k = 2 pi, oblique wave
    python3 bem_grading_scan.py k2          # k = 4 pi, standard wave
    python3 bem_grading_scan.py k2o         # k = 4 pi, oblique wave
    python3 bem_grading_scan.py report      # lower envelopes of all cases

Rows (g, regular, singular, n, dofs, time, error against the MFS reference)
are appended to bem_grading_scan_<case>.npz; finished rows are skipped.
"""
import os
import sys

import numpy as np

GS = (1.25, 1.5, 1.6, 1.75, 1.9, 2.0, 2.25)
QUAD = ((4, 6, (12, 16, 20, 24, 28, 32, 40)), (6, 8, (12, 16, 20, 24, 28, 32, 40)), (8, 10, (16, 20, 24, 28, 32, 40)))


def reference(case):
    if case == "k1":
        return 2 * np.pi, np.load("fig_cost_cube.npz")["dense_E11"], None
    if case == "k1o":
        d = np.load("cube_oblique.npz")
        return 2 * np.pi, d["Eref"], (d["d"], d["p"])
    if case == "k2":
        return 4 * np.pi, np.load("cube_wavenumber_k2.npz")["dense_E10"], None
    d = np.load("cube_wavenumber_k2_obl.npz")
    return 4 * np.pi, d["dense_E9"], (d["d"], d["p"])


def envelope(rows):
    """Lower envelope in (time, error): the runs that no faster run beats."""
    best, emin = [], np.inf
    for r in sorted(rows, key=lambda r: r[5]):
        if r[6] < 0.97 * emin:
            best.append(r)
            emin = r[6]
    return best


def scan(case):
    import bem_cube_sym as b
    import bempp_cl.api as bempp
    k, ref, wave = reference(case)
    b.K = k
    pts = np.load("bem_cube_h0.05.npz")["pts"]
    out = f"bem_grading_scan_{case}.npz"
    rows = [tuple(r) for r in np.load(out)["rows"]] if os.path.exists(out) else []
    done = {(r[0], int(r[1]), int(r[3])) for r in rows}
    args = (pts,) if wave is None else (pts, wave[0], wave[1])
    for reg, sing, ns in QUAD:
        bempp.GLOBAL_PARAMETERS.quadrature.regular = reg
        bempp.GLOBAL_PARAMETERS.quadrature.singular = sing
        b.solve(4, 1.0, *args)                              # compilation for this quadrature
        for g in GS:
            for n in ns:
                if (g, reg, n) in done:
                    continue
                E, N, sizes, ts, grid = b.solve(n, g, *args)
                err = np.linalg.norm(E - ref) / np.linalg.norm(ref)
                rows.append((g, reg, sing, n, N, sum(ts), err))
                print(f"{case} g = {g:4.2f} quadrature {reg}/{sing} n = {n:2d} dofs {N:5d} {sum(ts):6.1f} s  {err:.2e}", flush=True)
                np.savez(out, rows=np.array(rows))


def report():
    for case in ("k1", "k1o", "k2", "k2o"):
        f = f"bem_grading_scan_{case}.npz"
        if not os.path.exists(f):
            continue
        rows = [tuple(r) for r in np.load(f)["rows"]]
        print(f"{case}: {len(rows)} runs. Lower envelope, all settings:")
        for r in envelope(rows):
            print(f"   {r[5]:6.1f} s  {r[6]:.2e}   g = {r[0]:4.2f}  {int(r[1])}/{int(r[2])}  dofs {int(r[4])}")
        for g in sorted({r[0] for r in rows}):
            e = envelope([r for r in rows if r[0] == g])
            print(f"   g = {g:4.2f} alone: " + "  ".join(f"{r[5]:.1f}s {r[6]:.1e}" for r in e))
        e = envelope([r for r in rows if int(r[1]) >= 6])
        print("   quadrature 6/8 and 8/10 only: " + "  ".join(f"{r[5]:.1f}s {r[6]:.1e} (g {r[0]:g})" for r in e))


if __name__ == "__main__":
    if sys.argv[1] == "report":
        report()
    else:
        scan(sys.argv[1])
