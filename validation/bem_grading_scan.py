"""Bempp on the cube at its best: scan of the grading exponent g together
with the quadrature orders and the mesh size, using the mirror symmetry
(bem_cube_sym.py). NOTES.md items 58, 59.

    python3 bem_grading_scan.py k1          # k = 2 pi, standard wave
    python3 bem_grading_scan.py k1o         # k = 2 pi, oblique wave
    python3 bem_grading_scan.py k2          # k = 4 pi, standard wave
    python3 bem_grading_scan.py k3          # k = 6 pi, standard wave
    python3 bem_grading_scan.py report      # lower envelopes of all cases
    python3 bem_grading_scan.py retime k1   # re-time the runs near the envelope on a quiet machine (best of two);
                                            # writes bem_grading_scan_k1_final.npz, the envelope used by the figures

Rows (g, regular, singular, n, dofs, time, error against the MFS reference)
are appended to bem_grading_scan_<case>.npz; finished rows are skipped.
"""
import os
import sys

import numpy as np

GS = (1.25, 1.5, 1.6, 1.75, 1.9, 2.0, 2.25)
QUAD = ((4, 6, (12, 16, 20, 24, 28, 32, 40)), (6, 8, (12, 16, 20, 24, 28, 32, 40)), (8, 10, (16, 20, 24, 28, 32, 40)))


def reference(case):
    """Wavenumber, reference field and wave (None: the standard wave) of a case."""
    if case == "k1":
        return 2 * np.pi, np.load("fig_cost_cube.npz")["ref_E"], None
    if case == "k1o":
        from cube_mfs import wave
        return 2 * np.pi, np.load("cube_cases_k1o_ref.npz")["E"], wave(60.0, 25.0, 35.0)
    return {"k2": 4 * np.pi, "k3": 6 * np.pi}[case], np.load(f"cube_cases_{case}_ref.npz")["E"], None   # cube_cases.py <case> ref


def envelope(rows):
    """Lower envelope in (time, error): the runs that no faster run beats."""
    best, emin = [], np.inf
    for r in sorted(rows, key=lambda r: r[5]):
        if r[6] < 0.97 * emin:
            best.append(r)
            emin = r[6]
    return best


def scan(case, gs=GS):
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
        for g in gs:
            for n in ns:
                if (g, reg, n) in done:
                    continue
                E, N, sizes, ts, grid = b.solve(n, g, *args)
                err = np.linalg.norm(E - ref) / np.linalg.norm(ref)
                rows.append((g, reg, sing, n, N, sum(ts), err))
                print(f"{case} g = {g:4.2f} quadrature {reg}/{sing} n = {n:2d} dofs {N:5d} {sum(ts):6.1f} s  {err:.2e}", flush=True)
                np.savez(out, rows=np.array(rows))


def retime(case):
    """The scan runs next to other work, so its times are only indicative.
    Every run that is not beaten by a run at most 0.75 times as long is
    repeated twice here, and the envelope is taken from the faster time."""
    import bem_cube_sym as b
    import bempp_cl.api as bempp
    k, ref, wave = reference(case)
    b.K = k
    pts = np.load("bem_cube_h0.05.npz")["pts"]
    rows = [tuple(r) for r in np.load(f"bem_grading_scan_{case}.npz")["rows"]]
    cand = [r for r in rows if not any(q[5] <= 0.75 * r[5] and q[6] <= r[6] for q in rows)]
    args = (pts,) if wave is None else (pts, wave[0], wave[1])
    out = []
    for reg in sorted({int(r[1]) for r in cand}):
        sing = reg + 2
        bempp.GLOBAL_PARAMETERS.quadrature.regular = reg
        bempp.GLOBAL_PARAMETERS.quadrature.singular = sing
        b.solve(4, 1.0, *args)
        for r in sorted((r for r in cand if int(r[1]) == reg), key=lambda r: r[5]):
            ts = []
            for rep in range(2):
                E, N, sizes, t4, grid = b.solve(int(r[3]), r[0], *args)
                ts.append(sum(t4))
            err = np.linalg.norm(E - ref) / np.linalg.norm(ref)
            out.append((r[0], reg, sing, r[3], N, min(ts), err))
            print(f"{case} g = {r[0]:4.2f} quadrature {reg}/{sing} n = {int(r[3]):2d} dofs {N:5d} {min(ts):6.1f} s (scan: {r[5]:.1f} s)  {err:.2e}", flush=True)
    env = envelope(out)
    np.savez(f"bem_grading_scan_{case}_final.npz", rows=np.array(env), all=np.array(out))
    print("final envelope:")
    for r in env:
        print(f"   {r[5]:6.1f} s  {r[6]:.2e}   g = {r[0]:4.2f}  {int(r[1])}/{int(r[2])}  dofs {int(r[4])}")


def report():
    for case in ("k1", "k1o", "k2", "k3"):
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
    elif sys.argv[1] == "retime":
        retime(sys.argv[2])
    else:
        scan(sys.argv[1], tuple(float(v) for v in sys.argv[2].split(",")) if len(sys.argv) > 2 else GS)
