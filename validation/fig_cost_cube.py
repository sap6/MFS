"""Error against wall time for the unit cube (k = 2 pi): MFS with edge sources
and singular tail, the same without tail, plain MFS (uniform sources on an
inner cube), and Bempp (EFIE, RWG, dense LU).

    python3 fig_cost_cube.py run     # solve and time everything except Bempp; writes fig_cost_cube.npz
    python3 fig_cost_cube.py         # table (and figure, if matplotlib is present) from the saved data

Bempp fields and times come from bem_optical.py (bem_cube_h*_q4.npz), run on
the same machine. All times include assembly, solve and evaluation at the 400
test points on |x| = 2; the MFS uses the two mirror symmetries of this
incidence, Bempp does not. Run alone on the machine.
"""
import glob
import sys
import time

import numpy as np

import pec_cube_lightning3d as m

FAST = dict(q=0.5, gamma=2.0, kinds="em", n_base=4)
COL = dict(col_p=3, col_q=0.5, beta=0.5, hmax=0.25)
J_REF = 9


def run_mfs(J, pts, tail=9, n_base=4):
    t0 = time.time()
    kw = dict(FAST, n_base=n_base)
    src = m.SymSources(J, tail_n=tail, **kw)
    d_min = src.depths[-1] if J >= 0 else 0.2
    lv = max(2, int(np.ceil(np.log(d_min / (2 * m.H)) / np.log(COL["col_q"]))) + 1)
    c, info = m.solve(src, lv, COL["col_p"], COL["col_q"], COL["beta"], hmax=COL["hmax"] if J >= 0 else 0.1)
    E = (src.matrix(pts) @ c).reshape(-1, 3)
    t = time.time() - t0
    ext, sca = m.cross_sections(src, c)
    return E, t, info["unknowns"], info["rows"], abs(ext - sca) / sca


def run():
    pts = np.load("bem_cube_h0.05.npz")["pts"]
    run_mfs(2, pts)                                    # warm-up (compilation, threads)
    out = {"pts": pts}
    for name, Js, kw in (("tail", range(3, J_REF + 1), {}), ("notail", (4, 6, 8, 10), dict(tail=0)),
                         ("plain", (4, 6, 8, 10), None)):
        rows = []
        for J in Js:
            if kw is None:                             # plain MFS: J plays the role of the grid size
                E, t, unk, nrow, opt = run_mfs(-1, pts, n_base=J)
            else:
                E, t, unk, nrow, opt = run_mfs(J, pts, **kw)
            rows.append((J, t, unk, nrow, opt))
            out[f"{name}_E{J}"] = E
            print(f"{name:7s} J/n={J:2d} unknowns={unk:5d} rows={nrow:6d} time={t:6.1f}s optical={opt:.1e}", flush=True)
            out[name] = np.array(rows)
            np.savez("fig_cost_cube.npz", **out)


def rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


def report():
    d = np.load("fig_cost_cube.npz")
    ref = d[f"tail_E{J_REF}"]
    bem = sorted((np.load(f) for f in glob.glob("bem_cube_h*_q4.npz")), key=lambda b: int(b["ndof"]))
    bem = [b for b in bem if int(b["ndof"]) > 1000]
    # Bempp extrapolated from its two finest meshes, error ~ N^(-2/3) (h^(4/3), the edge-limited rate)
    f1, f2 = (float(b["ndof"]) ** (-2 / 3) for b in bem[-2:])
    ext = bem[-1]["E"] + (bem[-1]["E"] - bem[-2]["E"]) * f2 / (f1 - f2)
    print(f"reference: MFS with tail, J = {J_REF}. It differs from extrapolated Bempp by {rel(ref, ext):.2e}")
    curves = {}
    print("method        size   unknowns   time[s]   error vs MFS ref   error vs extrap. Bempp   optical")
    for name, label in (("tail", "MFS+tail"), ("notail", "MFS no tail"), ("plain", "plain MFS")):
        pts_ = []
        for J, t, unk, nrow, opt in d[name]:
            E = d[f"{name}_E{int(J)}"]
            if name == "tail" and int(J) == J_REF:
                continue
            pts_.append((t, rel(E, ref)))
            print(f"{label:12s} {int(J):5d} {int(unk):9d} {t:9.1f} {rel(E, ref):14.2e} {rel(E, ext):18.2e} {opt:16.1e}")
        curves[label] = np.array(pts_)
    pts_ = []
    for b in bem:
        pts_.append((float(b["t_solve"]), rel(b["E"], ref)))
        print(f"{'Bempp':12s} {float(b['h']):5.2f} {int(b['ndof']):9d} {float(b['t_solve']):9.1f} {rel(b['E'], ref):14.2e} "
              f"{rel(b['E'], ext):18.2e}")
    curves["Bempp"] = np.array(pts_)
    return curves


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        run()
    report()
