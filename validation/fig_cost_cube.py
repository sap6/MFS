"""Error against wall time for the unit cube (k = 2 pi): MFS with edge sources
and singular tail, the same without tail, plain MFS (uniform sources on an
inner cube), and Bempp (EFIE, RWG, dense LU) on quasi-uniform meshes and on meshes graded
towards the edges. Two MFS families: one source depth per octave (fast, stalls
near 1e-4) and two (converges further; its finest run is the reference).

    python3 fig_cost_cube.py run     # solve and time everything except Bempp; writes fig_cost_cube.npz
    python3 fig_cost_cube.py         # table (and figure, if matplotlib is present) from the saved data

Bempp fields and times come from bem_optical.py (bem_cube_h*_q4.npz) and
bem_cube_graded.py (bem_cube_graded_n*_g2.npz), run on the same machine. All times include assembly, solve and evaluation at the 400
test points on |x| = 2; the MFS uses the two mirror symmetries of this
incidence, Bempp does not. Run alone on the machine.
"""
import glob
import sys
import time

import numpy as np

import cube_mfs
from cube_mfs import COL, balanced  # noqa: F401  (re-exported: other scripts import them from here)

J_REF = 10          # reference: the family with two depths per octave at this level
J_FAST = 9          # deepest level of the family with one depth per octave


def run_mfs(J, pts, tail=9, n_base=None, family="balanced"):
    """Standard wave with the three-mirror solver (cube_mfs.py). Time covers
    set-up, assembly, the factorisations and the field at the test points."""
    kw = {} if n_base is None else dict(kw=dict(q=0.5, gamma=2.0, kinds="em", n_base=n_base), hmax=0.1)
    E, t, info = cube_mfs.run(J, pts=pts, tail=tail, family=family, **kw)
    return E, t, info["unknowns"], info["rows"], info["optical"]


def run():
    pts = np.load("bem_cube_h0.05.npz")["pts"]
    run_mfs(2, pts)                                    # warm-up (compilation, threads)
    run_mfs(2, pts, family="dense")
    out = {"pts": pts}
    for name, Js, kw in (("tail", range(3, J_FAST + 1), {}), ("dense", range(3, J_REF + 1), dict(family="dense")),
                         ("notail", (4, 6, 8), dict(tail=0)), ("plain", (4, 6, 8, 10), None)):
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
    ref = d[f"dense_E{J_REF}"]
    bem = sorted((np.load(f) for f in glob.glob("bem_cube_h*_q4.npz")), key=lambda b: int(b["ndof"]))
    bem = [b for b in bem if int(b["ndof"]) > 1000]
    # Bempp extrapolated from its two finest meshes, error ~ N^(-2/3) (h^(4/3), the edge-limited rate)
    f1, f2 = (float(b["ndof"]) ** (-2 / 3) for b in bem[-2:])
    ext = bem[-1]["E"] + (bem[-1]["E"] - bem[-2]["E"]) * f2 / (f1 - f2)
    print(f"reference: MFS with tails, two depths per octave, J = {J_REF}. It differs from extrapolated quasi-uniform Bempp by {rel(ref, ext):.2e}")
    curves = {}
    print("method        size   unknowns   time[s]   error vs MFS ref   error vs extrap. Bempp   optical")
    for name, label in (("tail", "MFS+tail"), ("dense", "MFS dense"), ("notail", "MFS no tail"), ("plain", "plain MFS")):
        pts_ = []
        for J, t, unk, nrow, opt in d[name]:
            E = d[f"{name}_E{int(J)}"]
            if name == "dense" and int(J) == J_REF:
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
    # Bempp on meshes graded towards the edges (bem_cube_graded.py, grading exponent 2, default quadrature)
    graded = sorted((np.load(f) for f in glob.glob("bem_cube_graded_n*_g2.npz")), key=lambda b: int(b["ndof"]))
    pts_ = []
    for b in graded:
        pts_.append((float(b["t_solve"]), rel(b["E"], ref)))
        print(f"{'Bempp graded':12s} {int(b['n']):5d} {int(b['ndof']):9d} {float(b['t_solve']):9.1f} {rel(b['E'], ref):14.2e} "
              f"{rel(b['E'], ext):18.2e}   (assembly {float(b['t_assembly']):.0f}s, LU {float(b['t_lu']):.0f}s, "
              f"evaluation {float(b['t_eval']):.0f}s)")
    if pts_:
        curves["Bempp graded"] = np.array(pts_)
    # the same meshes with higher quadrature orders: a second cost curve (6/8) and the independent reference (8/10)
    for q in ("6_8", "8_10"):
        fine = sorted((np.load(f) for f in glob.glob(f"bem_cube_graded_n*_g2_q{q}.npz")), key=lambda b: int(b["ndof"]))
        for b in fine:
            print(f"{'graded ' + q.replace('_', '/'):12s} {int(b['n']):5d} {int(b['ndof']):9d} {float(b['t_solve']):9.1f} {rel(b['E'], ref):14.2e} "
                  f"{'':18s}   (assembly {float(b['t_assembly']):.0f}s, LU {float(b['t_lu']):.0f}s); sigma_ext {float(b['ext']):.5f}")
        if q == "6_8" and fine:
            curves["Bempp graded 6/8"] = np.array([(float(b["t_solve"]), rel(b["E"], ref)) for b in fine])
        if len(fine) >= 3:                              # three-level extrapolation, E_n = E + C n^(-p)
            from scipy.optimize import brentq
            (n1, E1), (n2, E2), (n3, E3) = ((float(b["n"]), b["E"]) for b in fine[-3:])
            r = np.linalg.norm(E1 - E2) / np.linalg.norm(E2 - E3)
            p = brentq(lambda p: (n1 ** -p - n2 ** -p) / (n2 ** -p - n3 ** -p) - r, 0.2, 12)
            Ex = E3 + (E3 - E2) * n3 ** -p / (n2 ** -p - n3 ** -p)
            s1, s2, s3 = (float(b["ext"]) for b in fine[-3:])
            print(f"   extrapolated from the last three (fitted rate n^-{p:.2f}): differs from the MFS reference by "
                  f"{rel(Ex, ref):.2e}; sigma_ext {s3 + (s3 - s2) * n3 ** -p / (n2 ** -p - n3 ** -p):.5f}")
    # Bempp with the three mirror planes of the cube (bem_cube_sym.py): the comparison on an equal footing
    best = []
    for q, tag in (("4/6", ""), ("6/8", "_q6_8"), ("8/10", "_q8_10")):
        for b in sorted((np.load(f) for f in glob.glob(f"bem_cube_sym_n*_g2{tag}.npz")), key=lambda b: int(b["ndof"])):
            e = rel(b["E"], ref)
            best.append((float(b["t_solve"]), e))
            print(f"{'sym. ' + q:12s} {int(b['n']):5d} {int(b['ndof']):9d} {float(b['t_solve']):9.1f} {e:14.2e} "
                  f"{'':18s}   (set-up {float(b['t_setup']):.1f}s, assembly {float(b['t_assembly']):.1f}s, fold and LU {float(b['t_lu']):.1f}s)")
    if best:                                           # lower envelope over the three quadrature settings
        best.sort()
        env, emin = [], np.inf
        for t, e in best:
            if e < 0.97 * emin:
                env.append((t, e))
                emin = e
        curves["Bempp sym"] = np.array(env)
    sm = d["dense"][-1]
    print(f"MFS reference J = {int(sm[0])}: {int(sm[2])} unknowns, {sm[1]:.0f} s, optical mismatch {sm[4]:.1e}")
    return curves


# Fixed colour, marker and dash per method (colour is never the only cue, so
# the figure survives greyscale printing). Hues checked for colour-blind
# separation with the dataviz palette validator.
STYLE = {"MFS+tail": ("#2a78d6", "o", "-", "MFS with tails, one depth per octave"),
         "MFS dense": ("#184f96", "h", (0, (6, 1.5)), "MFS with tails, two depths per octave"),
         "MFS no tail": ("#1baf7a", "s", "--", "MFS, edge sources, no tail"),
         "plain MFS": ("#eda100", "^", ":", "plain MFS (inner cube)"),
         "Bempp": ("#eb6834", "D", "-.", "Bempp, quasi-uniform mesh"),
         "Bempp graded": ("#e87ba4", "v", (0, (5, 1, 1, 1)), "Bempp, graded mesh, default quadrature"),
         "Bempp graded 6/8": ("#8a5cc7", "P", (0, (1, 1)), "Bempp, graded mesh, quadrature 6/8"),
         "Bempp sym": ("#9c3b12", "X", (0, (3, 1)), "Bempp, graded mesh, using the mirror symmetry")}
INK, MUTED = "#1f1f1e", "#6b6a64"


def plot(curves, fname="fig_cost_cube"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.6})
    fig, ax = plt.subplots(figsize=(6.0, 4.6))
    for key, (col, mk, ls, label) in STYLE.items():
        if key not in curves:
            continue
        c = curves[key]
        ax.loglog(c[:, 0], c[:, 1], linestyle=ls, color=col, marker=mk, ms=6, lw=1.6, mec="white", mew=0.8, label=label)
    ax.set_xlabel("wall time [s]  (4 cores)")
    ax.set_ylabel("relative error of the field on |x| = 2")
    ax.grid(True, which="major", color="#e4e3df", lw=0.6)
    ax.grid(False, which="minor")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, columnspacing=1.2)
    ax.set_title("Unit cube, side = one wavelength", loc="left", fontsize=9, color=INK)
    fig.tight_layout()
    fig.savefig(fname + ".pdf", bbox_inches="tight")
    fig.savefig(fname + ".png", dpi=200, bbox_inches="tight")
    print("wrote", fname + ".pdf/.png")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        run()
    plot(report())
