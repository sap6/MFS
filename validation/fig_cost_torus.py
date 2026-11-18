"""Error against wall time for the torus R = 1, a = 0.5 (k = 2 pi, wave along
the axis): MFS with dipoles on the inner torus a_s = 0.5 a, and Bempp (EFIE,
RWG). No closed-form solution exists; the error is the relative difference of
the scattered field at 200 points on |x| = 2.5 from the finest MFS run
(13 440 dipoles, boundary residual printed by `run`).

    python3 fig_cost_torus.py run    # solve and time on this machine; writes fig_cost_torus.npz
    python3 fig_cost_torus.py        # table and figure from the saved data

"MFS" solves mode by mode about the axis (uses the rotational symmetry of the
body, any incident wave); "MFS, full matrix" is the same least-squares problem
without using the symmetry; Bempp uses the rotational symmetry of its
structured mesh (bem_mirror.Cyclic). Times are for the solve and the field
at the test points. Run alone on the machine.
"""
import sys
import time

import numpy as np

from pec_torus import Torus, circle, cross_sections, plane_wave, rel, residual, scatter, solve, solve_dense, sphere_rule, wave

K = 2 * np.pi
A = 0.5
LAYOUTS = ((6, 20), (8, 28), (12, 40), (16, 52), (20, 64), (24, 80), (28, 96), (36, 128), (44, 160))
DENSE = ((6, 20), (8, 28), (12, 40), (16, 52), (20, 64))
BEM_N = (12, 20, 28, 40, 56, 80, 112, 160)
REF = (56, 240)
STYLE = {"MFS": ("#2a78d6", "o", "-", "MFS, dipoles on an inner torus (mode by mode)"),
         "MFS dense": ("#1baf7a", "s", "--", "the same, full matrix (no symmetry used)"),
         "Bempp": ("#eb6834", "D", "-.", "Bempp (EFIE, RWG), rotational symmetry of the mesh")}
INK, MUTED = "#1f1f1e", "#6b6a64"


def pts():
    return 2.5 * sphere_rule(10, 20)[0]


def run():
    import subprocess
    tor = Torus(1.0, A)
    d, p = wave()
    x = pts()
    ref = scatter(tor, K, d, p, *REF, a_s=0.5 * A)
    Eref = ref.field(x)
    print(f"reference: {ref.info['unknowns']} unknowns, residual {ref.info['bc']:.1e}, optical {ref.info['opt']:.1e}, "
          f"sigma_ext {ref.info['sigma_ext']:.10f}, t={ref.info['t']:.1f}s", flush=True)
    out = {"ref": np.array([ref.info["unknowns"], ref.info["bc"], ref.info["opt"], ref.info["sigma_ext"]])}
    rows = []
    rhs = lambda y: -plane_wave(y, K, d, p)
    for ns, nphi in LAYOUTS:
        t0 = time.time()
        rho, z, _ = circle(ns, 1.0, 0.5 * A)
        s = solve(tor, K, rhs, rho, z, nphi)
        E = s.field(x)
        t = time.time() - t0                                    # solve and field at the test points; the checks below are not timed
        ext, sca = cross_sections(s, d, p)
        s.info.update(bc=residual(tor, s, rhs), opt=abs(ext - sca) / sca)
        rows.append((s.info["unknowns"], t, rel(E, Eref), s.info["bc"], s.info["opt"], s.info["alpha_norm"]))
        print(f"MFS       unknowns={rows[-1][0]:6d} time={t:7.2f}s error={rows[-1][2]:.2e} residual={rows[-1][3]:.1e} optical={rows[-1][4]:.1e}", flush=True)
    out["mfs"] = np.array(rows)
    rows = []
    for ns, nphi in DENSE:
        t0 = time.time()
        rho, z, _ = circle(ns, 1.0, 0.5 * A)
        s = solve_dense(tor, K, rhs, rho, z, nphi)
        E = s.field(x)
        t = time.time() - t0
        rows.append((s.info["unknowns"], t, rel(E, Eref)))
        print(f"MFS dense unknowns={rows[-1][0]:6d} time={t:7.2f}s error={rows[-1][2]:.2e}", flush=True)
    out["dense"] = np.array(rows)
    np.savez("fig_cost_torus.npz", **out)
    subprocess.run([sys.executable, "-W", "ignore", "bem_torus.py", ",".join(str(n) for n in BEM_N), str(A), "--rot"], check=True)
    rows = []
    for n in BEM_N:
        b = np.load(f"bem_torus_a{A:g}_n{n}_rot.npz")
        rows.append((int(b["ndof"]), float(b["t"]), rel(b["E"], Eref), float(b["ext"])))
        print(f"Bempp     unknowns={rows[-1][0]:6d} time={rows[-1][1]:7.2f}s error={rows[-1][2]:.2e}", flush=True)
    out["bem"] = np.array(rows)
    np.savez("fig_cost_torus.npz", **out)


def report():
    d = np.load("fig_cost_torus.npz")
    r = d["ref"]
    print(f"reference MFS: {int(r[0])} unknowns, boundary residual {r[1]:.1e}, optical mismatch {r[2]:.1e}, sigma_ext {r[3]:.8f}")
    print("method      unknowns   time[s]   error vs finest MFS   (MFS: boundary residual, optical mismatch)")
    for row in d["mfs"]:
        print(f"MFS        {int(row[0]):8d} {row[1]:9.2f} {row[2]:12.2e}          {row[3]:.1e}  {row[4]:.1e}")
    for row in d["dense"]:
        print(f"MFS dense  {int(row[0]):8d} {row[1]:9.2f} {row[2]:12.2e}")
    for row in d["bem"]:
        print(f"Bempp      {int(row[0]):8d} {row[1]:9.2f} {row[2]:12.2e}          sigma_ext {row[3]:.5f}")
    return {"MFS": d["mfs"][:, 1:3], "MFS dense": d["dense"][:, 1:3], "Bempp": d["bem"][:, 1:3]}


def plot(curves, fname="fig_cost_torus"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for name, c in curves.items():
        col, mk, ls, lab = STYLE[name]
        ax.loglog(c[:, 0], c[:, 1], color=col, marker=mk, ls=ls, lw=1.6, ms=5.5, label=lab)
    ax.set_xlabel("wall time [s]", color=INK)
    ax.set_ylabel("relative error of the scattered field", color=INK)
    ax.set_title("PEC torus (R = 1, a = 0.5, k = 2$\\pi$): error against time", color=INK, fontsize=10.5)
    ax.grid(True, which="major", color="#d9d8d2", lw=0.6)
    ax.tick_params(colors=MUTED)
    for s in ax.spines.values():
        s.set_color(MUTED)
    ax.legend(frameon=False, fontsize=8.5)
    fig.tight_layout()
    fig.savefig(fname + ".pdf")
    fig.savefig(fname + ".png", dpi=200)
    print("wrote", fname + ".pdf/.png")


if __name__ == "__main__":
    if "run" in sys.argv:
        run()
    plot(report())
