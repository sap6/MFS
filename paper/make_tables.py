"""Writes the bodies of the tables of the paper from the data files in
../validation, so that no number is typed by hand.

    python3 make_tables.py      # writes tables/*.tex and prints the numbers quoted in the text
"""
import glob
import os

import numpy as np

V = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "validation")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tables")
os.makedirs(OUT, exist_ok=True)


def sci(x):
    e = int(np.floor(np.log10(x)))
    m = x / 10.0 ** e
    if round(m, 1) >= 10:
        m, e = m / 10, e + 1
    return f"${m:.1f}\\times10^{{{e}}}$"


def tm(t):
    return f"{t:.2f}" if t < 0.995 else f"{t:.1f}" if t < 9.95 else f"{t:.0f}"


def rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


def write(name, rows):
    open(os.path.join(OUT, name + ".tex"), "w").write("\n".join(rows) + "\n")
    print(f"--- {name}")
    for r in rows:
        print("   ", r)


def sphere():
    d = np.load(os.path.join(V, "fig_cost_sphere.npz"))
    rows = [f"    MFS   & {int(u)} & {tm(t)} & {sci(e)} \\\\" for u, t, e in d["mfs"] if int(u) in (300, 600, 1200, 2400, 3600)]
    rows.append("    \\midrule")
    rows += [f"    Bempp & {int(u)} & {tm(t)} & {sci(e)} \\\\" for u, t, e in d["bem"]]
    write("sphere", rows)


def torus():
    d = np.load(os.path.join(V, "fig_cost_torus.npz"))
    rows = [f"    MFS, mode by mode & {int(r[0])} & {tm(r[1])} & {sci(r[2])} \\\\" for r in d["mfs"]
            if int(r[0]) in (1440, 3840, 8064, 13824, 21120)]
    rows.append("    \\midrule")
    rows += [f"    MFS, full matrix  & {int(r[0])} & {tm(r[1])} & {sci(r[2])} \\\\" for r in d["dense"] if int(r[0]) in (1440, 3840)]
    rows.append("    \\midrule")
    rows += [f"    Bempp & {int(r[0])} & {tm(r[1])} & {sci(r[2])} \\\\" for r in d["bem"]]
    write("torus", rows)
    print("    reference:", d["ref"], " Bempp sigma_ext:", d["bem"][:, 3])


def almond():
    d = np.load(os.path.join(V, "fig_cost_almond.npz"))
    names = ("3{:}0", "4{:}1", "5{:}2", "6{:}3", "8{:}5")
    ref = d[f"mfs_E{len(d['mfs']) - 1}"]
    rows = [f"    MFS ${names[i]}$ & {int(r[1])} & {tm(r[2])} & {sci(rel(d[f'mfs_E{i}'], ref))} \\\\"
            for i, r in enumerate(d["mfs"][:-1])]
    rows.append("    \\midrule")
    rows += [f"    Bempp & {int(r[0])} & {tm(r[1])} & {sci(rel(d[f'bem_E{i}'], ref))} \\\\" for i, r in enumerate(d["bem"])]
    write("almond", rows)
    print("    finest layout:", d["mfs"][-1])


def cube():
    d = np.load(os.path.join(V, "fig_cost_cube.npz"))
    ref = d["tail_E9"]
    rows = []
    for J, t, unk, nrow, opt in d["tail"]:
        J = int(J)
        err = "reference" if J == 9 else sci(rel(d[f"tail_E{J}"], ref))
        rows.append(f"    {J} & {int(unk) // 2} & {int(nrow)} & {tm(t)} & {err} & {sci(opt)} \\\\")
    write("cube_mfs", rows)
    for name in ("notail", "plain"):
        print(f"    {name}:", [(int(r[0]), int(r[2]), tm(r[1]), f"{rel(d[f'{name}_E{int(r[0])}'], ref):.2e}") for r in d[name]])
    rows = []
    bem = sorted((np.load(f) for f in glob.glob(os.path.join(V, "bem_cube_h*_q4.npz"))), key=lambda b: int(b["ndof"]))
    for b in bem:
        if int(b["ndof"]) > 1000:
            rows.append(f"    quasi-uniform & 4/6 & {int(b['ndof'])} & {tm(float(b['t_solve']))} & {sci(rel(b['E'], ref))} \\\\")
    for q, tag, keep in (("4/6", "", (2592, 4608, 10368)), ("6/8", "_q6_8", (4608, 7200, 10368, 14112)),
                         ("8/10", "_q8_10", (10368, 14112))):
        rows.append("    \\midrule")
        gr = sorted((np.load(f) for f in glob.glob(os.path.join(V, f"bem_cube_graded_n*_g2{tag}.npz"))), key=lambda b: int(b["ndof"]))
        for b in gr:
            if int(b["ndof"]) in keep:
                rows.append(f"    graded & {q} & {int(b['ndof'])} & {tm(float(b['t_solve']))} & {sci(rel(b['E'], ref))} \\\\")
            print(f"    graded {q} n={int(b['n'])} dofs={int(b['ndof'])} t={float(b['t_solve']):.1f} (assembly {float(b['t_assembly']):.1f}, "
                  f"LU {float(b['t_lu']):.1f}) err={rel(b['E'], ref):.2e} sigma_ext={float(b['ext']):.5f}")
    write("cube_bem", rows)


def cube_oblique():
    r = np.load(os.path.join(V, "cube_oblique.npz"))
    ref = r["E9"]
    rows = []
    for J in (4, 5, 6, 7, 8):
        t, unk = r[f"info{J}"][:2]
        rows.append(f"    MFS, $J = {J}$ & {int(unk)} & {tm(t)} & {sci(rel(r[f'E{J}'], ref))} \\\\")
    print("    oblique reference J=9:", r["info9"])
    rows.append("    \\midrule")
    files = sorted(glob.glob(os.path.join(V, "bem_cube_graded_*_obl.npz")), key=lambda f: ("_q" in f, int(np.load(f)["ndof"])))
    for f in files:
        b = np.load(f)
        if np.allclose(b["d"], r["d"]) and np.allclose(b["p"], r["p"]):
            q = os.path.basename(f).split("_q")[1].split("_obl")[0].replace("_", "/") if "_q" in f else "4/6"
            rows.append(f"    Bempp, quadrature {q} & {int(b['ndof'])} & {tm(float(b['t_solve']))} & {sci(rel(b['E'], ref))} \\\\")
    write("cube_oblique", rows)


if __name__ == "__main__":
    sphere()
    torus()
    almond()
    cube()
    cube_oblique()
