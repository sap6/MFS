"""
Tables 2, 4, 5 and C.1 of the paper: writes their bodies to paper/tables/*.tex
from data/, so that no number in a table is typed by hand.

    python figures/tables.py
"""
import os

import numpy as np
from style import HERE, data, envelope, rel

OUT = os.path.join(HERE, "..", "paper", "tables")


def sci(x):
    """x as m × 10^e in LaTeX, one decimal"""
    m, e = f"{x:.1e}".split("e")
    return f"${m}\\times10^{{{int(e)}}}$"


def tm(t):
    """A time in seconds, with two significant digits below ten"""
    return f"{t:.2f}" if t < 0.995 else f"{t:.1f}" if t < 9.95 else f"{t:.0f}"


def write(name, rows):
    os.makedirs(OUT, exist_ok=True)
    open(os.path.join(OUT, name + ".tex"), "w").write("\n".join(rows) + "\n")
    print(f"--- paper/tables/{name}.tex")
    print("\n".join(rows))


# Table 2: the cube, levels J = 4 to the reference J = 12 (unknowns and equations per class)
d = np.load(data("cube_convergence.npz"))
ref, jref = d["ref_E"], int(d["ref_info"][0])
write("cube_mfs_short", [f"    {int(J)} & {int(unk) // 2} & {int(rows)} & {tm(t)} & "
                         f"{'reference' if int(J) == jref else sci(rel(d[f'dense_E{int(J)}'], ref))} & {sci(opt)} \\\\"
                         for J, t, unk, rows, opt in d["dense"] if J >= 4])

# Table 4: the sphere (400 and 1200 sources), the torus and the almond (layout 6:3), with Bempp's largest mesh
# Columns: Shape, Method, whole-body unknowns, unknowns solved for the incident wave, wall time, error
s, t, a = np.load(data("cost_sphere.npz")), np.load(data("cost_torus.npz")), np.load(data("cost_almond.npz"))
aref = a[f"mfs_E{len(a['mfs']) - 1}"]
nsrc, ncls, talm = int(a["mfs"][3, 0]), int(a["mfs"][3, 1]), a["mfs"][3, 2]
alm_tot = 12 * nsrc + 2 * (ncls - 3 * nsrc) * 9 // 5     # four classes: two even in z (5 tail/row), two odd (4 tail/row)
nbem = int(a["bem"][-1, 0])
row = lambda sh, m, tot, sol, tm_, err: f"    {sh} & {m} & {int(tot)} & {int(sol)} & {tm(tm_)} & {sci(err)} \\\\"
write("summary", [*[row("Sphere" if int(r[0]) == 1200 else "", "MFS", r[0], r[0], r[1], r[2]) for r in s["mfs"] if int(r[0]) in (1200, 3600)],
                  row("", "Bempp", s["bem"][-1, 0], s["bem"][-1, 0], s["bem"][-1, 1], s["bem"][-1, 2]),
                  "    \\midrule", *[row("Torus", "MFS", r[0], 2 * r[0] // 160, r[1], r[2]) for r in t["mfs"] if int(r[0]) == 21120],
                  row("", "Bempp", t["bem"][-1, 0], 2 * t["bem"][-1, 0] // 448, t["bem"][-1, 1], t["bem"][-1, 2]), "    \\midrule",
                  row("Almond", "MFS $6{:}3$", alm_tot, ncls, talm, rel(a["mfs_E3"], aref)),
                  row("", "Bempp", nbem, nbem // 4, a["bem"][-1, 1], rel(a[f"bem_E{len(a['bem']) - 1}"], aref))])

# Table 5: both ends of the cube comparison. Bempp's most accurate run and the first MFS level
# below it; then wall time for an error of about 1e-3 (at most 1.2e-3)
rows = []
labels = (("k1", r"Standard ($k=2\pi$)"), ("k1o", r"General ($k=2\pi$)"),
          ("k2", r"Standard ($k=4\pi$)"), ("k3", r"Standard ($k=6\pi$)"))
for case, label in labels:
    c = np.load(data(f"cube_cases_{case}.npz"))
    mfs, bem = [(r[5], r[6]) for r in c["mfs"]], [(r[5], r[6]) for r in c["bem"]]
    tb, tm_ = min(t for t, e in bem if e <= 1.2e-3), min(t for t, e in mfs if e <= 1.2e-3)
    bt, be = envelope(bem)[-1]
    mt, me = min(((t, e) for t, e in mfs if e < be))
    rows.append(f"    {label} & {sci(be)} & {tm(bt)} & {sci(me)} & {tm(mt)} & {tm(tb)} & {tm(tm_)} \\\\")
write("cube_ends", rows)

# Table C.1: the fastest Bempp runs on the cube for each accuracy (the lower envelope), errors below 3e-3
write("cube_bem_env", [f"    {r[0]:g} & {int(r[1])}/{int(r[2])} & {int(r[4])} & {tm(r[5])} & {sci(r[6])} \\\\"
                       for r in np.load(data("cube_cases_k1.npz"))["env"] if r[6] < 3e-3])
