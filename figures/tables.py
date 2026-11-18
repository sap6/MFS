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
    e = int(np.floor(np.log10(x)))
    m = x / 10.0 ** e
    if round(m, 1) >= 10:
        m, e = m / 10, e + 1
    return f"${m:.1f}\\times10^{{{e}}}$"


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
s, t, a = np.load(data("cost_sphere.npz")), np.load(data("cost_torus.npz")), np.load(data("cost_almond.npz"))
aref = a[f"mfs_E{len(a['mfs']) - 1}"]
nb = len(a["bem"]) - 1
write("summary", [f"    Sphere & MFS & {int(r[0])} & {tm(r[1])} & {sci(r[2])} \\\\" if i == 0 else f"     & MFS & {int(r[0])} & {tm(r[1])} & {sci(r[2])} \\\\"
                  for i, r in enumerate(r for r in s["mfs"] if int(r[0]) in (1200, 3600))] +
      [f"     & Bempp & {int(s['bem'][-1][0])} & {tm(s['bem'][-1][1])} & {sci(s['bem'][-1][2])} \\\\", "    \\midrule"] +
      [f"    Torus & MFS & {int(r[0])} & {tm(r[1])} & {sci(r[2])} \\\\" for r in t["mfs"] if int(r[0]) == 21120] +
      [f"     & Bempp & {int(t['bem'][-1][0])} & {tm(t['bem'][-1][1])} & {sci(t['bem'][-1][2])} \\\\", "    \\midrule",
       f"    Almond & MFS $6{{:}}3$ & {int(a['mfs'][3][1])} & {tm(a['mfs'][3][2])} & {sci(rel(a['mfs_E3'], aref))} \\\\",
       f"     & Bempp & {int(a['bem'][nb][0])} & {tm(a['bem'][nb][1])} & {sci(rel(a[f'bem_E{nb}'], aref))} \\\\"])

# Table 5: both ends of the cube comparison. Time for an error of about 1e-3 (at most 1.2e-3);
# Bempp's most accurate run (the end of its lower envelope) and the first MFS level below it
rows = []
for case, label in (("k1", "Standard wave"), ("k1o", "General wave"), ("k2", "Two wavelengths"), ("k3", "Three wavelengths")):
    c = np.load(data(f"cube_cases_{case}.npz"))
    mfs, bem = [(r[5], r[6]) for r in c["mfs"]], [(r[5], r[6]) for r in c["bem"]]
    tb, tm_ = min(t for t, e in bem if e <= 1.2e-3), min(t for t, e in mfs if e <= 1.2e-3)
    bt, be = envelope(bem)[-1]
    mt, me = min(((t, e) for t, e in mfs if e < be))
    rows.append(f"    {label} & {tm(tb)} & {tm(tm_)} & {sci(be)} & {tm(bt)} & {sci(me)} & {tm(mt)} \\\\")
write("cube_ends", rows)

# Table C.1: the fastest Bempp runs on the cube for each accuracy (the lower envelope), errors below 3e-3
write("cube_bem_env", [f"    {r[0]:g} & {int(r[1])}/{int(r[2])} & {int(r[4])} & {tm(r[5])} & {sci(r[6])} \\\\"
                       for r in np.load(data("cube_cases_k1.npz"))["env"] if r[6] < 3e-3])
