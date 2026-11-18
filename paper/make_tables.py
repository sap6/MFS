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
    ref = d["ref_E"]                                    # the finest level, J = 12
    rows = []
    for name, label in (("tail", "one"), ("dense", "two")):
        if rows:
            rows.append("    \\midrule")
        for J, t, unk, nrow, opt in d[name]:
            J = int(J)
            err = "reference" if (name == "dense" and J == int(d["ref_info"][0])) else sci(rel(d[f"{name}_E{J}"], ref))
            rows.append(f"    {label} & {J} & {int(unk) // 2} & {int(nrow)} & {tm(t)} & {err} & {sci(opt)} \\\\")
    write("cube_mfs", rows)                             # both depth ratios, all levels (not in the paper)
    jref = int(d["ref_info"][0])
    main = []                                           # the table of the paper: two depths per octave, J = 4 to the reference
    for J, t, unk, nrow, opt in d["dense"]:
        J = int(J)
        if J >= 4:
            err = "reference" if J == jref else sci(rel(d[f"dense_E{J}"], ref))
            main.append(f"    {J} & {int(unk) // 2} & {int(nrow)} & {tm(t)} & {err} & {sci(opt)} \\\\")
    write("cube_mfs_short", main)
    print("    one depth per octave:", [(int(r[0]), int(r[2]) // 2, tm(r[1]), f"{rel(d[f'tail_E{int(r[0])}'], ref):.2e}") for r in d["tail"]])
    dn = d["dense"]
    ch = [rel(d[f"dense_E{int(dn[i][0])}"], d[f"dense_E{int(dn[i + 1][0])}"]) for i in range(len(dn) - 1)]
    print("    dense family, change to the next level:", " ".join(f"{c:.2e}" for c in ch), " ratios:",
          " ".join(f"{a / b:.2f}" for a, b in zip(ch[:-1], ch[1:])))
    for name in ("notail", "plain"):
        print(f"    {name}:", [(int(r[0]), int(r[2]), tm(r[1]), f"{rel(d[f'{name}_E{int(r[0])}'], ref):.2e}") for r in d[name]])
    rows = []
    bem = sorted((np.load(f) for f in glob.glob(os.path.join(V, "bem_cube_h*_q4.npz"))), key=lambda b: int(b["ndof"]))
    for b in bem:
        if int(b["ndof"]) > 1000:
            rows.append(f"    quasi-uniform & 4/6 & {int(b['ndof'])} & {sci(rel(b['E'], ref))} & {tm(float(b['t_solve']))} & -- \\\\")
    for q, tag, keep in (("4/6", "", (2592, 4608, 10368)), ("6/8", "_q6_8", (4608, 7200, 10368, 14112, 28800)),
                         ("8/10", "_q8_10", (7200, 10368, 14112, 18432, 28800))):
        rows.append("    \\midrule")
        full = {int(b["ndof"]): b for b in (np.load(f) for f in glob.glob(os.path.join(V, f"bem_cube_graded_n*_g2{tag}.npz")))}
        sym = {int(b["ndof"]): b for b in (np.load(f) for f in glob.glob(os.path.join(V, f"bem_cube_sym_n*_g2{tag}.npz")))}
        for N in sorted(set(full) | set(sym)):
            tf = tm(float(full[N]["t_solve"])) if N in full else "--"
            ts = tm(float(sym[N]["t_solve"])) if N in sym else "--"
            E = sym[N]["E"] if N in sym else full[N]["E"]
            if N in keep:
                rows.append(f"    graded & {q} & {N} & {sci(rel(E, ref))} & {tf} & {ts} \\\\")
            extra = f" (full mesh: {rel(full[N]['E'], ref):.2e})" if N in full and N in sym else ""
            print(f"    graded {q} dofs={N} err={rel(E, ref):.2e}{extra} t_full={tf} t_sym={ts}"
                  + (f" sigma_ext={float(full[N]['ext']):.6f}" if N in full else ""))
    # (the table of Bempp with g = 2 only is no longer in the paper; the rows are still printed above)
    fine = sorted((np.load(f) for f in glob.glob(os.path.join(V, "bem_cube_graded_n*_g2_q8_10.npz"))), key=lambda b: int(b["ndof"]))
    from scipy.optimize import brentq
    (n1, E1), (n2, E2), (n3, E3) = ((float(b["n"]), b["E"]) for b in fine[-3:])
    q_ = np.linalg.norm(E1 - E2) / np.linalg.norm(E2 - E3)
    p = brentq(lambda p: (n1 ** -p - n2 ** -p) / (n2 ** -p - n3 ** -p) - q_, 0.2, 12)
    Ex = E3 + (E3 - E2) * n3 ** -p / (n2 ** -p - n3 ** -p)
    s2, s3 = (float(b["ext"]) for b in fine[-2:])
    print(f"    graded Bempp 8/10 extrapolated (rate n^-{p:.2f}): error {rel(Ex, ref):.2e}, sigma_ext "
          f"{s3 + (s3 - s2) * n3 ** -p / (n2 ** -p - n3 ** -p):.6f}")


def cube_bem_env():
    """Bempp with symmetry at its best: the re-timed lower envelope of the scan over grading,
    quadrature and mesh size (../validation/bem_grading_scan.py)."""
    f = os.path.join(V, "bem_grading_scan_k1_final.npz")
    if not os.path.exists(f):
        print("--- cube_bem_env: no final scan data")
        return []
    env = [r for r in np.load(f)["rows"] if r[6] < 3e-3]
    write("cube_bem_env", [f"    {r[0]:g} & {int(r[1])}/{int(r[2])} & {int(r[4])} & {tm(r[5])} & {sci(r[6])} \\\\" for r in env])
    return env


def summary(env=None):
    """The table of selected runs of the smooth bodies and the almond in the cost section
    (the cube has its own tables: cube_mfs_short, cube_ends, cube_bem_env)."""
    rows = []
    d = np.load(os.path.join(V, "fig_cost_sphere.npz"))
    first = "Sphere"
    for u, t, e in d["mfs"]:
        if int(u) in (1200,):
            rows.append(f"    {first} & MFS & {int(u)} & {tm(t)} & {sci(e)} \\\\")
            first = ""
    u, t, e = d["bem"][-1][:3]
    rows.append(f"     & Bempp & {int(u)} & {tm(t)} & {sci(e)} \\\\")
    rows.append("    \\midrule")
    d = np.load(os.path.join(V, "fig_cost_torus.npz"))
    first = "Torus"
    for r in d["mfs"]:
        if int(r[0]) in (21120,):
            rows.append(f"    {first} & MFS & {int(r[0])} & {tm(r[1])} & {sci(r[2])} \\\\")
            first = ""
    r = d["bem"][-1]
    rows.append(f"     & Bempp & {int(r[0])} & {tm(r[1])} & {sci(r[2])} \\\\")
    rows.append("    \\midrule")
    d = np.load(os.path.join(V, "fig_cost_almond.npz"))
    ref = d[f"mfs_E{len(d['mfs']) - 1}"]
    names = ("3{:}0", "4{:}1", "5{:}2", "6{:}3")
    first = "Almond"
    for i in (3,):
        r = d["mfs"][i]
        rows.append(f"    {first} & MFS ${names[i]}$ & {int(r[1])} & {tm(r[2])} & {sci(rel(d[f'mfs_E{i}'], ref))} \\\\")
        first = ""
    i = len(d["bem"]) - 1
    r = d["bem"][i]
    rows.append(f"     & Bempp & {int(r[0])} & {tm(r[1])} & {sci(rel(d[f'bem_E{i}'], ref))} \\\\")
    write("summary", rows)


def cube_cases():
    """Cube beyond the standard case (../validation/cube_cases.py): selected runs of the lower
    envelopes of the MFS and of Bempp with symmetry, about a factor of two apart in error."""
    rows = []
    for case, label in (("k1o", "One wavelength, general wave"), ("k2", "Two wavelengths"), ("k3", "Three wavelengths")):
        f = os.path.join(V, f"cube_cases_{case}.npz")
        if not os.path.exists(f):
            print(f"--- cube_cases: no data for {case}")
            continue
        d = np.load(f)
        if rows:
            rows.append("    \\midrule")
        first = label
        mfs = sorted(d["mfs"], key=lambda r: r[5])
        env, emin = [], np.inf
        for r in mfs:                                   # MFS runs that no faster MFS run beats
            if r[6] < 0.97 * emin:
                env.append(r)
                emin = r[6]
        for name, runs, fmt in (("MFS", env, lambda r: f"MFS, $J = {int(r[1])}$ & {int(r[3])}"),
                                ("Bempp", list(d["env"]), lambda r: f"Bempp, $g = {r[0]:g}$, quadrature {int(r[1])}/{int(r[2])} & {int(r[4])}")):
            pick, last = [], np.inf
            for r in runs:
                if r[6] < 4e-3 and r[6] < last / 1.9:
                    pick.append(r)
                    last = r[6]
            if len(runs) and (not pick or pick[-1][6] != runs[-1][6]):
                pick.append(runs[-1])
            for r in pick:
                rows.append(f"    {first} & {fmt(r)} & {tm(r[5])} & {sci(r[6])} \\\\")
                first = ""
    write("cube_cases", rows)


def cube_ends():
    """Both ends of the comparison on the cube, for the four cases: the time each code needs for an error
    of about 1e-3 (its fastest run with at most 1.2e-3), the most accurate run of Bempp, and the first level of the
    MFS that is more accurate than that."""
    d = np.load(os.path.join(V, "fig_cost_cube.npz"))
    ref = d["ref_E"]
    jref = int(d["ref_info"][0])
    std_mfs = [(r[1], rel(d[f"dense_E{int(r[0])}"], ref), int(r[0])) for r in d["dense"] if int(r[0]) >= 4 and int(r[0]) != jref]
    std_bem = [(r[5], r[6]) for r in np.load(os.path.join(V, "bem_grading_scan_k1_final.npz"))["all"]]
    cases = [("Standard wave", std_mfs, std_bem)]
    for case, label in (("k1o", "General wave"), ("k2", "Two wavelengths"), ("k3", "Three wavelengths")):
        c = np.load(os.path.join(V, f"cube_cases_{case}.npz"))
        cases.append((label, [(r[5], r[6], int(r[1])) for r in c["mfs"]], [(r[5], r[6]) for r in c["bem"]]))
    rows = []
    for label, mfs, bem in cases:
        tb = min(t for t, e in bem if e <= 1.2e-3)                    # "an error of about 1e-3": at most 1.2e-3
        tm_ = min(t for t, e, J in mfs if e <= 1.2e-3)
        bt, be = min(bem, key=lambda r: r[1])
        better = [r for r in mfs if r[1] < be]
        mt, me, mJ = min(better, key=lambda r: r[0])
        rows.append(f"    {label} & {tm(tb)} & {tm(tm_)} & {sci(be)} & {tm(bt)} & {sci(me)} & {tm(mt)} \\\\")
        print(f"    cube_ends {label}: 1e-3: Bempp {tb:.1f} s, MFS {tm_:.1f} s (ratio {tm_ / tb:.1f}); Bempp best {be:.2e} in {bt:.1f} s; MFS J = {mJ}: {me:.2e} in {mt:.1f} s")
    write("cube_ends", rows)


if __name__ == "__main__":
    sphere()
    torus()
    almond()
    cube()
    summary(cube_bem_env())
    cube_cases()
    cube_ends()
