"""Every number of Sections 3, 5, 6, 7 and Appendix C that depends on stored data (python notes/paper_numbers.py)"""
import os

import numpy as np
from scipy.optimize import brentq

D = os.path.join(os.path.dirname(__file__), "..", "data", "")
rel = lambda a, b: float(np.linalg.norm(a - b) / np.linalg.norm(b))


def envelope(points, f=0.97):
    out, best = [], np.inf
    for t, e in sorted(points):
        if e < f * best:
            out.append((float(t), float(e)))
            best = float(e)
    return np.array(out)


def T(c, e):
    le, lt = np.log(c[:, 1])[::-1], np.log(c[:, 0])[::-1]
    return float(np.exp(np.interp(np.log(e), le, lt))) if c[:, 1].min() <= e <= c[:, 1].max() else np.nan


def T_step(points, e):
    ts = [float(t) for t, err in points if err <= e]
    return min(ts) if ts else np.nan


def cross(a, b):
    lo, hi = max(a[:, 1].min(), b[:, 1].min()), min(a[:, 1].max(), b[:, 1].max())
    if lo >= hi:
        return []
    es = np.exp(np.linspace(np.log(lo), np.log(hi), 6001))
    f = np.array([np.log(T(a, e) / T(b, e)) for e in es])
    return [float(np.exp(brentq(lambda x: np.log(T(a, np.exp(x)) / T(b, np.exp(x))), np.log(es[i]), np.log(es[i + 1]))))
            for i in range(len(es) - 1) if f[i] * f[i + 1] < 0]


def cross_staircase(m_pts, b_pts, emax=3e-3):
    bem_min = min(e for _, e in b_pts)
    thresholds = sorted({float(e) for _, e in list(m_pts) + list(b_pts) if bem_min <= e <= emax}, reverse=True)
    verdicts = [(e, T_step(m_pts, e), T_step(b_pts, e)) for e in thresholds]
    mfs_wins = [e for e, tm, tb in verdicts if tm <= tb]
    if not mfs_wins:
        return None, None, []
    e_first = max(mfs_wins)
    bem_below = [e for e, tm, tb in verdicts if e < e_first and tb < tm]
    e_strict = max([e for e in mfs_wins if not any(eb <= e for eb in bem_below)], default=None)
    return e_first, e_strict, bem_below


def local_orders(N_err):
    return [round(float(-np.log(N_err[i + 1, 1] / N_err[i, 1]) / np.log(N_err[i + 1, 0] / N_err[i, 0])), 2)
            for i in range(len(N_err) - 1)]


# Smooth bodies and Almond (Sections 3, 6, 7.1)
print("== Sections 3, 6 & 7.1")
import json  # noqa: E402
sc_json = [r for r in json.load(open(D + "sphere_convergence.json")) if abs(r["k"] - 2 * np.pi) < 1e-9]
print("Sphere Fig 1 (r=1.5, R_s=0.5):", [(r["N"], f"mid={r['mid']:.2e}", f"bc={r['bc']:.2e}") for r in sc_json if r["sweep"] == "N" and r["R"] == 0.5 and r["N"] in (100, 400, 1200)])
cs, ct, ca = np.load(D + "cost_sphere.npz"), np.load(D + "cost_torus.npz"), np.load(D + "cost_almond.npz")
aref = ca[f"mfs_E{len(ca['mfs']) - 1}"]
print("Sphere MFS (r=2.0):", [(int(r[0]), round(float(r[1]), 2), f"{r[2]:.2e}") for r in cs["mfs"] if int(r[0]) in (300, 1200, 3600)],
      "Bempp all:", [(int(r[0]), round(float(r[1]), 1), f"{r[2]:.2e}") for r in cs["bem"]])
print("Torus ref:", (int(ct["ref"][0]), f"bc={ct['ref'][1]:.2e}", f"opt={ct['ref'][2]:.2e}", f"ext={ct['ref'][3]:.6f}"))
print("Torus MFS:", [(int(r[0]), round(float(r[1]), 2), f"err={r[2]:.2e}", f"bc={r[3]:.2e}", f"opt={r[4]:.2e}")
                     for r in ct["mfs"] if int(r[0]) in (360, 1440, 21120)],
      "Bempp finest:", (int(ct["bem"][-1, 0]), round(float(ct["bem"][-1, 1]), 1), f"{ct['bem'][-1, 2]:.2e}"))
alm_mfs = [(int(r[1]), round(float(r[2]), 2), f"{rel(ca[f'mfs_E{i}'], aref):.2e}", f"opt={r[3]:.2e}") for i, r in enumerate(ca["mfs"][:-1])]
alm_bem = [(int(r[0]), round(float(r[1]), 2), f"{rel(ca[f'bem_E{i}'], aref):.2e}") for i, r in enumerate(ca["bem"])]
print("Almond 1.19 GHz MFS:", alm_mfs, "ref 8:5:", (int(ca["mfs"][-1, 1]), round(float(ca["mfs"][-1, 2]), 2), f"opt={ca['mfs'][-1, 3]:.2e}"), "Bempp:", alm_bem)
if os.path.exists(D + "cost_almond35.npz"):
    ca35 = np.load(D + "cost_almond35.npz")
    aref35 = ca35[f"mfs_E{len(ca35['mfs']) - 1}"]
    print("Almond 3.5 GHz MFS:", [(int(r[1]), round(float(r[2]), 2), f"{rel(ca35[f'mfs_E{i}'], aref35):.2e}") for i, r in enumerate(ca35["mfs"][:-1])],
          "Bempp:", [(int(r[0]), round(float(r[1]), 2), f"{rel(ca35[f'bem_E{i}'], aref35):.2e}") for i, r in enumerate(ca35["bem"])])
snf = np.load(D + "sphere_nearfield.npz")
print("Sphere nearfield max(mfs400), max(mfs1200):", f"{snf['mfs400'].max():.2e}", f"{snf['mfs1200'].max():.2e}",
      "bem9750 at delta ~0.03, 0.01, 1e-3:", [(round(float(d_), 4), f"{e_:.2e}") for d_, e_ in zip(snf["delta"], snf["bem9750"]) if any(abs(np.log10(d_) - x) < 0.15 for x in (-1.5, -2.0, -3.0))])

# Section 5.4: the standard wave
c = np.load(D + "cube_convergence.npz")
ref = c["ref_E"]
dense = {int(r[0]): r for r in c["dense"]}
err = {J: rel(c[f"dense_E{J}"], ref) for J in dense if J != 12}
print("\n== Section 5.4")
print("errors J = 3..11:", {J: f"{e:.2e}" for J, e in err.items()})
print("time J = 3..12:", {J: round(float(dense[J][1]), 1) for J in dense})
print("factor per level J = 4..10:", [round(err[J] / err[J + 1], 2) for J in range(4, 10)])
N = np.array([8 * dense[J][2] / 2 for J in range(4, 11)])            # unknowns of the whole cube (the entry is two classes)
e = np.array([err[J] for J in range(4, 11)])
b, a = np.polyfit(np.sqrt(N), np.log(e), 1)
print(f"fit eps = {np.exp(a):.3f} exp({b:.4f} sqrt(N_u)), largest deviation {np.max(np.abs(np.exp(a + b * np.sqrt(N)) / e - 1)):.1%}")
dmin = 0.2 * 2.0 ** -np.arange(4, 11)
print(f"decay like d_min^{np.polyfit(np.log(dmin), np.log(e), 1)[0]:.2f}")
print(f"J = 12: {int(dense[12][2]) // 2} unknowns per class, optical {dense[12][4]:.1e}; J = 11 against J = 12: {err[11]:.1e}; "
      f"error of the reference about {err[11] / 1.5:.1e}")
one = {int(r[0]): r for r in c["tail"]}
eo = {J: rel(c[f"tail_E{J}"], ref) for J in one}
print("one depth J=4..9: errors", {J: f"{eo[J]:.2e}" for J in range(4, 10)},
      "unknowns ratio", [round(float(one[J][2] / dense[J][2]), 2) for J in range(4, 10)],
      "time ratio", [round(float(one[J][1] / dense[J][1]), 2) for J in range(4, 10)],
      "error ratio J=4..7", [round(eo[J] / err[J], 2) for J in range(4, 8)])
nt = {int(r[0]): rel(c[f"notail_E{int(r[0])}"], ref) for r in c["notail"]}
print("no tails:", {J: f"{x:.1e}" for J, x in nt.items()}, "tails reduce J = 8 by", round(nt[8] / err[8]),
      f"decay like d_min^{np.polyfit(np.log(0.2 * 2.0 ** -np.array(sorted(nt))), np.log([nt[J] for J in sorted(nt)]), 1)[0]:.2f}")
pl = {int(r[0]): rel(c[f"plain_E{int(r[0])}"], ref) for r in c["plain"]}
print("plain layout:", {J: f"{x:.1e}" for J, x in pl.items()})
cu = np.load(D + "cube_bempp_uniform.npz")
print("Bempp quasi-uniform:", [(float(h), int(n), round(float(t), 1), f"{rel(E, ref):.2e}") for n, t, E, h in zip(cu["ndof"], cu["t"], cu["E"], cu["h"])])
c1 = np.load(D + "cube_cases_k1.npz")
best_mesh = np.array(sorted({int(N_): min(r[6] for r in c1["bem"] if int(r[4]) == int(N_)) for N_ in c1["bem"][:, 4]}.items()))
mfs_mesh = np.array([(4 * dense[J][2], err[J]) for J in range(4, 11)])
print("Orders on k1: Bempp =", local_orders(best_mesh), "MFS =", local_orders(mfs_mesh))

# Sections 5.5, 7.2 and Appendix C: the four cases
print("\n== Sections 5.5, 7.2 & Appendix C")
k1_two = {int(x[1]): x for x in c1["mfs"] if x[0] == 2}
for case in ("k1", "k1o", "k2", "k3"):
    d = np.load(D + f"cube_cases_{case}.npz")
    r = np.load(D + f"cube_cases_{case}_ref.npz")["info"]
    two = {int(x[1]): x for x in d["mfs"] if x[0] == 2}
    one = {int(x[1]): x for x in d["mfs"] if x[0] == 1}
    m_pts = [(float(x[5]), float(x[6])) for x in d["mfs"]]
    m2_pts = [(float(x[5]), float(x[6])) for x in two.values()]
    b_pts = [(float(x[5]), float(x[6])) for x in d["bem"]]
    bem = d["env"][np.argsort(d["env"][:, 5])][:, [5, 6]]
    jref = int(r[0])
    est_ref = two[jref - 1][6] / 1.5
    print(f"-- {case}: ref J={jref}, {r[1]:.0f} s, opt={r[4]:.1e}, ext={r[5]:.6f}, est ref err={est_ref:.1e} ({est_ref / two[8][6]:.0%} of J=8)")
    print("   two depths (J, +unk%, time, t/t_k1, error):",
          [(J, f"{x[2] / k1_two[J][2] - 1:+.0%}", round(float(x[5]), 1), round(float(x[5] / k1_two[J][5]), 2), f"{x[6]:.2e}")
           for J, x in two.items() if J in k1_two])
    print("   factors per level:", [round(float(two[J][6] / two[J + 1][6]), 2) for J in list(two)[:-1]])
    print("   one depth  (J, time, error):", [(J, round(float(x[5]), 1), f"{x[6]:.2e}") for J, x in one.items()])
    for N_chk in (28800, 41472):
        eb = [x[6] for x in d["bem"] if int(x[4]) == N_chk]
        if eb:
            print(f"   Bempp best on N={N_chk}: {min(eb):.2e}")
    m2, mall = envelope(m2_pts), envelope(m_pts)
    ef_all, es_all, gaps_all = cross_staircase(m_pts, b_pts)
    ef_2, es_2, gaps_2 = cross_staircase(m2_pts, b_pts)
    print(f"   crossover interp: two depths {['%.2e' % x for x in cross(m2, bem)]}, both {['%.2e' % x for x in cross(mall, bem)]}")
    print(f"   crossover step  : two depths first={ef_2:.2e} strict={es_2:.2e} gaps={[f'{g:.2e}' for g in gaps_2]}; "
          f"both first={ef_all:.2e} strict={es_all:.2e} gaps={[f'{g:.2e}' for g in gaps_all]}")
    for e_ in (1e-3, 5e-4, 3e-4, 2e-4, 1.5e-4, 1e-4):
        print(f"   at {e_:.1e}: interp (both {T(mall, e_):5.1f}s, 2dep {T(m2, e_):5.1f}s, BEM {T(bem, e_):5.1f}s) | "
              f"step (both {T_step(m_pts, e_):5.1f}s, 2dep {T_step(m2_pts, e_):5.1f}s, BEM {T_step(b_pts, e_):5.1f}s)")
    tb1, tm1 = min(t for t, e_ in b_pts if e_ <= 1.2e-3), min(t for t, e_ in m_pts if e_ <= 1.2e-3)
    bt, be = envelope(b_pts)[-1]
    mt, me, mnd, mJ = min((float(x[5]), float(x[6]), int(x[0]), int(x[1])) for x in d["mfs"] if x[6] < be)
    print(f"   Table 5 (<=1.2e-3): Bempp {tb1:.2f} s vs MFS {tm1:.2f} s (MFS/Bempp={tm1 / tb1:.2f}x); "
          f"most accurate: Bempp {be:.2e} in {bt:.1f} s vs MFS ({mnd},{mJ}) {me:.2e} in {mt:.1f} s ({mt / bt:.2f}x)")
    env_sub = [row for row in d["env"] if row[6] < 3e-3]
    print("   App C env (<3e-3) gradings:", sorted({float(row[0]) for row in env_sub}))

# Section 5.6: the regular tetrahedron
from mfs import tetrahedron as _tet  # noqa: E402
tc = np.load(D + "tetrahedron_convergence.npz")
tref = tc["ref_E"]
tdense = {int(r[0]): r for r in tc["dense"]}
tNu = {J: sum(_tet.Tetrahedron(J, classes=[]).unknowns.values()) for J in tdense}
terr = {J: rel(tc[f"dense_E{J}"], tref) for J in tdense if J != 10}
print("\n== Section 5.6 (Regular tetrahedron)")
print("errors J = 3..9:", {J: f"{e:.2e}" for J, e in terr.items()})
print("unk/cls & N_u J = 3..10:", {J: (int(tdense[J][2]), tNu[J]) for J in tdense})
print("time J = 3..10:", {J: round(float(tdense[J][1]), 2) for J in tdense})
print("optical J = 3..10:", {J: f"{tdense[J][4]:.2e}" for J in tdense})
print("ext J = 8..10:", {J: f"{tdense[J][5]:.6f}" for J in (8, 9, 10)})
print("factor per level J = 3..8:", [round(terr[J] / terr[J + 1], 2) for J in range(3, 9)])
tN = np.array([tNu[J] for J in range(3, 9)], float)
te = np.array([terr[J] for J in range(3, 9)])
tb, ta = np.polyfit(np.sqrt(tN), np.log(te), 1)
print(f"fit eps = {np.exp(ta):.3f} exp({tb:.4f} sqrt(N_u)), largest factor {np.exp(np.max(np.abs(ta + tb * np.sqrt(tN) - np.log(te)))):.2f} "
      f"(largest rel dev {np.max(np.abs(np.exp(ta + tb * np.sqrt(tN)) / te - 1)):.1%})")
tdmin = 0.14 * 2.0 ** -np.arange(3, 9)
print(f"decay like d_min^{np.polyfit(np.log(tdmin), np.log(te), 1)[0]:.2f}")
tone = {int(r[0]): r for r in tc["one"]}
teo = {J: rel(tc[f"one_E{J}"], tref) for J in tone}
print("one depth J=3..8: errors", {J: f"{teo[J]:.2e}" for J, r in tone.items()},
      "times", {J: round(float(r[1]), 2) for J, r in tone.items()})
tnt = {int(r[0]): rel(tc[f"notail_E{int(r[0])}"], tref) for r in tc["notail"]}
print("no tails:", {J: f"{x:.2e}" for J, x in tnt.items()}, "tails reduce J = 8 by", round(tnt[8] / terr[8]),
      f"decay like d_min^{np.polyfit(np.log(0.14 * 2.0 ** -np.array(sorted(tnt))), np.log([tnt[J] for J in sorted(tnt)]), 1)[0]:.2f}")
tpl = {int(r[0]): (sum(_tet.Tetrahedron(None, classes=[], n_base=int(r[0])).unknowns.values()), rel(tc[f"plain_E{int(r[0])}"], tref))
       for r in tc["plain"]}
print("plain layout (N_u, err):", {nb: (nu, f"{e_:.2e}") for nb, (nu, e_) in tpl.items()})
print("Bempp graded (n, q, N_u, t, err):",
      [(int(r[0]), f"{int(r[1])}/{int(r[1])+2}", int(r[2]), round(float(r[3]), 1), f"{r[4]:.2e}") for r in tc["bem"]])
tbem_env = envelope([(float(r[3]), float(r[4])) for r in tc["bem"]])
tm2_env = envelope([(float(tdense[J][1]), terr[J]) for J in range(3, 10)])
tmall_env = envelope([(float(tdense[J][1]), terr[J]) for J in range(3, 10)] + [(float(tone[J][1]), teo[J]) for J in tone])
print("Bempp time envelope:", [(round(float(t), 1), f"{e:.2e}") for t, e in tbem_env])
print("Tet crossover interp: two depths", [f"{x:.2e}" for x in cross(tm2_env, tbem_env)],
      "both", [f"{x:.2e}" for x in cross(tmall_env, tbem_env)])

# Section 5.6: zero-symmetry scalene convex tetrahedron
if os.path.exists(D + "scalene_tetrahedron.npz"):
    import itertools
    sc = np.load(D + "scalene_tetrahedron.npz")
    sv = _tet.SCALENE_VERTS
    pr_sc = _tet.ScaleneTetrahedron(2, verts=sv)
    elens = [float(np.linalg.norm(sv[b] - sv[a])) for a, b in itertools.combinations(range(4), 2)]
    dih_deg, gammas = [], []
    for a, b in itertools.combinations(range(4), 2):
        fc, fd = [f for f in range(4) if f not in (a, b)]
        dih = float(np.arccos(np.clip(-pr_sc.fnorms[fc] @ pr_sc.fnorms[fd], -1.0, 1.0)))
        dih_deg.append(np.degrees(dih))
        gammas.append(np.pi / (2 * np.pi - dih))
    nu1s = [e[0][0] for e in _tet.SCALENE_EXP]
    ri = sc["ref_info"]
    sm2 = sc["mfs"][sc["mfs"][:, 0] == 2]
    sp = np.polyfit(np.sqrt(sm2[:, 2]), np.log(sm2[:, 5]), 1)
    print("\n== Section 5.6 (Scalene zero-symmetry tetrahedron)")
    print(f"edge lengths: [{min(elens):.3f}, {max(elens):.3f}], dihedral deg: [{min(dih_deg):.1f}, {max(dih_deg):.1f}], "
          f"gamma_e: [{min(gammas):.3f}, {max(gammas):.3f}], nu_1: [{min(nu1s):.3f}, {max(nu1s):.3f}]")
    print(f"ref J={int(ri[0])}: N_u={int(ri[1])}, t={ri[3]:.1f}s, opt={ri[4]:.2e}")
    print("m=2 (J, N_u, t, err, opt):", [(int(r[1]), int(r[2]), round(float(r[4]), 1), f"{r[5]:.2e}", f"{r[6]:.2e}") for r in sm2])
    print("m=1 (J, N_u, t, err, opt):", [(int(r[1]), int(r[2]), round(float(r[4]), 1), f"{r[5]:.2e}", f"{r[6]:.2e}") for r in sc["mfs"] if r[0] == 1])
    print(f"fit eps = {np.exp(sp[1]):.2f} exp({sp[0]:.4f} sqrt(N_u)), factors: {[round(float(sm2[i,5]/sm2[i+1,5]), 2) for i in range(len(sm2)-1)]}")
    print("Bempp n=40:", [(int(r[0]), f"{int(r[1])}/{int(r[1])+2}", int(r[2]), round(float(r[3]), 1), f"{r[4]:.2e}") for r in sc["bem"] if int(r[0]) == 40])

# Section 8: matrix-free block-QR preconditioned LSMR, non-convex L-block, and k = 8*pi cube
if os.path.exists(D + "iterative_check.npz"):
    ic = np.load(D + "iterative_check.npz")
    ex = {int(r[0]): r for r in ic.get("scalene_extra", np.zeros((0, 6)))}
    print("\n== Section 8 (Matrix-free Block-QR LSMR)")
    for r in ic["cube"]:
        print(f"Cube J={int(r[0])} k={int(r[1])}*2pi: A={r[4]:.1f} MiB, P={r[5]:.2f} MiB ({r[4]/r[5]:.1f}x), "
              f"Blk itn={int(r[6])}, tot={r[7]+r[8]:.1f}s, rel_E={r[9]:.2e} | Diag itn={int(r[10])}, rel_E={r[11]:.2e}")
    for r in ic["scalene"]:
        J = int(r[0])
        it1e4 = int(ex[J][5]) if (J in ex and len(ex[J]) > 5) else int(r[5])
        opt_s = f", opt={ex[J][1]:.2e}, vs_J7={ex[J][2]:.2e}" if J in ex else ""
        print(f"Scalene J={J} (N={int(r[2])}): A={r[3]:.1f} MiB ({r[3]/1024:.2f} GiB), P={r[4]:.1f} MiB ({r[3]/r[4]:.1f}x), "
              f"itn={int(r[5])} (itn@1e-4={it1e4}), tot={(r[6]+r[7])/60:.1f} min, rel_vs_dir={r[8]:.2e}{opt_s}")

if os.path.exists(D + "nonconvex_check.npz"):
    nc = np.load(D + "nonconvex_check.npz")
    print("\n== Section 8 (Non-convex 3D L-block, k = 2*pi)")
    print(f"inner copy min err: {nc['mfs_inner'][:, 5].min():.3f}, full J=5 err vs J=6: {nc['mfs_full'][-1, 5]:.2e} "
          f"(opt J=6: {nc['ref_info'][7]:.2e}), Bempp best diff vs J=6: {nc['bem'][:, 9].min():.2e}")

if os.path.exists(D + "cube_k4_check.npz"):
    k4 = np.load(D + "cube_k4_check.npz")
    m1_8 = [r for r in k4["mfs"] if int(r[0]) == 1 and int(r[1]) == 8][0]
    m2_8 = [r for r in k4["mfs"] if int(r[0]) == 2 and int(r[1]) == 8][0]
    b48 = [b for b in k4["bem"] if int(b[3]) == 48][0]
    ng6 = k4["ng6_J8"]
    print("\n== Section 8 (Four-wavelength cube, k = 8*pi)")
    print(f"m=1 J=8: err={m1_8[6]:.2e} in {m1_8[5]:.1f}s | m=2 J=8: err={m2_8[6]:.2e} in {m2_8[5]:.1f}s (ng=6 shift {ng6[1]:.2e}) | "
          f"Bempp n=48 6/8: diff={b48[6]:.2e} in {b48[5]:.1f}s")



