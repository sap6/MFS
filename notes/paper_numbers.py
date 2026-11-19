"""Every number of Sections 5 and 7 that depends on the cube runs, computed from data/ (python notes/paper_numbers.py)"""
import os

import numpy as np
from scipy.optimize import brentq

D = os.path.join(os.path.dirname(__file__), "..", "data", "")
rel = lambda a, b: float(np.linalg.norm(a - b) / np.linalg.norm(b))


def envelope(points, f=0.97):
    out, best = [], np.inf
    for t, e in sorted(points):
        if e < f * best:
            out.append((t, e))
            best = e
    return np.array(out)


def T(c, e):
    le, lt = np.log(c[:, 1])[::-1], np.log(c[:, 0])[::-1]
    return np.exp(np.interp(np.log(e), le, lt)) if c[:, 1].min() <= e <= c[:, 1].max() else np.nan


def cross(a, b):
    lo, hi = max(a[:, 1].min(), b[:, 1].min()), min(a[:, 1].max(), b[:, 1].max())
    es = np.exp(np.linspace(np.log(lo), np.log(hi), 6001))
    f = np.array([np.log(T(a, e) / T(b, e)) for e in es])
    return [np.exp(brentq(lambda x: np.log(T(a, np.exp(x)) / T(b, np.exp(x))), np.log(es[i]), np.log(es[i + 1])))
            for i in range(len(es) - 1) if f[i] * f[i + 1] < 0]


# Section 5.4: the standard wave
c = np.load(D + "cube_convergence.npz")
ref = c["ref_E"]
dense = {int(r[0]): r for r in c["dense"]}
err = {J: rel(c[f"dense_E{J}"], ref) for J in dense if J != 12}
print("== Section 5.4")
print("errors J = 3..11:", {J: f"{e:.2e}" for J, e in err.items()})
print("time J = 4..12:", {J: round(dense[J][1], 1) for J in dense})
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
print("one depth: errors", {J: f"{x:.2e}" for J, x in eo.items()},
      "unknowns ratio", [round(one[J][2] / dense[J][2], 2) for J in one if J in dense],
      "time ratio", [round(one[J][1] / dense[J][1], 2) for J in one if J in dense],
      "error ratio", [round(eo[J] / err[J], 2) for J in one if J in err])
nt = {int(r[0]): rel(c[f"notail_E{int(r[0])}"], ref) for r in c["notail"]}
print("no tails:", {J: f"{x:.1e}" for J, x in nt.items()}, "tails reduce J = 8 by", round(nt[8] / err[8]),
      f"decay like d_min^{np.polyfit(np.log(0.2 * 2.0 ** -np.array(sorted(nt))), np.log([nt[J] for J in sorted(nt)]), 1)[0]:.2f}")
pl = {int(r[0]): rel(c[f"plain_E{int(r[0])}"], ref) for r in c["plain"]}
print("plain layout:", {J: f"{x:.1e}" for J, x in pl.items()})

# Section 5.5 and 7: the four cases
print("== Sections 5.5 and 7")
for case in ("k1", "k1o", "k2", "k3"):
    d = np.load(D + f"cube_cases_{case}.npz")
    r = np.load(D + f"cube_cases_{case}_ref.npz")["info"]
    two = {int(x[1]): x for x in d["mfs"] if x[0] == 2}
    one = {int(x[1]): x for x in d["mfs"] if x[0] == 1}
    bem = d["env"][np.argsort(d["env"][:, 5])][:, [5, 6]]
    print(f"-- {case}: reference J = {int(r[0])}, {r[1]:.0f} s, optical {r[4]:.1e}, ext {r[5]:.6f}")
    print("   two depths (J, time, error):", [(J, round(x[5], 1), f"{x[6]:.2e}") for J, x in two.items()])
    print("   factors per level:", [round(two[J][6] / two[J + 1][6], 2) for J in list(two)[:-1]])
    print("   one depth  (J, time, error):", [(J, round(x[5], 1), f"{x[6]:.2e}") for J, x in one.items()])
    m2 = envelope([(x[5], x[6]) for x in two.values()])
    mall = envelope([(x[5], x[6]) for x in d["mfs"]])
    print(f"   crossover: two depths {['%.1e' % x for x in cross(m2, bem)]}, both layouts {['%.1e' % x for x in cross(mall, bem)]}")
    for e in (1e-3, 5e-4, 3e-4, 2e-4, 1e-4):
        print(f"   at {e:.0e}: MFS both {T(mall, e):6.1f} s, two depths {T(m2, e):6.1f} s, Bempp {T(bem, e):6.1f} s")
    bt, be = envelope([(x[5], x[6]) for x in d["bem"]])[-1]
    mt, me = min((x[5], x[6]) for x in d["mfs"] if x[6] < be)
    print(f"   Bempp best {be:.1e} in {bt:.0f} s; MFS first below: {me:.1e} in {mt:.0f} s")
