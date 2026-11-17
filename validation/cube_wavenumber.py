"""The cube layout at another wavenumber (side = k / (2 pi) wavelengths).
The interior collocation panels are scaled with the wavelength and each line
gets 4 more modes per extra wavelength; everything else is the family of
cube_mfs.py. (Without the finer interior collocation the family stalls near
3e-4 at twice the wavenumber; a finer base grid, 12 or 16 points per face,
makes the coefficients grow and is not needed.)

    python3 cube_wavenumber.py 2 explore      # k = 2 * 2 pi: which scaling is needed at J = 6
    python3 cube_wavenumber.py 2              # convergence table with the adopted scaling
    python3 cube_wavenumber.py 2 bem          # Bempp with the mirror symmetry on graded meshes, against the MFS field
"""
import sys

import numpy as np

import pec_cube_lightning3d as m

S = float(sys.argv[1]) if len(sys.argv) > 1 else 2.0
m.K = 2 * np.pi * S                       # before any solver object is built
import cube_mfs as c                      # noqa: E402


def rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


def params(J, family, n_base, hmax, dn, depth0=0.2):
    kw = c.dense(J) if family == "dense" else dict(zip(("kw", "hmax"), c.balanced(J)))
    kw = dict(kw)
    kw["kw"] = dict(kw["kw"], n_base=n_base, n0=kw["kw"]["n0"] + dn, depth0=depth0)
    kw["hmax"] = hmax
    return kw


def run(J, family, n_base, hmax, dn, pts, depth0=0.2):
    kw = params(J, family, n_base, hmax, dn, depth0)
    E, t, info = c.run(J, pts=pts, family=family, **kw)
    return E, t, info


if __name__ == "__main__":
    pts = np.load("bem_cube_h0.05.npz")["pts"]
    c.run(2, pts=pts)
    if "explore" in sys.argv:
        out = {}
        for n_base, hmax, dn in ((6, 0.15, 0), (8, 0.1, 0), (8, 0.1, 4), (12, 0.075, 4), (12, 0.075, 8), (16, 0.05, 8)):
            E, t, info = run(6, "balanced", n_base, hmax, dn, pts)
            out[(n_base, hmax, dn)] = E
            print(f"k = {S:g} * 2 pi, J = 6, base {n_base} / {hmax}, modes + {dn}: unknowns {info['per_class']} rows {info['rows']} t = {t:.1f} s "
                  f"optical {info['optical']:.1e} sigma_ext {info['ext']:.5f} |c| {info['coef_norm']:.1f}", flush=True)
        keys = list(out)
        for a in keys[:-1]:
            print(f"  {a} against {keys[-1]}: {rel(out[a], out[keys[-1]]):.2e}")
    elif "bem" in sys.argv:
        import bem_cube_sym as b
        import bempp_cl.api as bempp
        b.K = m.K
        d = np.load(f"cube_wavenumber_k{S:g}.npz")
        ref = d["dense_E10"] if "dense_E10" in d else d["dense_E8"]
        rows = []
        for reg, sing, ns in ((4, 6, (12, 16, 24, 32, 40)), (6, 8, (16, 24, 32, 40)), (8, 10, (16, 24, 32, 40))):
            bempp.GLOBAL_PARAMETERS.quadrature.regular = reg
            bempp.GLOBAL_PARAMETERS.quadrature.singular = sing
            b.solve(4, 1.0, pts)                           # compilation for this quadrature
            for n in ns:
                E, N, sizes, ts, grid = b.solve(n, 2.0, pts)
                rows.append((reg, sing, n, N, sum(ts), rel(E, ref)))
                print(f"k = {S:g} * 2 pi, Bempp with symmetry, graded n = {n}, quadrature {reg}/{sing}: dofs {N}, {sum(ts):.1f} s, "
                      f"difference from the MFS reference: {rel(E, ref):.2e}", flush=True)
                np.savez(f"cube_wavenumber_k{S:g}_bem.npz", rows=np.array(rows))
    else:
        dn = int(round(4 * (S - 1)))
        base = lambda J: (4, 0.25 / S) if J <= 5 else (6, 0.15 / S)     # interior panels scaled with the wavelength
        depth0 = 0.2
        if S > 2:                                          # beyond two wavelengths the whole base layout follows the wavelength
            base = lambda J: (int(4 * S), 0.25 / S) if J <= 5 else (int(6 * S), 0.15 / S)
            depth0 = 0.2 / S
        res = {}
        JREF = 10 if S <= 2 else 8                         # memory: one class of the reference must stay below about 5 GB
        for fam, Js in (("balanced", range(4, 9)), ("dense", range(5, JREF + 1))):
            for J in Js:
                E, t, info = run(J, fam, *base(J), dn, pts, depth0)
                res[(fam, J)] = (E, t, info)
                print(f"{fam:8s} J = {J:2d}: unknowns/class {max(info['per_class'])} rows {info['rows']} t = {t:5.1f} s optical {info['optical']:.1e} "
                      f"sigma_ext {info['ext']:.6f} |c| {info['coef_norm']:.1f}", flush=True)
        ref = res[("dense", JREF)][0]
        if S <= 2:
            E, t, info = run(8, "dense", 8, 0.15 / S, dn + 4, pts)
            print(f"check: dense J = 8 with base 8 and {dn + 4} more modes differs from the family's J = 8 by {rel(E, res[('dense', 8)][0]):.2e} "
                  f"(optical {info['optical']:.1e})")
        print(f"errors against dense J = {JREF}:")
        for (fam, J), (E, t, info) in res.items():
            print(f"  {fam:8s} J = {J:2d}  t = {t:5.1f} s  {rel(E, ref):.2e}")
        np.savez(f"cube_wavenumber_k{S:g}.npz", pts=pts, k=m.K, **{f"{fam}_E{J}": v[0] for (fam, J), v in res.items()},
                 **{f"{fam}_info{J}": np.array([v[1], max(v[2]["per_class"]), v[2]["rows"], v[2]["optical"], v[2]["ext"]]) for (fam, J), v in res.items()})
