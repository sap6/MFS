"""Tuning of the cube layout at three wavelengths (k = 6 pi), NOTES.md item 64.

    python3 cube_k3_tune.py size                 # matrix sizes of the candidate references
    python3 cube_k3_tune.py ref J                # two depths per octave at level J -> cube_k3_ref_J<J>.npz
    python3 cube_k3_tune.py try fam,J,base,hmax[,dn[,colp[,pm-ng]]] ...   # variants, errors against the best reference on disk
"""
import glob
import sys

import numpy as np

import pec_cube_lightning3d as m

S = 3
m.K = 2 * np.pi * S
import cube_mfs as c  # noqa: E402
import line_tree as lt  # noqa: E402

NODES = lt.nodes


def kwargs(J, fam, base, hmax, dn=8, colp=None):
    kw = c.dense(J) if fam == "dense" else dict(zip(("kw", "hmax"), c.balanced(J)))
    kw = dict(kw)
    kw["kw"] = dict(kw["kw"], n_base=base, n0=kw["kw"]["n0"] + dn, depth0=0.2 / S)
    kw["hmax"] = hmax
    if colp is not None:
        kw["col"] = dict(kw.get("col", c.COL), col_p=colp)
    return kw


def best_ref():
    f = sorted(glob.glob("cube_k3_ref_J*.npz"), key=lambda f: int(f.split("_J")[1].split(".")[0]))[-1]
    return f, np.load(f)["E"]


if __name__ == "__main__":
    pts = np.load("bem_cube_h0.05.npz")["pts"]
    if sys.argv[1] == "size":
        for J in (8, 9, 10, 11):
            kw = kwargs(J, "dense", 18, 0.05)
            rows, unk = c._size(2 * J, c.classes_of(m.D, m.P)[0], kw)
            print(f"dense J = {J}: rows {rows} unknowns/class about {unk}: {16.0 * rows * unk / 2 ** 30:.1f} GB per class")
        sys.exit()
    c.run(2, pts=pts, **kwargs(2, "balanced", 18, 0.05))
    c.run(2, pts=pts, family="dense", **kwargs(2, "dense", 18, 0.05))
    if sys.argv[1] == "ref":
        J = int(sys.argv[2])
        E, t, info = c.run(J, pts=pts, family="dense", **kwargs(J, "dense", 18, 0.05))
        print(f"reference dense J = {J}: unknowns/class {info['per_class']} rows {info['rows']} t = {t:.0f} s optical {info['optical']:.1e} "
              f"sigma_ext {info['ext']:.6f} |c| {info['coef_norm']:.1f}", flush=True)
        np.savez(f"cube_k3_ref_J{J}.npz", E=E, pts=pts, info=np.array([t, max(info["per_class"]), info["rows"], info["optical"], info["ext"]]))
        for f in sorted(glob.glob("cube_k3_ref_J*.npz")):
            print(f"  against {f}: {np.linalg.norm(E - np.load(f)['E']) / np.linalg.norm(E):.2e}")
    else:
        fref, ref = best_ref()
        import os
        store = dict(np.load("cube_k3_t_tries.npz")) if os.path.exists("cube_k3_t_tries.npz") else {}
        if sys.argv[1] == "diff":                                       # mutual differences of the stored fields
            keys = [k for k in store if not k.startswith("info:")]
            names = keys + ["ref"]
            F = [store[k] for k in keys] + [ref]
            for i, a in enumerate(names):
                print(f"{a:28s} " + " ".join(f"{np.linalg.norm(F[i] - F[j]) / np.linalg.norm(ref):8.1e}" for j in range(len(F))))
            sys.exit()
        for spec in sys.argv[2:]:
            a = spec.split(",")
            fam, J, base, hmax = a[0], int(a[1]), int(a[2]), float(a[3])
            dn = int(a[4]) if len(a) > 4 else 8
            colp = int(a[5]) if len(a) > 5 and a[5] else None
            if len(a) > 6:                                              # line quadrature nodes, e.g. 20-10
                pm, ng = (int(v) for v in a[6].split("-"))
                lt.nodes = lambda k, H, pm=pm, ng=ng: (pm, ng)
            else:
                lt.nodes = NODES
            E, t, info = c.run(J, pts=pts, family=fam, **kwargs(J, fam, base, hmax, dn, colp))
            print(f"{fam:8s} J = {J:2d} base {base} hmax {hmax} dn {dn} colp {colp} nodes {a[6] if len(a) > 6 else 'rule'}: unknowns/class {max(info['per_class'])} rows {info['rows']} "
                  f"t = {t:5.1f} s (assembly {info['t_assemble']:.1f}, factor {info['t_factor']:.1f})  error {np.linalg.norm(E - ref) / np.linalg.norm(ref):.2e} "
                  f"optical {info['optical']:.1e} |c| {info['coef_norm']:.1f}", flush=True)
            store[spec] = E
            np.savez("cube_k3_t_tries.npz", **store)
        print("reference:", fref)
