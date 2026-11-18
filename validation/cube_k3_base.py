"""Three wavelengths (k = 6 pi): does the base grid limit the cube layout?
(NOTES.md item 79.) One level of the family with one depth per octave, with
the first depth d0 (which is also the depth of the base grid below the
faces) and the number of base points per side varied. Printed: difference
from Bempp (28800 unknowns) and from the reference of cube_cases.py, optical
mismatch, size of the coefficients.

    python3 cube_k3_base.py J d0:n,n,... d0:n,...       e.g.  7 0.0667:16,20,24 0.2:12,16,20
"""
import os
import sys

import numpy as np

import pec_cube_lightning3d as m

m.K = 6 * np.pi
import cube_mfs as c  # noqa: E402
import line_tree as lt  # noqa: E402

lt.nodes = lambda k, H: (20, 10)

if __name__ == "__main__":
    J = int(sys.argv[1])
    fam = os.environ.get("FAMILY", "balanced")
    dn = int(os.environ.get("DN", "8"))
    hmax = float(os.environ.get("HMAX", "0.065"))
    pts = np.load("bem_cube_h0.05.npz")["pts"]
    bem = np.load("cube_k3_t_bem_n40.npy")
    ref = np.load("cube_cases_k3_ref.npz")["E"]
    store = dict(np.load("cube_k3_t_base.npz")) if os.path.exists("cube_k3_t_base.npz") else {}
    rel = lambda a, b: np.linalg.norm(a - b) / np.linalg.norm(b)
    c.run(2, pts=pts)
    for spec in sys.argv[2:]:
        d0 = float(spec.split(":")[0])
        for nb in (int(v) for v in spec.split(":")[1].split(",")):
            kw = c.dense(J) if fam == "dense" else dict(zip(("kw", "hmax"), c.balanced(J)))
            kw = dict(kw, hmax=hmax)
            kw["kw"] = dict(kw["kw"], n_base=nb, n0=kw["kw"]["n0"] + dn, depth0=d0)
            E, t, info = c.run(J, pts=pts, family=fam, **kw)
            name = f"{fam},J{J},d0={d0:g},nb={nb},dn={dn},h={hmax:g}"
            store[name] = E
            np.savez("cube_k3_t_base.npz", **store)
            print(f"{name}: unknowns/class {max(info['per_class'])} rows {info['rows']} t = {t:5.1f} s  from Bempp {rel(E, bem):.2e}  "
                  f"from the reference {rel(E, ref):.2e}  optical {info['optical']:.1e}  |c| {info['coef_norm']:.1f}", flush=True)
