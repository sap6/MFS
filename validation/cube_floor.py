"""What sets the floor of the cube family near 5e-5 (NOTES item 42)?
Variants at a fixed deepest depth, compared with the stored J = 9 field.

    python3 cube_floor.py NAME [NAME ...]      # see VARIANTS
"""
import sys
import time

import numpy as np

import cube_mfs as cm
import pec_cube_lightning3d as m

S = 2 ** -0.5
VARIANTS = {   # name: (J, source keywords overriding balanced(9), hmax, collocation)
    "J9": (9, {}, 0.1, {}),
    "J10": (10, {}, 0.1, {}),
    "dense_p4": (18, dict(q=S, n0=14), 0.1, dict(col_p=4)),
    "dense_colq": (18, dict(q=S, n0=14), 0.1, dict(col_q=S)),
    "dense_tail": (18, dict(q=S, n0=14, tail_q=S), 0.1, dict(col_q=S)),
    "J9_colq": (9, {}, 0.1, dict(col_q=S)),
    "J9_vrays": (9, dict(vrays=4), 0.1, {}),
    "J9_vexp2": (9, dict(vtail_exp=2), 0.1, {}),
    "J9_beta03": (9, {}, 0.1, dict(beta=0.3)),
    "J9_beta07": (9, {}, 0.1, dict(beta=0.7)),
    "J9_modes": (9, dict(n0=20), 0.1, {}),
    "J9_gamma3": (9, dict(gamma=3.0), 0.1, {}),
    "J9_base12": (9, dict(n_base=12), 0.07, {}),
}


def rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


if __name__ == "__main__":
    d = np.load("fig_cost_cube.npz")
    pts, ref = d["pts"], d["tail_E9"]
    cm.run(2, pts=pts)
    fields = {}
    for name in sys.argv[1:]:
        J, over, hmax, colo = VARIANTS[name]
        kw, _ = cm.balanced(9)
        tail = 18 if over.get("tail_q") else 9
        kw.update(over)
        if "tail_q" in over:
            kw["vtail_n"] = 18
        col = dict(cm.COL, **colo)
        t0 = time.time()
        cube = cm.Cube(J, kw=kw, hmax=hmax, tail=tail, classes=cm.classes_of(m.D, m.P), col=col)
        coef, res = cube.solve(m.D, m.P)
        E = cube.field(coef, pts)
        t = time.time() - t0
        ext, sca = cube.cross_sections(coef, m.D, m.P)
        fields[name] = E
        np.save(f"cube_floor_t_{name}.npy", E)
        others = " ".join(f"vs {k}: {rel(E, v):.1e}" for k, v in fields.items() if k != name)
        print(f"{name:11s} unknowns/class={max(cube.unknowns.values())} rows={cube.rows} t={t:.0f}s ({type(next(iter(cube.fac.values()))).__name__}) "
              f"residual={res:.2e} optical={abs(ext - sca) / sca:.1e} sigma_ext={ext:.6f} |c|={cube.coef_norm:.1f} vs stored J9: {rel(E, ref):.2e} {others}", flush=True)
