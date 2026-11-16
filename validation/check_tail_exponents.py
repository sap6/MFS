"""Cube: are the tail exponents the derived ones? (DERIVATIONS.md, sections B and D)

Solves the standard case at one depth J of the balanced family with different
exponent sets for the edge tail, and with and without the vertex tail, and
prints the error of each against the stored J = 9 field (fig_cost_cube.npz)
and the optical-theorem mismatch.

    python3 check_tail_exponents.py J        # about 6 solves; J = 5 takes a few minutes

Result of 16 Nov 2026 (J = 5 / J = 7, no vertex tail): the derived leading
exponents (2/3 | 5/3 | 4/3 | 7/3) give 0.314% / 0.049%, the earlier set with a
second, integer exponent per component 0.327% / 0.049%, the derived set with
its second terms 0.316% at J = 5. The vertex tail takes 0.31% to 0.20% (J = 5)
and 0.049% to 0.029% (J = 7).
"""
import sys
import time

import numpy as np

import pec_cube_lightning3d as m
from fig_cost_cube import COL, balanced

VARIANTS = (
    ("derived, leading term only (2/3 | 5/3 | 4/3 | 7/3), no vertex tail", dict(vtail_n=0)),
    ("the same with the vertex tail (the default)", dict()),
    ("vertex tail with two exponents per kind", dict(vtail_exp=2)),
    ("derived, two terms (2/3, 8/3 | 5/3, 11/3 | 4/3, 8/3 | 7/3, 11/3), no vertex tail",
     dict(vtail_n=0, tail_e=(2 / 3, 8 / 3), tail_er=(5 / 3, 11 / 3), tail_m=(4 / 3, 8 / 3), tail_mr=(7 / 3, 11 / 3))),
    ("earlier set with integer second exponents (2/3, 1 | 5/3, 2 | 1, 4/3 | 2, 7/3), no vertex tail",
     dict(vtail_n=0, tail_e=(2 / 3, 1.0), tail_er=(5 / 3, 2.0), tail_m=(1.0, 4 / 3), tail_mr=(2.0, 7 / 3))),
    ("one exponent for all axis (2/3) and all radial (5/3) components, no vertex tail",
     dict(vtail_n=0, tail_e=(2 / 3,), tail_er=(5 / 3,), tail_m=(2 / 3,), tail_mr=(5 / 3,))),
)

if __name__ == "__main__":
    J = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    pts = np.load("bem_cube_h0.05.npz")["pts"]
    ref = np.load("fig_cost_cube.npz")["tail_E9"]
    for name, extra in VARIANTS:
        kw, hmax = balanced(J)
        kw.update(extra)
        t0 = time.time()
        src = m.SymSources(J, tail_n=9, **kw)
        lv = max(2, int(np.ceil(np.log(src.depths[-1] / (2 * m.H)) / np.log(COL["col_q"]))) + 1)
        c, info = m.solve(src, lv, COL["col_p"], COL["col_q"], COL["beta"], hmax=hmax)
        E = (src.matrix(pts) @ c).reshape(-1, 3)
        ext, sca = m.cross_sections(src, c)
        print(f"J={J} {name}:\n    unknowns {info['unknowns']}, error against J = 9: {np.linalg.norm(E - ref) / np.linalg.norm(ref):.2e}, "
              f"optical mismatch {abs(ext - sca) / sca:.1e}, coefficient norm {info['coef_norm']:.1e} ({time.time() - t0:.0f}s)",
              flush=True)
