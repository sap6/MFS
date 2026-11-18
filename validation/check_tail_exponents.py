"""Cube: are the tail exponents the derived ones? (DERIVATIONS.md, sections B and D)

Solves the standard case at one level J of the layout of the paper (two
depths per octave) with different exponent sets for the edge tail, and with
and without the vertex tail, and prints the error of each against the
reference field (fig_cost_cube.npz) and the optical-theorem mismatch.

    python3 check_tail_exponents.py [J]      # default 6; six solves, two minutes

The derived leading exponents are 2/3 (electric, along the edge), 5/3
(electric, along the bisector), 4/3 and 7/3 (magnetic).
"""
import sys

import numpy as np

import cube_mfs as c

VARIANTS = (
    ("derived exponents, vertex tail (the layout of the paper)", dict()),
    ("the same without the vertex tail", dict(vtail_n=0)),
    ("vertex tail with two exponents per kind", dict(vtail_exp=2)),
    ("edge tail with the next term of each expansion (2/3, 8/3 | 5/3, 11/3 | 4/3, 8/3 | 7/3, 11/3)",
     dict(tail_e=(2 / 3, 8 / 3), tail_er=(5 / 3, 11 / 3), tail_m=(4 / 3, 8 / 3), tail_mr=(7 / 3, 11 / 3))),
    ("exponents of the two components along the bisector lowered by one (2/3 | 2/3 | 4/3 | 4/3)",
     dict(tail_er=(2 / 3,), tail_mr=(4 / 3,))),
    ("one exponent for all components along the edge (2/3) and along the bisector (5/3)",
     dict(tail_m=(2 / 3,), tail_mr=(5 / 3,))),
)

if __name__ == "__main__":
    J = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    pts = np.load("bem_cube_h0.05.npz")["pts"]
    ref = np.load("fig_cost_cube.npz")["ref_E"]
    c.run(2, pts=pts, family="dense")                      # compilation
    for name, extra in VARIANTS:
        E, t, info = c.run(J, pts=pts, family="dense", kw=dict(c.dense(J)["kw"], **extra))
        print(f"J={J} {name}:\n    unknowns per class {max(info['per_class'])}, error {np.linalg.norm(E - ref) / np.linalg.norm(ref):.3e}, "
              f"optical mismatch {info['optical']:.1e}, coefficient norm {info['coef_norm']:.1e} ({t:.0f} s)", flush=True)
