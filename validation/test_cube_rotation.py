"""Check of the rotation about the body diagonal in cube_mfs.run: a general
wave solved with four class matrices (run) against the direct solve with all
eight (Cube), and the standard wave against the stored field.

    python3 test_cube_rotation.py [J]
"""
import sys
import time

import numpy as np

import cube_mfs as c
import pec_cube_lightning3d as m

J = int(sys.argv[1]) if len(sys.argv) > 1 else 4
pts = np.load("bem_cube_h0.05.npz")["pts"]
d, p = c.wave(60.0, 25.0, 35.0)
c.run(2, d, p, pts, family="dense")
E, t, info = c.run(J, d, p, pts, family="dense")
t0 = time.time()
cube = c.Cube(2 * J, **c.dense(J))
coef, res = cube.solve(d, p)
E8 = cube.field(coef, pts)
t8 = time.time() - t0
print(f"general wave, J = {J}: four matrices {t:.1f} s, eight {t8:.1f} s, fields differ by {np.linalg.norm(E - E8) / np.linalg.norm(E8):.1e}, "
      f"residual {info['residual']:.3e} against {res:.3e}, optical {info['optical']:.1e}")
E, t, info = c.run(J, pts=pts, family="dense")
ref = np.load("fig_cost_cube.npz")[f"dense_E{J}"]
print(f"standard wave, J = {J}: {t:.1f} s, differs from the stored field by {np.linalg.norm(E - ref) / np.linalg.norm(ref):.1e}")
