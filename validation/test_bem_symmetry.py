"""Checks of the boundary element solves with symmetry (bem_mirror.py)
against full solves on small meshes: the cube (three mirror planes and the
rotation about the body diagonal, standard and general wave), the torus
(rotation of the mesh about the axis) and the almond (two mirror planes).
The mirror solves agree with the full ones to the quadrature error of the
library (about 3e-5 with the default orders), the rotation to rounding.

    python3 test_bem_symmetry.py
"""
import bem_env  # noqa: F401
import bempp_cl.api as bempp
import numpy as np

import bem_cube_sym as bc
import bem_reference as br
import bem_torus as bt
import pec_almond_skeleton as a
from cube_mfs import wave as cube_wave
from pec_sphere_convergence import fibonacci_sphere
from pec_torus import wave as torus_wave

rel = lambda x, y: np.linalg.norm(x - y) / np.linalg.norm(y)

if __name__ == "__main__":
    pts = 2.0 * fibonacci_sphere(100)
    worst = 0.0
    for name, (d, p) in (("standard wave", (bc.D, bc.P)), ("general wave", cube_wave(60.0, 25.0, 35.0))):
        E, N, sizes, ts, grid = bc.solve(8, 2.0, pts, d, p)
        lam, div, _ = br.efie_solve(grid, bc.K, d, p)
        Ef = -(bempp.operators.potential.maxwell.electric_field(div, pts.T, bc.K) * lam).T
        e = rel(E, Ef)
        worst = max(worst, e)
        print(f"cube, {name}, {N} unknowns: symmetric against full solve {e:.1e}")
    for ang in ((0.0, 0.0, 90.0), (60.0, 25.0, 35.0)):
        d, p = torus_wave(*ang)
        E1 = bt.solve(bt.torus_grid(1.0, 0.5, 12), d, p, pts)[0]
        E2 = bt.solve_rot(1.0, 0.5, 12, d, p, pts)[0]
        E3 = bt.solve_sym(1.0, 0.5, 12, d, p, pts)[0]
        print(f"torus, wave {ang}: rotation against full solve {rel(E2, E1):.1e}; three mirrors (another mesh) against it {rel(E3, E1):.1e}")
        assert rel(E2, E1) < 1e-10
    k = 2 * np.pi * 1.19e9 / 299792458.0
    x = 0.3 * fibonacci_sphere(100)
    (n1, Es, _), = a.bem_reference(k, x, ((20, 20),), sym=True)
    (n2, Ef, _), = a.bem_reference(k, x, ((20, 20),), sym=False)
    print(f"almond, {n1} unknowns: two mirrors against the full solve on the unsymmetric mesh {rel(Es, Ef):.1e} (different triangulations)")
    assert worst < 2e-4
    print("ok")
