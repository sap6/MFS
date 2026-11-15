"""Boundary-element reference solutions with Bempp (EFIE, RWG elements) for
scattering by perfectly conducting bodies, used where no exact solution exists.

Usage:
  python bem_reference.py sphere   # check the Bempp set-up against Mie
  python bem_reference.py cube     # write reference field for the unit cube
"""
import sys
import time

import bempp_cl.api as bempp
import numpy as np

import mie

D = np.array([0.0, 0.0, -1.0])
P = np.array([0.0, 1.0, 0.0])


def scattered_field(grid, k, points):
    """E_sc at points (P,3) for E_inc = P exp(i k D.x), via the EFIE."""
    div = bempp.function_space(grid, "RWG", 0)
    curl = bempp.function_space(grid, "SNC", 0)
    efie = bempp.operators.boundary.maxwell.electric_field(div, div, curl, k)

    @bempp.complex_callable
    def tangential_trace(x, n, domain_index, result):
        inc = P * np.exp(1j * k * (x[0] * D[0] + x[1] * D[1] + x[2] * D[2]))
        result[:] = np.cross(inc, n)

    rhs = bempp.GridFunction(div, fun=tangential_trace, dual_space=curl)
    lam = bempp.linalg.lu(efie, rhs)
    pot = bempp.operators.potential.maxwell.electric_field(div, points.T, k)
    return -(pot * lam).T, div.global_dof_count


def test_points(n=400, r=2.0):
    from pec_sphere_convergence import fibonacci_sphere
    return r * fibonacci_sphere(n)


def rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


if __name__ == "__main__":
    which = sys.argv[1]
    k = 2 * np.pi
    pts = test_points()
    if which == "sphere":
        exact = mie.scattered_field(pts, k, D, P)
        for h in (0.2, 0.1, 0.07):
            t0 = time.time()
            E, ndof = scattered_field(bempp.shapes.sphere(h=h), k, pts)
            print(f"sphere h={h}: dofs={ndof} rel err vs Mie on r=2: {rel(E, exact):.2e} "
                  f"({time.time() - t0:.0f}s)", flush=True)
    else:
        prev = None
        for h in (0.1, 0.07, 0.05):
            t0 = time.time()
            grid = bempp.shapes.cube(length=1, origin=(-0.5, -0.5, -0.5), h=h)
            E, ndof = scattered_field(grid, k, pts)
            msg = f"cube h={h}: dofs={ndof} ({time.time() - t0:.0f}s)"
            if prev is not None:
                msg += f" change from previous mesh: {rel(prev, E):.2e}"
            print(msg, flush=True)
            np.savez(f"bem_cube_h{h}.npz", pts=pts, E=E, k=k, h=h, ndof=ndof)
            prev = E
