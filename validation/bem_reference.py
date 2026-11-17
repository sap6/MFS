"""Boundary-element reference solutions with Bempp (EFIE, RWG elements) for
scattering by perfectly conducting bodies, used where no exact solution exists.

Usage:
  python bem_reference.py sphere   # check the Bempp set-up against Mie
  python bem_reference.py cube     # write reference field for the unit cube
"""
import sys
import time

import bem_env  # noqa: F401  (memory for the OpenCL driver; must come before bempp)
import bempp_cl.api as bempp
import numpy as np
import scipy.linalg as sla

import mie

D = np.array([0.0, 0.0, -1.0])
P = np.array([0.0, 1.0, 0.0])


@bempp.callable(complex=True, parameterized=True)
def plane_wave_trace(x, n, domain_index, result, parameters):
    """E_inc x n for E_inc = p exp(i k d.x); parameters = (k, d, p). Defined
    once at module level, so that Numba compiles it once per process (a
    closure defined inside the solver would be recompiled at every call,
    which costs about 3 s each time)."""
    ph = np.exp(1j * parameters[0] * (x[0] * parameters[1] + x[1] * parameters[2] + x[2] * parameters[3]))
    e0, e1, e2 = parameters[4] * ph, parameters[5] * ph, parameters[6] * ph
    result[0] = e1 * n[2] - e2 * n[1]
    result[1] = e2 * n[0] - e0 * n[2]
    result[2] = e0 * n[1] - e1 * n[0]


def efie_solve(grid, k, d, p):
    """EFIE with RWG trial and SNC test functions for E_inc = p exp(i k d.x),
    solved by dense LU in place. Returns the surface current, its space and the
    times (assembly, LU). Bempp assembles with OpenCL when pyopencl and a CPU
    driver are installed and with Numba otherwise (about 20 times slower)."""
    div = bempp.function_space(grid, "RWG", 0)
    curl = bempp.function_space(grid, "SNC", 0)
    t0 = time.time()
    efie = bempp.operators.boundary.maxwell.electric_field(div, div, curl, k)
    A = np.asarray(efie.weak_form().to_dense())
    par = np.concatenate(([k], np.asarray(d, float), np.asarray(p, float)))
    rhs = bempp.GridFunction(div, fun=plane_wave_trace, dual_space=curl, function_parameters=par).projections(curl)
    t1 = time.time()
    lu = sla.lu_factor(A, overwrite_a=True, check_finite=False)
    lam = bempp.GridFunction(div, coefficients=sla.lu_solve(lu, rhs, check_finite=False))
    return lam, div, (t1 - t0, time.time() - t1)


def scattered_field(grid, k, points, times=None):
    """E_sc at points (P,3) for E_inc = P exp(i k D.x), via the EFIE. If a
    list is passed as `times`, (assembly, LU, evaluation) is appended to it."""
    lam, div, (ta, tl) = efie_solve(grid, k, D, P)
    t0 = time.time()
    pot = bempp.operators.potential.maxwell.electric_field(div, points.T, k)
    E = -(pot * lam).T
    if times is not None:
        times.append((ta, tl, time.time() - t0))
    return E, div.global_dof_count


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
