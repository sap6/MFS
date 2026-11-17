"""Is Bempp run at its best?  Sphere against Mie with each assembly backend
(Numba, OpenCL) and precision (double, single), time split into assembly,
solve and field evaluation.

    python3 bem_backend.py 0.14,0.1            # mesh sizes
    python3 bem_backend.py 0.1 opencl:single   # only one configuration

Needs pyopencl and an OpenCL CPU driver (pip install pyopencl
pocl-binary-distribution) for the OpenCL rows. Run alone on the machine.
"""
import sys
import time

import bempp_cl.api as bempp
import numpy as np
from scipy.linalg import lu_factor, lu_solve

import mie
from pec_sphere_convergence import fibonacci_sphere

D = np.array([0.0, 0.0, -1.0])
P = np.array([0.0, 1.0, 0.0])
K = 2 * np.pi


def solve(grid, k, pts, device, precision, d=D, p=P):
    """E_sc at pts; returns field, dofs and the three times."""
    div = bempp.function_space(grid, "RWG", 0)
    curl = bempp.function_space(grid, "SNC", 0)
    t0 = time.time()
    efie = bempp.operators.boundary.maxwell.electric_field(div, div, curl, k, device_interface=device,
                                                           precision=precision)
    A = efie.weak_form().to_dense()
    t1 = time.time()

    @bempp.complex_callable
    def trace(x, n, domain_index, result):
        inc = p * np.exp(1j * k * (x[0] * d[0] + x[1] * d[1] + x[2] * d[2]))
        result[:] = np.cross(inc, n)

    rhs = bempp.GridFunction(div, fun=trace, dual_space=curl)
    lam = lu_solve(lu_factor(A, overwrite_a=True, check_finite=False), rhs.projections(curl).astype(A.dtype))
    t2 = time.time()
    pot = bempp.operators.potential.maxwell.electric_field(div, pts.T, k, device_interface=device,
                                                           precision=precision)
    E = -(pot * bempp.GridFunction(div, coefficients=lam)).T
    t3 = time.time()
    return E, div.global_dof_count, (t1 - t0, t2 - t1, t3 - t2)


if __name__ == "__main__":
    hs = [float(s) for s in sys.argv[1].split(",")]
    configs = [("numba", "double"), ("numba", "single"), ("opencl", "double"), ("opencl", "single")]
    if len(sys.argv) > 2:
        configs = [tuple(c.split(":")) for c in sys.argv[2].split(",")]
    pts = 2.0 * fibonacci_sphere(400)
    exact = mie.scattered_field(pts, K, D, P)
    for device, precision in configs:
        solve(bempp.shapes.sphere(h=0.4), K, pts, device, precision)          # warm-up (compilation)
        for h in hs:
            grid = bempp.shapes.sphere(h=h)
            E, ndof, (ta, ts, te) = solve(grid, K, pts, device, precision)
            err = np.linalg.norm(E - exact) / np.linalg.norm(exact)
            print(f"{device:6s} {precision:6s} h={h:<5} dofs={ndof:6d} assembly={ta:7.1f}s solve={ts:6.1f}s "
                  f"eval={te:5.1f}s total={ta + ts + te:7.1f}s error={err:.3e}", flush=True)
