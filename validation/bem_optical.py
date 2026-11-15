"""Independent accuracy check of the Bempp cube reference: the optical theorem.

For a perfect conductor the extinction cross-section (from the forward
scattering amplitude) equals the scattering cross-section (integral of the
far-field intensity). The relative mismatch is an error measure that needs no
second solver. Also varies the regular quadrature order, to see whether the
reference is limited by the mesh or by the quadrature.

Usage: python3 bem_optical.py h1,h2,... [regular_order] [singular_order]
Writes bem_cube_h{h}_q{order}.npz (field at the 400 test points and far field).
"""
import os
import sys
import time

os.environ["PATH"] = os.path.join(sys.prefix, "bin") + os.pathsep + os.environ.get("PATH", "")

import bempp_cl.api as bempp  # noqa: E402
import numpy as np  # noqa: E402

D = np.array([0.0, 0.0, -1.0])
P = np.array([0.0, 1.0, 0.0])
K = 2 * np.pi
RFAR = 1.0e6


def sphere_rule(nt=40, nphi=80):
    """Gauss-Legendre in cos(theta) times trapezoid in phi."""
    g, w = np.polynomial.legendre.leggauss(nt)
    ph = (np.arange(nphi) + 0.5) / nphi * 2 * np.pi
    ct, phh = np.meshgrid(g, ph, indexing="ij")
    st = np.sqrt(1 - ct ** 2)
    dirs = np.column_stack(((st * np.cos(phh)).ravel(), (st * np.sin(phh)).ravel(), ct.ravel()))
    return dirs, np.repeat(w, nphi) * 2 * np.pi / nphi


def optical(far_fun):
    """far_fun(dirs) -> F with E_sc ~ F exp(ikr)/r. Returns sigma_ext, sigma_sca."""
    dirs, w = sphere_rule()
    F = far_fun(np.vstack((dirs, D[None, :])))
    ext = 4 * np.pi / K * np.imag(np.conj(P) @ F[-1])
    sca = np.sum(w * np.sum(np.abs(F[:-1]) ** 2, axis=1))
    return ext, sca, F[:-1]


def solve(h, points):
    grid = bempp.shapes.cube(length=1, origin=(-0.5, -0.5, -0.5), h=h)
    div = bempp.function_space(grid, "RWG", 0)
    curl = bempp.function_space(grid, "SNC", 0)
    efie = bempp.operators.boundary.maxwell.electric_field(div, div, curl, K)

    @bempp.complex_callable
    def trace(x, n, domain_index, result):
        inc = P * np.exp(1j * K * (x[0] * D[0] + x[1] * D[1] + x[2] * D[2]))
        result[:] = np.cross(inc, n)

    lam = bempp.linalg.lu(efie, bempp.GridFunction(div, fun=trace, dual_space=curl))

    def field(x):
        return -(bempp.operators.potential.maxwell.electric_field(div, x.T, K) * lam).T

    def far(dirs):
        return field(RFAR * dirs) * RFAR * np.exp(-1j * K * RFAR)

    return field(points), far, div.global_dof_count


if __name__ == "__main__":
    hs = [float(a) for a in sys.argv[1].split(",")]
    reg = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    sing = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    bempp.GLOBAL_PARAMETERS.quadrature.regular = reg
    bempp.GLOBAL_PARAMETERS.quadrature.singular = sing
    pts = np.load("bem_cube_h0.05.npz")["pts"]
    for h in hs:
        t0 = time.time()
        E, far, ndof = solve(h, pts)
        t1 = time.time()
        ext, sca, F = optical(far)
        print(f"cube h={h} quad={reg}/{sing} dofs={ndof} solve={t1 - t0:.0f}s  sigma_ext={ext:.6f} "
              f"sigma_sca={sca:.6f} mismatch={abs(ext - sca) / sca:.2e}", flush=True)
        np.savez(f"bem_cube_h{h}_q{reg}.npz", pts=pts, E=E, F=F, k=K, h=h, ndof=ndof, ext=ext, sca=sca,
                 t_solve=t1 - t0)
