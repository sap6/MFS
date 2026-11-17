"""Bempp reference for the PEC torus (EFIE, RWG), structured triangle mesh.

Usage: python3 bem_torus.py n1,n2,...  [a] [theta phi pol]
n = number of mesh cells around the tube; round(n R / a) cells around the axis.
Writes bem_torus_a{a}_n{n}.npz: field at the test points, far field in the
standard directions, cross-sections, time.
"""
import os
import sys
import time

os.environ["PATH"] = os.path.join(sys.prefix, "bin") + os.pathsep + os.environ.get("PATH", "")

import bem_env  # noqa: F401  (memory for the OpenCL driver; must come before bempp)
import bempp_cl.api as bempp  # noqa: E402

import bem_reference as br  # noqa: E402
import numpy as np  # noqa: E402

from pec_torus import sphere_rule, wave  # noqa: E402

K = 2 * np.pi
RFAR = 1.0e6


def test_points():
    return 2.5 * sphere_rule(10, 20)[0]


def torus_grid(R, a, n):
    """Quasi-uniform structured mesh: n cells around the tube, m around the axis."""
    m = int(round(n * R / a))
    t = 2 * np.pi * np.arange(n) / n
    ph = 2 * np.pi * np.arange(m) / m
    T, P = np.meshgrid(t, ph, indexing="ij")
    P = P + 0.5 * (2 * np.pi / m) * (np.arange(n)[:, None] % 2) * 0     # no stagger: keeps the mesh simple
    rho = R + a * np.cos(T)
    v = np.stack((rho * np.cos(P), rho * np.sin(P), a * np.sin(T)), axis=0).reshape(3, -1)
    idx = lambda i, j: (i % n) * m + (j % m)
    el = []
    for i in range(n):
        for j in range(m):
            p00, p10, p01, p11 = idx(i, j), idx(i + 1, j), idx(i, j + 1), idx(i + 1, j + 1)
            # outward normal: (d/dphi) x (d/dtheta) points outwards
            el.append((p00, p01, p11))
            el.append((p00, p11, p10))
    return bempp.Grid(v, np.array(el, dtype=np.uint32).T)


def solve(grid, d, p, points, dirs):
    lam, div, _ = br.efie_solve(grid, K, d, p)

    def field(x):
        return -(bempp.operators.potential.maxwell.electric_field(div, x.T, K) * lam).T

    F = field(RFAR * dirs) * RFAR * np.exp(-1j * K * RFAR)
    return field(points), F, div.global_dof_count


def main():
    ns = [int(v) for v in sys.argv[1].split(",")]
    a = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
    ang = [float(v) for v in sys.argv[3:6]] if len(sys.argv) > 5 else [0.0, 0.0, 90.0]
    d, p = wave(*ang)
    tag = "" if ang == [0.0, 0.0, 90.0] else "_w" + "_".join(f"{v:g}" for v in ang)
    dirs, w = sphere_rule()
    pts = test_points()
    solve(torus_grid(1.0, a, 6), d, p, pts[:4], dirs[:4])                # warm-up (compilation)
    for n in ns:
        t0 = time.time()
        g = torus_grid(1.0, a, n)
        E, F, ndof = solve(g, d, p, pts, np.vstack((dirs, d[None])))
        t = time.time() - t0
        ext = 4 * np.pi / K * np.imag(np.conj(p) @ F[-1])
        sca = np.sum(w * np.sum(np.abs(F[:-1]) ** 2, axis=1))
        np.savez(f"bem_torus_a{a:g}_n{n}{tag}.npz", E=E, F=F[:-1], ndof=ndof, t=t, ext=ext, sca=sca)
        print(f"Bempp torus a={a} n={n} dofs={ndof} time={t:.1f}s sigma_ext={ext:.6f} sigma_sca={sca:.6f}", flush=True)


if __name__ == "__main__":
    main()
