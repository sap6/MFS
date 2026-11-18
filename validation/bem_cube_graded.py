"""Bempp (EFIE, RWG, dense LU) for the unit cube on meshes graded towards the
edges: an independent reference that does not suffer from the h^(4/3) rate of
quasi-uniform meshes, and a fairer ("tuned") boundary-element competitor.

Each face carries a tensor-product grid with n intervals per side whose nodes
are t_i = +-0.5 (1 - (1 - |i|/(n/2))^g), i.e. clustered towards the edges with
grading exponent g (g = 1: uniform); every rectangle is cut into two
triangles. The elements near the edges are anisotropic (long along the edge).

    python3 bem_cube_graded.py n1,n2,... g [regular_order] [singular_order] [--wave=theta,phi,psi]

--wave gives an oblique plane wave in the convention of cube_mfs.wave
(degrees); the files then end in _obl. Without it the wave travels along -z
with E along y, as everywhere else.
Writes bem_cube_graded_n{n}_g{g}.npz with the field at the 400 test points,
the far field, the cross-sections and the assembly and factorisation times
separately. Compilation is excluded by a small warm-up solve.
"""
import os
import sys
import time

os.environ["PATH"] = os.path.join(sys.prefix, "bin") + os.pathsep + os.environ.get("PATH", "")

import bem_env  # noqa: F401  (memory for the OpenCL driver; must come before bempp)
import bempp_cl.api as bempp  # noqa: E402

import bem_reference as br  # noqa: E402
import numpy as np  # noqa: E402
import scipy.linalg as sla  # noqa: E402

from bem_optical import D, P, K, RFAR, sphere_rule  # noqa: E402


def graded_nodes(n, g):
    s = np.linspace(-1.0, 1.0, n + 1)
    return 0.5 * np.sign(s) * (1.0 - (1.0 - np.abs(s)) ** g)


def cube_grid(n, g):
    """Tensor graded triangulation of the surface of [-1/2, 1/2]^3, outward normals."""
    t = graded_nodes(n, g)
    index = {}
    verts = []

    def vid(p):
        key = tuple(np.round(p, 12))
        if key not in index:
            index[key] = len(verts)
            verts.append(p)
        return index[key]

    tris = []
    for axis in range(3):
        for side in (-0.5, 0.5):
            a, b = (axis + 1) % 3, (axis + 2) % 3     # (a, b, axis) is right-handed
            for i in range(n):
                for j in range(n):
                    q = []
                    for (ii, jj) in ((i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)):
                        p = np.zeros(3)
                        p[axis], p[a], p[b] = side, t[ii], t[jj]
                        q.append(vid(p))
                    # diagonal chosen towards the nearest face corner so the pattern is mirror symmetric
                    if (i < n // 2) == (j < n // 2):
                        cells = [(q[0], q[1], q[2]), (q[0], q[2], q[3])]
                    else:
                        cells = [(q[0], q[1], q[3]), (q[1], q[2], q[3])]
                    for c in cells:
                        tris.append(c if side > 0 else c[::-1])
    return bempp.Grid(np.array(verts).T, np.array(tris, dtype=np.uint32).T)


def optical(far_fun, d, p):
    """Extinction and scattering cross-sections from the far-field function."""
    dirs, w = sphere_rule()
    F = far_fun(np.vstack((dirs, d[None, :])))
    return 4 * np.pi / K * np.imag(np.conj(p) @ F[-1]), np.sum(w * np.sum(np.abs(F[:-1]) ** 2, axis=1)), F[:-1]


def solve(grid, points, D=D, P=P):
    lam, div, (ta, tl) = br.efie_solve(grid, K, D, P)
    t2 = time.time()

    def field(x):
        return -(bempp.operators.potential.maxwell.electric_field(div, x.T, K) * lam).T

    def far(dirs):
        return field(RFAR * dirs) * RFAR * np.exp(-1j * K * RFAR)

    E = field(points)
    t3 = time.time()
    return E, far, div.global_dof_count, (ta, tl, t3 - t2)


if __name__ == "__main__":
    wave = [a[7:] for a in sys.argv if a.startswith("--wave=")]
    args = [a for a in sys.argv if not a.startswith("--")]
    ns = [int(a) for a in args[1].split(",")]
    g = float(args[2])
    reg = int(args[3]) if len(args) > 3 else 4
    sing = int(args[4]) if len(args) > 4 else 6
    if wave:
        from cube_mfs import wave as wave_dp
        d, p = wave_dp(*(float(v) for v in wave[0].split(",")))
    else:
        d, p = D, P
    bempp.GLOBAL_PARAMETERS.quadrature.regular = reg
    bempp.GLOBAL_PARAMETERS.quadrature.singular = sing
    pts = np.load("bem_cube_h0.05.npz")["pts"]
    ref = np.load("fig_cost_cube.npz")["ref_E"] if not wave else None       # the MFS reference field
    solve(cube_grid(4, 1.0), pts)                       # warm-up: compilation
    for n in ns:
        E, far, ndof, (ta, tl, te) = solve(cube_grid(n, g), pts, d, p)
        ext, sca, F = optical(far, d, p)
        err = np.linalg.norm(E - ref) / np.linalg.norm(ref) if ref is not None else np.nan
        tag = f"n{n}_g{g:g}" + ("" if (reg, sing) == (4, 6) else f"_q{reg}_{sing}") + ("_obl" if wave else "")
        print(f"cube graded n={n} g={g:g} quad={reg}/{sing} dofs={ndof} assembly={ta:.0f}s lu={tl:.0f}s eval={te:.0f}s "
              f"sigma_ext={ext:.6f} sigma_sca={sca:.6f} diff from the MFS reference: {err:.2e}", flush=True)
        np.savez(f"bem_cube_graded_{tag}.npz", pts=pts, E=E, F=F, k=K, n=n, g=g, ndof=ndof, ext=ext, sca=sca, d=d, p=p,
                 t_assembly=ta, t_lu=tl, t_eval=te, t_solve=ta + tl + te)
