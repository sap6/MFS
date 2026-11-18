"""MFS for exterior scattering by a perfectly conducting NASA almond,
checked against a Bempp (EFIE/RWG) reference.

Source placement. Every cross-section of the almond is a 3:1 ellipse, whose
scattered field continues analytically inside only up to the focal set, which
lies within 6% of the rim. A uniformly shrunk almond therefore cannot serve as
the source surface. Sources are placed instead on a *confocal* surface:
cross-section semi-axes (sqrt(a^2 - (1 - tau^2) b^2), tau b), 0 < tau < 1.

Usage: python pec_almond.py [freq_GHz]
"""
import json
import sys
import time

import numpy as np
from scipy.linalg import lstsq

import mie
from pec_sphere_convergence import dyadic, evaluate, fibonacci_sphere, tangent_basis

L = 0.2524
T0, T1 = -0.416667, 0.58333
D = np.array([0.0, 0.0, -1.0])
P = np.array([0.0, 1.0, 0.0])


def semi_axes(t):
    """Cross-section semi-axes (a, b) in metres at parameter t."""
    t = np.asarray(t, float)
    f1 = np.sqrt(np.clip(1 - (t / 0.416667) ** 2, 0, None))
    f2 = np.sqrt(np.clip(1 - (t / 2.08335) ** 2, 0, None)) - 0.96
    a = np.where(t < 0, 0.193333 * f1, 4.83345 * f2) * L
    b = np.where(t < 0, 0.064444 * f1, 1.61115 * f2) * L
    return np.clip(a, 0, None), np.clip(b, 0, None)


def surface(t, s, tau=None):
    """Points on the almond (tau=None) or on the confocal source surface."""
    a, b = semi_axes(t)
    x = L * t
    if tau is not None:
        lam = (1 - tau ** 2) * b ** 2
        if np.any(t < 0):                       # rear half: confocal ellipsoid, x-axis shrinks too
            A = 0.416667 * L
            b0 = 0.064444 * L
            A1 = np.sqrt(A ** 2 - (1 - tau ** 2) * b0 ** 2)
            x = np.where(t < 0, x * A1 / A, x)
        a, b = np.sqrt(np.clip(a ** 2 - lam, 0, None)), tau * b
    return np.stack((x * np.ones_like(s), a * np.cos(s), b * np.sin(s)), axis=-1)


def t_nodes(h, end=0.004):
    """Parameter values giving rings uniformly spaced by h (metres) along x,
    stopping a distance `end` (in t) short of the rounded end and of the tip."""
    n = int(np.ceil((T1 - T0 - 2 * end) * L / h))
    return T0 + end + (np.arange(n) + 0.5) / n * (T1 - T0 - 2 * end)


def point_cloud(h, n_s, tau=None, b_min=0.0):
    """Points and (for the surface) quadrature-like weights and normals.
    Rings with minor semi-axis below b_min * max are skipped (used for sources)."""
    t = t_nodes(h)
    t = t[semi_axes(t)[1] >= b_min * 0.064444 * L]
    pts, wts, nrm = [], [], []
    dt = np.gradient(t)
    for ti, dti in zip(t, dt):
        a, b = semi_axes(ti)
        m = max(6, int(round(n_s * np.sqrt(a / (0.193333 * L)))))   # fewer points on small sections
        s = (np.arange(m) + 0.5 * (int(round(ti * 1e4)) % 2)) / m * 2 * np.pi
        p = surface(np.full(m, ti), s, tau)
        pts.append(p)
        if tau is None:
            e = 1e-6
            ps = (surface(np.full(m, ti), s + e) - surface(np.full(m, ti), s - e)) / (2 * e)
            pt = (surface(np.full(m, ti + e), s) - surface(np.full(m, ti - e), s)) / (2 * e)
            n = np.cross(ps, pt)
            jac = np.linalg.norm(n, axis=1)
            wts.append(jac * dti * 2 * np.pi / m)
            nrm.append(n / jac[:, None])
    pts = np.concatenate(pts)
    if tau is not None:
        return pts
    nrm = np.concatenate(nrm)
    nrm *= np.sign(np.sum(nrm * (pts - np.array([0.02, 0, 0])), axis=1))[:, None]   # outward
    return pts, np.concatenate(wts), nrm


def mesh(n_t, n_s):
    """Structured triangulation of the almond for Bempp (outward normals)."""
    th = (np.arange(1, n_t) / n_t) * np.pi
    t = 0.5 * (T0 + T1) - 0.5 * (T1 - T0) * np.cos(th)
    s = np.arange(n_s) / n_s * 2 * np.pi
    rings = [surface(np.full(n_s, ti), s) for ti in t]
    verts = np.concatenate([np.array([[L * T0, 0, 0]])] + rings + [np.array([[L * T1, 0, 0]])])
    idx = lambda i, j: 1 + i * n_s + (j % n_s)
    last = len(verts) - 1
    tris = []
    for j in range(n_s):
        tris.append((0, idx(0, j + 1), idx(0, j)))
        tris.append((last, idx(n_t - 2, j), idx(n_t - 2, j + 1)))
        for i in range(n_t - 2):
            tris.append((idx(i, j), idx(i, j + 1), idx(i + 1, j)))
            tris.append((idx(i, j + 1), idx(i + 1, j + 1), idx(i + 1, j)))
    tris = np.array(tris)
    # make orientation outward
    c = verts[tris].mean(axis=1)
    n = np.cross(verts[tris[:, 1]] - verts[tris[:, 0]], verts[tris[:, 2]] - verts[tris[:, 0]])
    flip = np.sum(n * (c - np.array([0.02, 0, 0])), axis=1) < 0
    tris[flip] = tris[flip][:, [0, 2, 1]]
    return verts, tris


def mesh_quarter(n_t, n_s):
    """The part of mesh(n_t, n_s) with y, z >= 0 (n_s a multiple of 4), for
    the boundary element solve with the two mirror planes (bem_mirror.py)."""
    assert n_s % 4 == 0
    q = n_s // 4
    th = (np.arange(1, n_t) / n_t) * np.pi
    t = 0.5 * (T0 + T1) - 0.5 * (T1 - T0) * np.cos(th)
    s = np.arange(q + 1) / n_s * 2 * np.pi
    rings = [surface(np.full(q + 1, ti), s) for ti in t]
    verts = np.concatenate([np.array([[L * T0, 0, 0]])] + rings + [np.array([[L * T1, 0, 0]])])
    verts[np.abs(verts) < 1e-15] = 0.0
    idx = lambda i, j: 1 + i * (q + 1) + j
    last = len(verts) - 1
    tris = []
    for j in range(q):
        tris.append((0, idx(0, j + 1), idx(0, j)))
        tris.append((last, idx(n_t - 2, j), idx(n_t - 2, j + 1)))
        for i in range(n_t - 2):
            tris.append((idx(i, j), idx(i, j + 1), idx(i + 1, j)))
            tris.append((idx(i, j + 1), idx(i + 1, j + 1), idx(i + 1, j)))
    tris = np.array(tris)
    c = verts[tris].mean(axis=1)
    n = np.cross(verts[tris[:, 1]] - verts[tris[:, 0]], verts[tris[:, 2]] - verts[tris[:, 0]])
    flip = np.sum(n * (c - np.array([0.02, 0, 0])), axis=1) < 0
    tris[flip] = tris[flip][:, [0, 2, 1]]
    return verts, tris


def mfs(k, h, n_s, tau, src_frac=0.65, b_min=0.35, inc=None, cond=1e-13):
    inc = inc or (lambda p: mie.plane_wave(p, k, D, P))
    x, w, nrm = point_cloud(h, n_s)
    y = point_cloud(2 * h, max(6, int(n_s * src_frac)), tau, b_min)
    T = tangent_basis(nrm)
    sw = np.sqrt(w / w.sum())
    A = np.einsum("mab,mnbc->manc", T, dyadic(x, y, k)) * sw[:, None, None, None]
    A = A.reshape(2 * len(x), 3 * len(y))
    b = (-np.einsum("mab,mb->ma", T, inc(x)) * sw[:, None]).reshape(-1)
    scale = np.linalg.norm(A, axis=0)
    t0 = time.time()
    c, _, rank, _ = lstsq(A / scale, b, lapack_driver="gelsy", cond=cond)
    alpha = (c / scale).reshape(len(y), 3)
    info = {"M": len(x), "N": len(y), "rank": int(rank), "t_solve": time.time() - t0,
            "coef_norm": float(np.linalg.norm(c)), "tau": tau}
    return y, alpha, info


def boundary_residual(k, y, alpha, inc=None, h=0.0013, n_s=70):
    inc = inc or (lambda p: mie.plane_wave(p, k, D, P))
    x, w, nrm = point_cloud(h, n_s)
    T = tangent_basis(nrm)
    rhs = -np.einsum("mab,mb->ma", T, inc(x))
    res = np.einsum("mab,mb->ma", T, evaluate(x, y, alpha, k)) - rhs
    return float(np.sqrt(np.sum(w * np.sum(np.abs(res) ** 2, 1)) / np.sum(w * np.sum(np.abs(rhs) ** 2, 1))))


def bem(k, n_t, n_s, pts):
    import bem_env  # noqa: F401  (memory for the OpenCL driver; must come before bempp)
    import bempp_cl.api as bempp
    from bem_reference import scattered_field
    v, e = mesh(n_t, n_s)
    grid = bempp.Grid(v.T, e.T.astype("uint32"))
    return scattered_field(grid, k, pts)


if __name__ == "__main__":
    freq = float(sys.argv[1]) if len(sys.argv) > 1 else 1.19
    k = 2 * np.pi * freq * 1e9 / 299792458.0
    pts = 0.3 * fibonacci_sphere(400)
    print(f"f = {freq} GHz, k = {k:.3f}, almond length = {L * k / (2 * np.pi):.2f} wavelengths", flush=True)

    # sanity check: data from an interior dipole, exact exterior field known
    y0, a0 = np.array([[0.01, 0.005, 0.002]]), np.array([[1.0, 0.5j, -0.7]])
    inc = lambda p: -np.einsum("pnab,nb->pa", dyadic(p, y0, k), a0)
    exact = -inc(pts)
    for h, n_s in ((0.004, 30), (0.0025, 40)):
        y, alpha, info = mfs(k, h, n_s, 0.5, inc=inc)
        print(f"dipole check h={h}: M={info['M']} N={info['N']} bc={boundary_residual(k, y, alpha, inc):.2e} "
              f"error on r=0.3: {np.linalg.norm(evaluate(pts, y, alpha, k) - exact) / np.linalg.norm(exact):.2e}",
              flush=True)

    results, fields = [], []
    configs = [(0.004, 30, 0.5), (0.003, 36, 0.5), (0.0022, 44, 0.5), (0.0022, 44, 0.35), (0.0022, 44, 0.7)]
    if "--fine" in sys.argv:
        configs.append((0.0017, 52, 0.5))
    for h, n_s, tau in configs:
        y, alpha, info = mfs(k, h, n_s, tau)
        info["bc"] = boundary_residual(k, y, alpha)
        info.update(h=h, n_s=n_s)
        E = evaluate(pts, y, alpha, k)
        results.append(info)
        fields.append(E)
        print(f"MFS h={h} n_s={n_s} tau={tau}: M={info['M']} N={info['N']} rank={info['rank']} "
              f"bc={info['bc']:.2e} |c|={info['coef_norm']:.1e} t={info['t_solve']:.1f}s", flush=True)

    best = fields[int(np.argmin([r["bc"] for r in results]))]
    for r, E in zip(results, fields):
        r["vs_best"] = float(np.linalg.norm(E - best) / np.linalg.norm(best))

    np.savez(f"almond_mfs_fields_{freq}GHz.npz", pts=pts, fields=np.array(fields))
    if "--nobem" in sys.argv:
        sys.exit()
    for n_t, n_s in ((40, 36), (60, 54)):
        t0 = time.time()
        Eb, ndof = bem(k, n_t, n_s, pts)
        d = float(np.linalg.norm(Eb - best) / np.linalg.norm(best))
        results.append({"bem_dofs": int(ndof), "bem_vs_best_mfs": d, "t": time.time() - t0})
        print(f"BEM dofs={ndof}: difference from best MFS on r=0.3: {d:.2e} ({time.time() - t0:.0f}s)", flush=True)

    json.dump(results, open(f"results_pec_almond_{freq}GHz.json", "w"), indent=1)
