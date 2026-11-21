"""
Quick checks of the package against exact solutions and against the stored
results of the paper (data/). About a minute, most of it compiling the kernels:
python -m pytest tests
"""
import os

import numpy as np
import pytest

from mfs import almond, cone, core, cube, mie, sphere, tetrahedron, torus

DATA = os.path.join(os.path.dirname(__file__), "..", "data")
K = 2 * np.pi
D, P = core.wave()                                       # the standard wave: along -z, polarised along y
DW, PW = core.wave(60, 25, 35)                           # the general wave of Section 5.5


def test_mie():
    # The total tangential field vanishes on the sphere, for any wave
    u = core.fibonacci(500)
    for d, p in ((D, P), (DW, PW)):
        tot = mie.scattered_field(u, K, d, p) + core.plane_wave(u, K, d, p)
        assert np.abs(tot - np.sum(tot * u, axis=1)[:, None] * u).max() < 1e-12


def test_sphere():
    # Section 3.1: 400 dipoles on the sphere of radius 0.5, error 2.0e-9 against the Mie series
    pts = 2.0 * core.fibonacci(400)
    y, alpha = sphere.solve(400)
    assert core.rel(core.field(pts, y, alpha, K), mie.scattered_field(pts, K, D, P)) < 2.1e-9


def test_torus():
    # A field with a known exterior: dipoles inside the tube, away from the source torus
    yk = np.array([[1.08, 0.0, 0.05], [-0.3, 0.95, -0.06]])
    mk = np.array([[1.0, 0.5j, -0.3], [0.2, 1.0, 1.0j]])
    known = lambda x: np.einsum("pnab,nb->pa", core.dyadic(x, yk, K), mk)
    pts = 2.5 * core.sphere_rule(8, 16)[0]
    errors = [core.rel(torus.solve(torus.Torus(), K, known, ns, nphi, 0.25).field(pts), known(pts))
              for ns, nphi in ((20, 64), (28, 96))]
    assert errors[1] < 1e-6 and errors[1] < 1e-2 * errors[0]


def test_torus_residual_nondefault_a():
    # torus.residual uses sol.R and sol.a rather than defaulting to a = 0.5
    sol = torus.scatter(K, D, P, 24, 96, a=0.8, a_s=0.52)
    assert torus.residual(sol, D, P) < 1e-3 and torus.residual(sol, D, P, a=0.5) > 0.5


def test_cube_against_paper():
    # Section 5.4: the level J = 4 reproduces the stored field of the paper
    d = np.load(os.path.join(DATA, "cube_convergence.npz"))
    E, t, info = cube.solve(4, pts=d["pts"])
    assert core.rel(E, d["dense_E4"]) < 1e-10
    assert info["unknowns"] == [1128, 1128]


def test_cube_wavenumbers_and_empty_pts():
    # Unified mode and boundary rules at k = 4pi, 6pi match stored counts, and empty target arrays return cleanly
    for case, kw in (("k2", 4 * np.pi), ("k3", 6 * np.pi)):
        row = [r for r in np.load(os.path.join(DATA, f"cube_cases_{case}.npz"))["mfs"] if r[0] == 2 and r[1] == 5][0]
        c = cube.Cube(5, k=kw, classes=())
        assert c.unknowns[(1, -1, 1)] == int(row[3]) // 2 and c.rows == int(row[4])
    assert cube.solve(2, pts=np.empty((0, 3)))[0].shape == (0, 3)


def test_cube_rotation():
    # A general wave: four class matrices and the rotation give the same field as all eight classes
    pts = 2.0 * core.fibonacci(50)
    E, _, info = cube.solve(2, DW, PW, pts=pts)
    full = cube.Cube(2, classes=core.classes_of(K, DW, PW, cube.REFL))
    coef, _, _ = full.solve(DW, PW)
    assert len(info["unknowns"]) == 4 and len(coef) == 8
    assert core.rel(E, full.field(coef, pts)) < 1e-10


def test_almond_against_paper():
    # Section 6.3: layout 3:0 at 1.19 GHz, nose-on (with tip tail = 815 unknowns), reproduces stored field
    d = np.load(os.path.join(DATA, "cost_almond.npz"))
    k = 2 * np.pi * 1.19e9 / almond.C0
    src, c, E, shape = almond.nose_on(k, 3, 0, 0.3 * core.fibonacci(400))
    assert isinstance(src, tuple) and len(src) == 2 and shape[1] == 815
    assert core.rel(E, d["mfs_E0"]) < 1e-10
    src4, c4, _, _ = almond.nose_on(k, 4, 1, 0.3 * core.fibonacci(10))
    ext, sca = almond.optical(k, src4, c4)
    assert abs(ext - sca) / sca < 2e-2
    r = almond.rcs(3.5, 4, 1, step=90.0)
    ref = np.load(os.path.join(DATA, "almond_rcs_3.5GHz.npz"))
    for pol in ("sigma_HH", "sigma_VV"):
        assert np.abs(10 * np.log10(r[pol] / ref[pol][[0, 180, 360]])).max() < 0.25


def test_csne_qr_fallback():
    # Well-conditioned C-ordered CSNE and ill-conditioned (cond > 1e8) automatic CSNE -> QR fallback
    rng = np.random.default_rng(0)
    Ac = np.ascontiguousarray(rng.normal(size=(30, 8)) + 1j * rng.normal(size=(30, 8)))
    bc = rng.normal(size=30) + 1j * rng.normal(size=30)
    assert core.rel(core.CSNE(Ac).solve(bc), core.QR(Ac.copy()).solve(bc)) < 1e-13
    U, _ = np.linalg.qr(rng.normal(size=(40, 12)) + 1j * rng.normal(size=(40, 12)))
    V, _ = np.linalg.qr(rng.normal(size=(12, 12)) + 1j * rng.normal(size=(12, 12)))
    s = np.logspace(0, -9, 12)
    A = np.asfortranarray((U * s) @ V.conj().T)
    b = U @ np.ones(12, complex) + 1e-4 * rng.normal(size=40)
    sol = core.factor(A.copy(), overwrite=True)
    assert isinstance(sol, core.QR)
    assert core.rel(sol.solve(b), core.QR(A.copy()).solve(b)) < 1e-12


def test_line_tree_cache_cap():
    # Sweeping multiple cube levels keeps lines._TREES bounded in memory
    from mfs import lines
    for J in (3, 4, 5, 6):
        cube.Cube(J, classes=[(1, -1, 1)])
    assert len(lines._TREES) <= 32


def test_tetrahedron_and_rotation():
    # Section 5.6: level J = 4 reproduces the stored field, and 2-matrix ROT matches all 4 classes
    d = np.load(os.path.join(DATA, "tetrahedron_convergence.npz"))
    E, _, info = tetrahedron.solve(4, pts=d["pts"])
    assert core.rel(E, d["dense_E4"]) < 1e-10
    assert info["unknowns"] == [1376]
    assert tetrahedron.solve(2, pts=np.empty((0, 3)))[0].shape == (0, 3)
    pts = 2.0 * core.fibonacci(50)
    Er, _, ir = tetrahedron.solve(2, DW, PW, pts=pts)
    full = tetrahedron.Tetrahedron(2, classes=core.classes_of(K, DW, PW, tetrahedron.REFL))
    cf, _, _ = full.solve(DW, PW)
    assert len(ir["unknowns"]) == 2 and len(cf) == 4
    assert core.rel(Er, full.field(cf, pts)) < 1e-10


def test_cone():
    # Appendix A: smallest electric exponents at the almond tip, cube vertex, and tetrahedron vertex
    nu = cone.EllipticCone().exponents(1, 1, hi=0.5)
    assert abs(nu[0] - 0.29959) < 1e-5
    nu_cube = cone.CubeVertex(n=10).exponents()
    assert abs(nu_cube[0] - 0.45417) < 2e-4
    nu_tet = cone.CubeVertex(n=10, V=cone.TET_EDGES).exponents()
    assert abs(nu_tet[0] - 0.32221) < 2e-4


@pytest.mark.skipif(not os.environ.get("MFS_TEST_BEMPP"), reason="set MFS_TEST_BEMPP=1 to test the Bempp comparison")
def test_bempp_symmetry():
    # Symmetric boundary element solves (Cyclic torus, Mirror cube, tetrahedron & almond), F-ordered fold, chunked assembly
    from mfs import bem
    F = bem.fold(np.ones((4, 4), complex), np.array([0, 1]), np.array([[0, 1], [2, 3]]), np.ones((2, 2)))
    assert F.flags.f_contiguous
    pts = 2.0 * core.fibonacci(50)
    E_sym, _ = bem.torus(8, K, D, P, pts)
    grid_t = bem.torus_grid(1.0, 0.5, 8)[0]
    E_full, _ = bem.efie(grid_t, K, D, P, pts)
    assert core.rel(E_sym, E_full) < 1e-10
    div_t = bem.bempp.function_space(grid_t, "RWG", 0)
    assert core.rel(bem.assemble_rows(div_t, K, np.array([0, 3, 7]), max_gb=1e-6), bem.assemble_rows(div_t, K, np.array([0, 3, 7]))) < 1e-14
    Ec_sym, _ = bem.cube(4, 1.5, K, DW, PW, pts)
    Ec_full, _ = bem.efie(bem.cube_grid(4, 1.5)[0], K, DW, PW, pts)
    assert core.rel(Ec_sym, Ec_full) < 5e-5
    Et_sym, _ = bem.tetrahedron(6, 1.8, K, DW, PW, pts)
    Et_full, _ = bem.efie(bem.tetrahedron_grid(6, 1.8)[0], K, DW, PW, pts)
    assert core.rel(Et_sym, Et_full) < 1e-12
    ka, pt_a = 2 * np.pi * 1.19e9 / almond.C0, 0.3 * core.fibonacci(30)
    Ea_sym, _ = bem.almond(12, 12, ka, DW, PW, pt_a)
    q = 3
    th, phi = np.arange(1, 12) / 12 * np.pi, np.arange(q + 1) / 12 * 2 * np.pi
    xs = almond.L * (0.5 * (almond.T0 + almond.T1) - 0.5 * (almond.T1 - almond.T0) * np.cos(th))
    a, b = almond.semi_axes(xs)
    rings = np.stack((np.repeat(xs[:, None], q + 1, 1), a[:, None] * np.cos(phi), b[:, None] * np.sin(phi)), axis=-1).reshape(-1, 3)
    verts = np.concatenate(([[almond.L * almond.T0, 0, 0]], rings, [[almond.L * almond.T1, 0, 0]]))
    verts[np.abs(verts) < 1e-15] = 0.0
    idx = lambda i, j: 1 + i * (q + 1) + j
    last, tris = len(verts) - 1, []
    for j in range(q):
        tris += [(0, idx(0, j + 1), idx(0, j)), (last, idx(10, j), idx(10, j + 1))]
        for i in range(10):
            tris += [(idx(i, j), idx(i, j + 1), idx(i + 1, j)), (idx(i, j + 1), idx(i + 1, j + 1), idx(i + 1, j))]
    Ea_full, _ = bem.efie(bem.mirror_grid(verts, np.array(tris), bem.group((0, 1, 1)))[0], ka, DW, PW, pt_a)
    assert core.rel(Ea_sym, Ea_full) < 5e-5
