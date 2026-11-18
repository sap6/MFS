"""
Quick checks of the package against exact solutions and against the stored
results of the paper (data/). About two minutes: python -m pytest tests
"""
import os

import numpy as np
import pytest

from mfs import almond, cone, core, cube, mie, sphere, torus

DATA = os.path.join(os.path.dirname(__file__), "..", "data")
K = 2 * np.pi
D, P = core.wave()                                       # the standard wave: along -z, polarised along y


def test_mie():
    # The total tangential field vanishes on the sphere, for any wave
    u = core.fibonacci(500)
    for d, p in ((D, P), core.wave(60, 25, 35)):
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


def test_cube_against_paper():
    # Section 5.4: the level J = 4 reproduces the stored field of the paper
    d = np.load(os.path.join(DATA, "cube_convergence.npz"))
    E, t, info = cube.solve(4, pts=d["pts"])
    assert core.rel(E, d["dense_E4"]) < 1e-12
    assert info["unknowns"] == [1128, 1128]


def test_cube_rotation():
    # A general wave: four class matrices and the rotation give the same field as all eight classes
    d, p = core.wave(60, 25, 35)
    pts = 2.0 * core.fibonacci(50)
    E, _, info = cube.solve(2, d, p, pts=pts)
    full = cube.Cube(2, classes=cube.classes_of(K, d, p))
    coef, _, _ = full.solve(d, p)
    assert len(info["unknowns"]) == 4 and len(coef) == 8
    assert core.rel(E, full.field(coef, pts)) < 1e-10


def test_almond_against_paper():
    # Section 6.3: the layout 3:0 at 1.19 GHz, nose-on, reproduces the stored field
    d = np.load(os.path.join(DATA, "cost_almond.npz"))
    k = 2 * np.pi * 1.19e9 / almond.C0
    E = almond.nose_on(k, 3, 0, 0.3 * core.fibonacci(400))[2]
    assert core.rel(E, d["mfs_E0"]) < 1e-12


def test_cone():
    # Appendix A: the smallest exponent at the tip of the almond
    nu = cone.EllipticCone().exponents(1, 1, hi=0.5)
    assert abs(nu[0] - 0.29959) < 1e-5


@pytest.mark.skipif(not os.environ.get("MFS_TEST_BEMPP"), reason="set MFS_TEST_BEMPP=1 to test the Bempp comparison")
def test_bempp_symmetry():
    # The symmetric boundary element solve equals the full one (Appendix C)
    import bempp_cl.api as bempp
    from mfs import bem
    pts = 2.0 * core.fibonacci(50)
    E_sym, N, _ = bem.torus(8, K, D, P, pts)
    E_full, _ = bem.efie(bem.torus_grid(1.0, 0.5, 8)[0], K, D, P, pts)
    assert core.rel(E_sym, E_full) < 1e-10
