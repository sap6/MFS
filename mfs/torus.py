"""
Scattering by the perfectly conducting torus (sqrt(x² + y²) - R)² + z² = a²
(Sections 3.2 and 3.3 of the paper). Electric dipoles on rings of an inner
torus with minor radius a_s, boundary points on rings of the surface.

Both point sets are invariant under the rotation by 2π/n_phi about the axis,
so the matrix is block-circulant and the least-squares problem splits into
n_phi small problems, one for each Fourier mode about the axis, for any
incident wave (as for acoustic scattering in Liu and Barnett 2016).
"""
import numpy as np
from scipy.linalg import lstsq

from .core import cross_sections, dyadic, field, plane_wave, sphere_rule


def circle(n, centre, radius):
    """n points on a circle in the meridian half-plane: rho, z and the angle"""
    t = 2 * np.pi * (np.arange(n) + 0.5) / n
    return centre + radius * np.cos(t), radius * np.sin(t), t


def ring_points(rho, z, nphi):
    """Rings of nphi points about the z-axis, shape (len(rho), nphi, 3), and the angles"""
    ph = 2 * np.pi * np.arange(nphi) / nphi
    x = np.stack((rho[:, None] * np.cos(ph), rho[:, None] * np.sin(ph), z[:, None] + 0 * ph), axis=2)
    return x, ph


def frame(ph):
    """Unit vectors e_rho, e_phi, e_z at the angles ph, shape (len(ph), 3, 3)"""
    c, s, o = np.cos(ph), np.sin(ph), 0 * ph
    return np.stack((np.stack((c, s, o), 1), np.stack((-s, c, o), 1), np.stack((o, o, o + 1), 1)), 1)


class Torus:
    """The torus (sqrt(x² + y²) - R)² + z² = a²: boundary points and test points"""
    def __init__(self, R=1.0, a=0.5):
        self.R, self.a = R, a

    def surface(self, ntheta, nphi):
        """
        Boundary points on ntheta rings of nphi points
          Returns points (ntheta, nphi, 3), tangents (ntheta, nphi, 2, 3) and area weights (ntheta,)
        """
        rho, z, t = circle(ntheta, self.R, self.a)
        x, ph = ring_points(rho, z, nphi)
        F = frame(ph)
        t_pol = -np.sin(t)[:, None, None] * F[None, :, 0] + np.cos(t)[:, None, None] * F[None, :, 2]
        t_tor = np.broadcast_to(F[None, :, 1], t_pol.shape)
        w = self.a * rho * (2 * np.pi / ntheta) * (2 * np.pi / nphi)
        return x, np.stack((t_pol, t_tor), axis=2), w

    def test_points(self, n=1500, seed=3):
        """Random points, uniform in area, and their tangents (to measure the boundary residual)"""
        rng = np.random.default_rng(seed)
        t = np.empty(0)
        while len(t) < n:                          # rejection sampling: density proportional to rho
            c = rng.uniform(0, 2 * np.pi, 2 * n)
            t = np.concatenate((t, c[rng.uniform(0, self.R + self.a, 2 * n) < self.R + self.a * np.cos(c)]))
        t, ph = t[:n], rng.uniform(0, 2 * np.pi, n)
        rho = self.R + self.a * np.cos(t)
        x = np.column_stack((rho * np.cos(ph), rho * np.sin(ph), self.a * np.sin(t)))
        F = frame(ph)
        return x, np.stack((-np.sin(t)[:, None] * F[:, 0] + np.cos(t)[:, None] * F[:, 2], F[:, 1]), axis=1)


class Solution:
    """Dipole positions y (N, 3) and moments alpha (N, 3)"""

    def __init__(self, y, alpha, k):
        self.y, self.alpha, self.k = y, alpha, k

    def field(self, pts):
        """Scattered field at pts (P, 3)"""
        return field(pts, self.y, self.alpha, self.k)

    def far(self, dirs):
        """Far-field pattern F, with E ~ F e^{ikr}/r, in the directions dirs"""
        s = np.exp(-1j * self.k * (dirs @ self.y.T)) @ self.alpha
        return s - dirs * np.sum(dirs * s, axis=1, keepdims=True)


def solve(torus, k, rhs, ns, nphi, a_s, rings=2):
    """
    Least squares, one small problem per Fourier mode about the axis
      rhs  : function x (M, 3) → the field to match on the boundary, (M, 3) (minus the incident field for scattering)
      ns   : number of source rings, nphi: dipoles per ring (3 ns nphi unknowns)
      a_s  : minor radius of the source torus
      rings: boundary rings per source ring (rings · ns rings of nphi points, 2 rings ns nphi equations)
    """
    rho, z, _ = circle(ns, torus.R, a_s)
    nc = rings * ns
    x, T, w = torus.surface(nc, nphi)
    sw = np.sqrt(w)

    # Blocks of one source per ring (at angle 0) against all boundary points
    y0 = np.column_stack((rho, 0 * rho, z))
    G = dyadic(x.reshape(-1, 3), y0, k).reshape(nc, nphi, ns, 3, 3)
    B = np.einsum("ipca,ipjab->ipcjb", T, G) * sw[:, None, None, None, None]
    b = np.einsum("ipca,ipa->ipc", T, rhs(x.reshape(-1, 3)).reshape(nc, nphi, 3)) * sw[:, None, None]

    # Fourier transform over the rotation, then one least-squares problem per mode
    Bh, bh = np.fft.fft(B, axis=1), np.fft.fft(b, axis=1)
    ah = np.array([lstsq(Bh[:, m].reshape(2 * nc, 3 * ns), bh[:, m].reshape(-1))[0].reshape(ns, 3)
                   for m in range(nphi)])

    # Back to the dipoles: alpha_q = (1/nphi) Σ_m ah_m e^{2πi mq/nphi}, from cylindrical to Cartesian components
    y, ph = ring_points(rho, z, nphi)
    alpha = np.einsum("qjc,qca->jqa", np.fft.ifft(ah, axis=0), frame(ph))
    return Solution(y.reshape(-1, 3), alpha.reshape(-1, 3), k)


def scatter(k, d, p, ns, nphi, R=1.0, a=0.5, a_s=None, rings=2):
    """Plane-wave scattering with sources on the inner torus of minor radius a_s (default a / 2)"""
    torus = Torus(R, a)
    return solve(torus, k, lambda x: -plane_wave(x, k, d, p), ns, nphi, 0.5 * a if a_s is None else a_s, rings)


def residual(sol, d, p, R=1.0, a=0.5, n=1500):
    """Relative tangential boundary residual at random points of the surface"""
    xt, Tt = Torus(R, a).test_points(n)
    got = np.einsum("mca,ma->mc", Tt, sol.field(xt))
    want = -np.einsum("mca,ma->mc", Tt, plane_wave(xt, sol.k, d, p))
    return float(np.linalg.norm(got - want) / np.linalg.norm(want))


def optical(sol, d, p):
    """Extinction cross-section and the optical-theorem mismatch |ext - sca| / sca"""
    ext, sca = cross_sections(sol.far, sol.k, d, p, sphere_rule(48, 96))
    return ext, abs(ext - sca) / sca
