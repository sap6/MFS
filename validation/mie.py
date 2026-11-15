"""Mie series for plane-wave scattering by a perfectly conducting sphere.

Time convention e^{-i omega t}. Notation follows Bohren & Huffman,
"Absorption and Scattering of Light by Small Particles" (1983), Ch. 4.

The series is derived for the standard configuration (incident wave travelling
along +z, polarised along x) and rotated to an arbitrary direction d and
polarisation p.
"""
import numpy as np
from scipy.special import spherical_jn, spherical_yn


def _radial(n, rho, kind):
    """Return z_n(rho) and [rho z_n(rho)]'/rho for kind 'j' or 'h' (h^(1))."""
    z = spherical_jn(n, rho)
    dz = spherical_jn(n, rho, derivative=True)
    if kind == "h":
        z = z + 1j * spherical_yn(n, rho)
        dz = dz + 1j * spherical_yn(n, rho, derivative=True)
    return z, z / rho + dz


def _angular(nmax, mu):
    """pi_n(mu), tau_n(mu) for n = 0..nmax; arrays of shape (nmax+1, P)."""
    pi = np.zeros((nmax + 1,) + mu.shape)
    tau = np.zeros_like(pi)
    pi[1] = 1.0
    tau[1] = mu
    for n in range(2, nmax + 1):
        pi[n] = ((2 * n - 1) * mu * pi[n - 1] - n * pi[n - 2]) / (n - 1)
        tau[n] = n * mu * pi[n] - (n + 1) * pi[n - 1]
    return pi, tau


def default_nmax(x):
    return int(np.ceil(x + 4.05 * x ** (1 / 3) + 25))


def _series(pts, k, a, which, nmax):
    """Field in the standard frame. which: 'inc' (series form) or 'sc'."""
    r = np.linalg.norm(pts, axis=1)
    mu = np.clip(pts[:, 2] / r, -1.0, 1.0)
    theta = np.arccos(mu)
    phi = np.arctan2(pts[:, 1], pts[:, 0])
    st, ct, sp, cp = np.sin(theta), np.cos(theta), np.sin(phi), np.cos(phi)
    rho, x = k * r, k * a
    pi, tau = _angular(nmax, mu)

    Er = np.zeros(len(pts), complex)
    Et = np.zeros(len(pts), complex)
    Ep = np.zeros(len(pts), complex)
    for n in range(1, nmax + 1):
        En = 1j ** n * (2 * n + 1) / (n * (n + 1))
        if which == "inc":
            z, dz = _radial(n, rho, "j")
            cM, cN = En, -1j * En
        else:
            z, dz = _radial(n, rho, "h")
            jx, djx = _radial(n, x, "j")
            hx, dhx = _radial(n, x, "h")
            an = djx / dhx  # psi_n'(x) / xi_n'(x)
            bn = jx / hx    # psi_n(x)  / xi_n(x)
            cM, cN = -bn * En, 1j * an * En
        # M_o1n
        Et += cM * cp * pi[n] * z
        Ep += -cM * sp * tau[n] * z
        # N_e1n
        Er += cN * cp * n * (n + 1) * st * pi[n] * z / rho
        Et += cN * cp * tau[n] * dz
        Ep += -cN * sp * pi[n] * dz

    rhat = np.column_stack((st * cp, st * sp, ct))
    that = np.column_stack((ct * cp, ct * sp, -st))
    phat = np.column_stack((-sp, cp, np.zeros_like(sp)))
    return Er[:, None] * rhat + Et[:, None] * that + Ep[:, None] * phat


def _rotation(d, p):
    """Proper rotation Q with Q z = d and Q x = p."""
    d = np.asarray(d, float) / np.linalg.norm(d)
    p = np.asarray(p, float) / np.linalg.norm(p)
    assert abs(d @ p) < 1e-14, "polarisation must be orthogonal to direction"
    return np.column_stack((p, np.cross(d, p), d))


def scattered_field(pts, k, d, p, a=1.0, nmax=None):
    """Scattered E at pts (P,3), |pts| >= a, for E_inc = p exp(i k d.x)."""
    Q = _rotation(d, p)
    nmax = nmax or default_nmax(k * a)
    return _series(pts @ Q, k, a, "sc", nmax) @ Q.T


def incident_series(pts, k, d, p, nmax=None):
    """Series form of the incident plane wave (used only for self-checks)."""
    Q = _rotation(d, p)
    nmax = nmax or default_nmax(k * np.linalg.norm(pts, axis=1).max())
    return _series(pts @ Q, k, 1.0, "inc", nmax) @ Q.T


def plane_wave(pts, k, d, p):
    return np.asarray(p)[None, :] * np.exp(1j * k * (pts @ np.asarray(d)))[:, None]


if __name__ == "__main__":
    # Self-checks: (1) series reproduces the plane wave; (2) total tangential
    # field vanishes on the sphere; (3) radiation: |E_sc| ~ 1/r in the far field.
    rng = np.random.default_rng(0)
    k = 4 * np.pi
    d = np.array([0.0, 0.0, -1.0])
    p = np.array([0.0, 1.0, 0.0])
    u = rng.normal(size=(2000, 3))
    u /= np.linalg.norm(u, axis=1, keepdims=True)

    pts = 1.3 * u
    err = np.abs(incident_series(pts, k, d, p) - plane_wave(pts, k, d, p)).max()
    print(f"plane-wave expansion max error at r=1.3: {err:.2e}")

    tot = scattered_field(u, k, d, p) + plane_wave(u, k, d, p)
    tang = tot - np.sum(tot * u, axis=1, keepdims=True) * u
    print(f"max |tangential E_tot| on r=1:           {np.abs(tang).max():.2e}")

    e1 = np.linalg.norm(scattered_field(1e3 * u, k, d, p), axis=1)
    e2 = np.linalg.norm(scattered_field(2e3 * u, k, d, p), axis=1)
    print(f"far-field 1/r decay, max |2 E(2r)/E(r) - 1|: {np.abs(2 * e2 / e1 - 1).max():.2e}")
