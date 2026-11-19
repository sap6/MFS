"""
Mie series: the exact field scattered by a perfectly conducting sphere of
radius a, for the incident wave p e^{ik d·x} (time dependence e^{-iωt}).
Notation of Bohren and Huffman (1983), Chapter 4. The series is summed for a
wave along +z polarised along x, and the result is rotated to (d, p).
"""
import numpy as np
from scipy.special import spherical_jn, spherical_yn


def radial(n, rho, hankel):
    """z_n(rho) and [rho z_n(rho)]'/rho, with z = j (Bessel) or h (Hankel, first kind)"""
    z, dz = spherical_jn(n, rho), spherical_jn(n, rho, derivative=True)
    if hankel:
        z, dz = z + 1j * spherical_yn(n, rho), dz + 1j * spherical_yn(n, rho, derivative=True)
    return z, z / rho + dz


def scattered_field(pts, k, d, p, a=1.0, nmax=None):
    """
    Scattered field of the sphere
      pts : points outside the sphere, shape (P, 3)
      a   : radius of the sphere; nmax: number of terms (default ka + 4.05 (ka)^{1/3} + 25)
      Returns an array of shape (P, 3)
    """
    d, p = np.asarray(d, float), np.asarray(p, float)
    Q = np.column_stack((p, np.cross(d, p), d))         # rotation with Q z = d and Q x = p
    nmax = nmax or int(np.ceil(k * a + 4.05 * (k * a) ** (1 / 3) + 25))
    x = pts @ Q

    # Spherical coordinates of the points in the frame of the wave
    r, phi = np.linalg.norm(x, axis=1), np.arctan2(x[:, 1], x[:, 0])
    ct = np.clip(x[:, 2] / r, -1, 1)
    st, sp, cp = np.sqrt(1 - ct ** 2), np.sin(phi), np.cos(phi)

    # Sum the series, with the angular functions pi_n and tau_n by recurrence
    Er, Et, Ep = (np.zeros(len(x), complex) for _ in range(3))
    pi0, pi1 = np.zeros_like(ct), np.ones_like(ct)
    for n in range(1, nmax + 1):
        tau = n * ct * pi1 - (n + 1) * pi0
        z, dz = radial(n, k * r, True)
        jx, djx = radial(n, k * a, False)
        hx, dhx = radial(n, k * a, True)
        En = 1j ** n * (2 * n + 1) / (n * (n + 1))
        cM, cN = -jx / hx * En, 1j * djx / dhx * En    # coefficients -b_n E_n and i a_n E_n
        Et += cM * cp * pi1 * z + cN * cp * tau * dz
        Ep += -cM * sp * tau * z - cN * sp * pi1 * dz
        Er += cN * cp * n * (n + 1) * st * pi1 * z / (k * r)
        pi0, pi1 = pi1, ((2 * n + 1) * ct * pi1 - (n + 1) * pi0) / n

    # Back to Cartesian components, and to the frame of (d, p)
    rhat = np.column_stack((st * cp, st * sp, ct))
    that = np.column_stack((ct * cp, ct * sp, -st))
    phat = np.column_stack((-sp, cp, 0 * sp))
    return (Er[:, None] * rhat + Et[:, None] * that + Ep[:, None] * phat) @ Q.T
