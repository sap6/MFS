"""Where are the singularities of the continued scattered field of the torus?

Sources on a circle of radius b about (rho0, 0) in the meridian half-plane.
If the field continued into the tube is analytic outside a disc of radius
delta about (rho0, 0), the n-th Fourier coefficient of the dipole density
around the source circle decays like (delta / b)^n. The script prints delta
estimated from the ratio of successive coefficients over eight windows of n
(4-8, 9-13, ... , 39-43; the later windows reach the rounding floor when the
decay is fast) next to the distance from (rho0, 0) to the focal ring of the
toroidal coordinates, rho = c = sqrt(R^2 - a^2), z = 0.

Prediction (DERIVATIONS.md, section E): the singular set is the centre circle
rho = R and the focal ring rho = c, so delta = max(|rho0 - c|, |R - rho0|).

Usage: python3 torus_singularities.py
"""
import numpy as np

from pec_torus import Torus, circle, plane_wave, residual, solve, wave


def density_decay(tor, k, d, p, rho0, b, ns=96, nphi=192):
    rho, z, _ = circle(ns, rho0, b)
    rhs = lambda x: -plane_wave(x, k, d, p)
    s = solve(tor, k, rhs, rho, z, nphi, ntheta_col=4 * ns)
    amp = 0.0
    for m in np.where(s.rhs_modes > 1e-3 * s.rhs_modes.max())[0]:   # modes about the axis that the wave excites
        h = np.fft.fft(s.modes[m], axis=0) / ns
        amp = amp + np.sum(np.abs(h) ** 2, axis=1)
    v = np.sqrt(amp[1:ns // 2] + amp[:ns // 2:-1])                   # |n| = 1 .. ns/2 - 1
    r = v[1:] / v[:-1]
    est = [b * np.exp(np.mean(np.log(r[i:i + 5]))) for i in range(3, 43, 5)]
    return est, residual(tor, s, rhs), v[42] / v.max()


def main():
    cases = [(0.8, 2, wave(), 1.0, 0.60, "along the axis"),
             (0.8, 2, wave(90, 0, 0), 1.0, 0.60, "in the plane of the ring"),
             (0.8, 2, wave(60, 25, 35), 1.0, 0.60, "oblique"),
             (0.8, 4, wave(), 1.0, 0.60, "along the axis, k = 4 pi"),
             (0.8, 2, wave(), 1.1, 0.62, "along the axis"),
             (0.5, 2, wave(), 1.0, 0.30, "along the axis"),
             (0.25, 2, wave(), 1.0, 0.11, "along the axis")]
    print(" a/R  wave                      source circle   |rho0 - c|  a^2/2R   delta in eight windows of n                                 residual")
    for a, kk, (d, p), rho0, b, name in cases:
        tor = Torus(1.0, a)
        est, bc, last = density_decay(tor, kk * np.pi, d, p, rho0, b)
        print(f" {a:4.2f} {name:25s} ({rho0:4.2f}, {b:4.2f})    {abs(rho0 - tor.c):6.4f}    {a * a / 2:6.4f}   "
              + " ".join(f"{e:6.3f}" for e in est) + f"   {last:.0e}      {bc:.0e}")


if __name__ == "__main__":
    main()
