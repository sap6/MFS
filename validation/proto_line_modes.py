"""One line of the cube layout at three wavelengths: how well do the density
modes resolve an oscillating density? (NOTES.md item 67, task T30.)

The modes are Legendre polynomials in the stretched variable
u = log((H + t)/(H - t)) / log((2H - d)/d). Test densities on |t| <= H - d:
exp(i k s t) times a vertex factor ((H - t)^2 + d^2)^(a/2), for direction
cosines s and exponents a. Printed: relative L2(dt) error of the best fit
with n stretched modes, and with n - m stretched modes plus m Legendre
polynomials in t.

    python3 proto_line_modes.py
"""
import numpy as np

H = 0.5
LEG = np.polynomial.legendre


def fit(d, k, n_u, n_t, s, a, nq=4000):
    lam = np.log((2 * H - d) / d)
    g, w = LEG.leggauss(nq)                                  # Gauss nodes in u: fine near both ends
    t = H * np.tanh(0.5 * g * lam)
    wt = w * H * 0.5 * lam / np.cosh(0.5 * g * lam) ** 2
    f = np.exp(1j * k * s * t) * ((H - t) ** 2 + d ** 2) ** (a / 2)
    cols = [LEG.legvander(g, n_u - 1)]
    if n_t:
        cols.append(LEG.legvander(t / (H - d), n_t - 1))
    A = np.hstack(cols) * np.sqrt(wt)[:, None]
    b = f * np.sqrt(wt)
    c = np.linalg.lstsq(A, b, rcond=1e-13)[0]
    return np.linalg.norm(A @ c - b) / np.linalg.norm(b)


if __name__ == "__main__":
    for S in (1, 2, 3):
        k = 2 * np.pi * S
        print(f"side = {S} wavelengths")
        for d in (0.2, 0.05, 0.0125, 0.003125):
            n = 8 + int(np.ceil(2 * np.log(2 * H / d))) + 4 * (S - 1) + (4 if S == 3 else 0)
            for s, a in ((1.0, 0.0), (0.6, -0.3)):
                row = [f"  d = {d:8.5f}  n = {n:2d}  s = {s}  a = {a:4.1f}:  stretched only {fit(d, k, n, 0, s, a):.1e}"]
                for m in (4, 8, 12, 16):
                    if m < n - 4:
                        row.append(f"{n - m}+{m}: {fit(d, k, n - m, m, s, a):.1e}")
                print("   ".join(row))
