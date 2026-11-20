"""
Figure 1: the unit sphere at k = 2π, error against the Mie series.
Left: against √N for four source radii. Right: against the source radius for N = 400.

    python figures/sphere.py          # the figure, from data/sphere_convergence.json
    python figures/sphere.py --run    # solve again (about two minutes)

The stored runs used the SVD; the package uses QR, which gives the same errors except
where rounding dominates (R_s <= 0.2 with N = 400), where they differ by up to a factor 2.
"""
import json
import sys

import numpy as np
from style import BLUE, GREEN, MUTED, ORANGE, PURPLE, curve, data, finish, plt, rel

from mfs import core, mie, sphere

K = 2 * np.pi

if "--run" in sys.argv:
    # Errors on the sphere of radius 1.5, the boundary residual, and the coefficient norm
    d, p = core.wave()
    u = core.fibonacci(1500)
    exact, T = mie.scattered_field(1.5 * u, K, d, p), core.tangent_basis(u)
    inc = np.einsum("mab,mb->ma", T, core.plane_wave(u, K, d, p))
    res = []
    Ns = (50, 100, 150, 200, 300, 400, 500, 600, 800, 1000, 1200)
    runs = [("N", R, N) for R in (0.3, 0.5, 0.7, 0.9) for N in Ns] + [("R", R, 400) for R in np.round(np.linspace(0.1, 0.95, 18), 3)]
    for sweep, R, N in runs:
        y, alpha = sphere.solve(N, K, R)
        res.append(dict(k=K, sweep=sweep, R=float(R), N=N, mid=rel(core.field(1.5 * u, y, alpha, K), exact),
                        bc=rel(np.einsum("mab,mb->ma", T, core.field(u, y, alpha, K)), -inc),
                        alpha_norm=float(np.linalg.norm(alpha))))
        print(sweep, R, N, f"{res[-1]['mid']:.2e}", flush=True)
    json.dump(res, open(data("sphere_convergence.json"), "w"), indent=1)

# The figure, from the k = 2π runs
res = [r for r in json.load(open(data("sphere_convergence.json"))) if abs(r["k"] - K) < 1e-9]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.5, 2.7))
for (col, mk, ls), R in zip(((BLUE, "o", "-"), (ORANGE, "D", "-."), (GREEN, "s", "--"), (PURPLE, "^", ":")), (0.3, 0.5, 0.7, 0.9)):
    s = sorted((r for r in res if r["sweep"] == "N" and r["R"] == R), key=lambda r: r["N"])
    curve(ax1, [np.sqrt(r["N"]) for r in s], [r["mid"] for r in s], col, mk, ls, f"$R_s = {R}$", 4.5, 1.4, 0.6, logx=False)
ax1.set(xlabel=r"$\sqrt{N}$,  $N$ = number of dipoles ($3N$ unknowns)", ylabel="relative error against Mie")
ax1.set_title(r"(a) $k = 2\pi$: error against $\sqrt{N}$", loc="left", fontsize=8.5)
ax1.legend(frameon=False, fontsize=8, loc="lower left")

s = sorted((r for r in res if r["sweep"] == "R"), key=lambda r: r["R"])
R = [r["R"] for r in s]
curve(ax2, R, [r["mid"] for r in s], BLUE, "o", "-", "error of the field", 4.5, 1.4, 0.6, logx=False)
curve(ax2, R, [r["bc"] for r in s], ORANGE, "D", "-.", "boundary residual", 4.5, 1.4, 0.6, logx=False)
curve(ax2, R, [r["alpha_norm"] for r in s], MUTED, "s", ":", r"coefficient norm $\|\boldsymbol{\alpha}\|$", 4, 1.2, 0.5, logx=False)
ax2.set(xlabel="source radius $R_s$", ylabel=r"error, residual, $\|\boldsymbol{\alpha}\|$")
ax2.set_title(rf"(b) $k = 2\pi$, $N = {s[0]['N']}$: dependence on $R_s$", loc="left", fontsize=8.5)
ax2.legend(frameon=False, fontsize=7.8, loc="upper right")
finish(fig, (ax1, ax2), "fig_sphere")
