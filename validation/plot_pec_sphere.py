"""Plot results_pec_sphere.json -> pec_sphere_convergence.png"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

res = json.load(open("results_pec_sphere.json"))
ks = sorted({r["k"] for r in res})
fig, axes = plt.subplots(len(ks), 3, figsize=(16, 4.6 * len(ks)), squeeze=False)

for row, k in enumerate(ks):
    rn = [r for r in res if r["k"] == k and r["sweep"] == "N"]
    rr = sorted((r for r in res if r["k"] == k and r["sweep"] == "R"), key=lambda r: r["R"])
    Rs = sorted({r["R"] for r in rn})

    # error vs N and vs sqrt(N)
    for ax, xf, xl in ((axes[row, 0], lambda n: n, "$N$"),
                       (axes[row, 1], np.sqrt, r"$\sqrt{N}$")):
        for R in Rs:
            s = sorted((r for r in rn if r["R"] == R), key=lambda r: r["N"])
            ax.semilogy([xf(r["N"]) for r in s], [r["mid"] for r in s], "o-", ms=4, label=f"$R={R}$")
        ax.set_xlabel(xl)
        ax.set_ylabel("relative error vs Mie on $r=1.5$")
        ax.set_title(f"$k={k / np.pi:.0f}\\pi$, $M=2N$")
        ax.grid(True, which="both", ls=":", alpha=0.6)
        ax.legend()

    ax = axes[row, 2]
    R = [r["R"] for r in rr]
    for key, lab in (("bc", "boundary residual"), ("mid", "error on $r=1.5$"), ("far", "far-field error")):
        ax.semilogy(R, [r[key] for r in rr], "o-", ms=4, label=lab)
    ax.semilogy(R, [r["alpha_norm"] for r in rr], "s--", ms=4, color="gray", label=r"$\|\alpha\|$")
    ax.set_xlabel("source radius $R$")
    ax.set_title(f"$k={k / np.pi:.0f}\\pi$, $N={rr[0]['N']}$, $M=2N$")
    ax.grid(True, which="both", ls=":", alpha=0.6)
    ax.legend()

fig.suptitle("MFS for PEC unit sphere, validated against the Mie series")
fig.tight_layout()
fig.savefig("pec_sphere_convergence.png", dpi=150)
print("wrote pec_sphere_convergence.png")
