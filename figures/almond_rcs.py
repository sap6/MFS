"""
Figure 9 and Table 3: monostatic radar cross-section of the NASA almond against
the azimuth (0 = nose-on), HH and VV, at 3.5, 5.125, 7 and 10.25 GHz: the MFS
(layout 6:3), and the measurements and reference simulation of the Austin RCS
Benchmark Suite (Kelley et al. 2020, problem III-A; data in data/austin_almond).

    python figures/almond_rcs.py          # the figure and the table, from data/
    python figures/almond_rcs.py --run    # solve again (2 to 6 minutes per frequency)
"""
import sys

import numpy as np
from style import BLUE, ORANGE, data, finish, plt

from mfs import almond

FREQS = {"x1": 3.5, "x2": 5.125, "x3": 7.0, "x4": 10.25}     # the suite's file names and frequencies [GHz]


def austin(kind, fx, pol):
    """The suite's curve (azimuth, dBsm) on the half-degree grid 0 to 180"""
    d = np.loadtxt(data(f"austin_almond/{kind}_rcs.III.A.s1.f{fx}.{pol[0]}.txt"))
    phi, db = np.round(d[:, 2], 2), d[:, 3]
    keep = np.isclose(phi * 2, np.round(phi * 2))
    return phi[keep], db[keep]


if "--run" in sys.argv:
    for f in FREQS.values():
        np.savez(data(f"almond_rcs_{f}GHz.npz"), **almond.rcs(f))

# Mean absolute difference in dB from the reference simulation, after raising both curves
# to the suite's threshold (largest reference value - 80 dB)
fig, axes = plt.subplots(2, 4, figsize=(6.6, 3.7), sharex=True, sharey=True)
print("frequency  pol  mean |MFS - reference| [dB]")
for j, (fx, f) in enumerate(FREQS.items()):
    res = np.load(data(f"almond_rcs_{f}GHz.npz"))
    for i, pol in enumerate(("HH", "VV")):
        mfs = 10 * np.log10(res["sigma_" + pol])
        phi, ref = austin("ref", fx, pol)
        th = ref.max() - 80.0
        print(f"{f:8.3f}   {pol}   {np.mean(np.abs(np.maximum(mfs, th) - np.maximum(ref, th))):.3f}")
        ax = axes[i, j]
        ax.plot(phi, austin("meas.m1", fx, pol)[1], "-", color="#9a988f", lw=0.9, label="measured (Kelley et al. 2020)", zorder=1)
        ax.plot(phi, mfs, "-", color=BLUE, lw=1.3, label="MFS (this work)", zorder=3)
        ax.plot(phi[::12], ref[::12], "o", color=ORANGE, ms=3.2, mfc="white", mew=0.9, label="reference simulation", zorder=4)
        ax.set_title(f"{f:g} GHz, {pol}", loc="left", fontsize=8.5)
        ax.set_xticks([0, 60, 120, 180])
        if i == 1:
            ax.set_xlabel("azimuth [deg]")
for ax in axes[:, 0]:
    ax.set_ylabel("RCS [dBsm]")
h, lab = axes[0, 0].get_legend_handles_labels()
fig.legend(h, lab, frameon=False, fontsize=8.2, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.005))
finish(fig, axes, "fig_almond_benchmark", rect=(0, 0.06, 1, 1))
