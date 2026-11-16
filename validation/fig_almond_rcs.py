"""Monostatic RCS of the NASA almond against azimuth at 1.19, 7 and 9.92 GHz,
HH and VV: MFS (lines, half-degree steps) against Bempp (markers, 10-degree
steps, from bem_almond_rcs_ref.json).

    python3 fig_almond_rcs.py run    # compute the three MFS curves (about 2 min each), then plot
    python3 fig_almond_rcs.py        # plot and print the differences from saved data
"""
import json
import sys

import numpy as np

import almond_rcs

FREQS = (1.19, 7.0, 9.92)
# colour, dash and marker per polarisation: never colour alone (greyscale print)
STYLE = {"HH": ("#2a78d6", "-", "o"), "VV": ("#eb6834", "--", "s")}
INK, MUTED = "#1f1f1e", "#6b6a64"


def run():
    for f in FREQS:
        print(f"{f} GHz", flush=True)
        np.savez(f"fig_almond_rcs_{f}GHz.npz", **almond_rcs.rcs(f, 6, 3, step=0.5, s=0.5))


def plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.6})
    ref = json.load(open("bem_almond_rcs_ref.json"))
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.3), sharey=True)
    for ax, f in zip(axes, FREQS):
        d = np.load(f"fig_almond_rcs_{f}GHz.npz")
        r = ref[str(f)]
        idx = [int(np.argmin(np.abs(d["phi_deg"] - p))) for p in ref["phi_deg"]]
        for pol, (col, ls, mk) in STYLE.items():
            ax.plot(d["phi_deg"], 10 * np.log10(d["sigma_" + pol]), ls, color=col, lw=1.6, label=f"{pol}, MFS")
            ax.plot(ref["phi_deg"], 10 * np.log10(r[pol]), mk, color=col, ms=5, mfc="white", mew=1.2,
                    label=f"{pol}, Bempp")
            db = 10 * np.log10(d["sigma_" + pol][idx] / np.array(r[pol]))
            print(f"{f} GHz {pol}: MFS - Bempp ({r['ndof']} dofs) rms {np.sqrt(np.mean(db ** 2)):.3f} dB, max {np.abs(db).max():.3f} dB")
        ax.set_title(f"{f:g} GHz", loc="left", fontsize=9, color=INK)
        ax.set_xlabel("azimuth [deg]  (0 = nose-on)")
        ax.set_xticks(range(0, 181, 30))
        ax.grid(True, color="#e4e3df", lw=0.6)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel("monostatic RCS [dBsm]")
    axes[0].legend(frameon=False, fontsize=8, ncol=2, loc="lower center")
    fig.tight_layout()
    fig.savefig("fig_almond_rcs.pdf")
    fig.savefig("fig_almond_rcs.png", dpi=200)
    print("wrote fig_almond_rcs.pdf/.png")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        run()
    plot()
