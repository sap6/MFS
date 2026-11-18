"""NASA almond (L = 9.936 in, perfectly conducting): MFS monostatic RCS against
the public Austin RCS Benchmark Suite, problem set III-A: the Rye Canyon
measurements (Kelley, Chamulak, Courtney, Yilmaz, IEEE Antennas Propag. Mag.
62(1), 2020, doi 10.1109/MAP.2019.2955702) and the suite's reference
simulation (ARCHIE-AIM, finest mesh), at 3.5, 5.125, 7 and 10.25 GHz, HH and
VV, azimuth 0 to 180 degrees.

    python3 fig_almond_benchmark.py run    # compute the four MFS curves (2-4 min each), then report
    python3 fig_almond_benchmark.py        # table and figure from saved data

The data files are in austin_almond/ (unmodified copies, with the suite's
licence files; commit 300010c of github.com/UTAustinCEMGroup/AustinCEMBenchmarks).
Conventions agree with almond_rcs.py: theta = 90 degrees, phi = azimuth from
the tip, HH = phi-phi, VV = theta-theta. Our almond is 0.2524 m long, the
suite's 9.936 in = 0.25237 m (0.01% shorter).

The error measure is the suite's: mean over the angles of |difference in dB|,
after raising both curves to the threshold (largest reference value - 80 dB).
"m1" is the measurement of the 9.936 in almond, "m2" that of the almond twice
as large at half the frequency, shifted by -10 log10(4) dB.
"""
import sys

import numpy as np

import almond_rcs

FREQS = {"x1": 3.5, "x2": 5.125, "x3": 7.0, "x4": 10.25}
DATA = "austin_almond/"
INK, MUTED = "#1f1f1e", "#6b6a64"
BLUE, ORANGE, GREY = "#2a78d6", "#eb6834", "#9a988f"


def load(kind, fx, pol):
    """(phi, dB) on the half-degree grid 0..180."""
    d = np.loadtxt(f"{DATA}{kind}_rcs.III.A.s1.f{fx}.{pol[0]}.txt")
    phi, db = np.round(d[:, 2], 2), d[:, 3]
    keep = np.isclose(phi * 2, np.round(phi * 2))
    return phi[keep], db[keep]


def suite_error(db, ref):
    th = ref.max() - 80.0
    return np.mean(np.abs(np.maximum(db, th) - np.maximum(ref, th)))


def run(layouts=((6, 3, 0.5),)):
    for fx, f in FREQS.items():
        for n_mid, J, s in layouts:
            print(f"{f} GHz, layout {n_mid}:{J}, spacing {s}", flush=True)
            out = almond_rcs.rcs(f, n_mid, J, step=0.5, s=s)
            name = f"fig_almond_benchmark_{f}GHz" + ("" if (n_mid, J, s) == (6, 3, 0.5) else f"_{n_mid}_{J}_s{s}")
            np.savez(name + ".npz", **out)


def report():
    print("freq[GHz] pol   MFS-reference: mean|dB|  max|dB| at[deg]   MFS-m1   MFS-m2   reference-m1   reference-m2   (mean |dB|)")
    curves = {}
    for fx, f in FREQS.items():
        d = np.load(f"fig_almond_benchmark_{f}GHz.npz")
        for pol in ("HH", "VV"):
            mfs = 10 * np.log10(d["sigma_" + pol])
            phi, ref = load("ref", fx, pol)
            assert np.allclose(phi, d["phi_deg"])
            m1, m2 = load("meas.m1", fx, pol)[1], load("meas.m2", fx, pol)[1]
            i = int(np.argmax(np.abs(mfs - ref)))
            print(f"{f:8.3f}  {pol}   {suite_error(mfs, ref):14.3f} {np.abs(mfs - ref)[i]:9.2f} {phi[i]:7.1f}   "
                  f"{suite_error(mfs, m1):7.2f} {suite_error(mfs, m2):8.2f} {suite_error(ref, m1):12.2f} {suite_error(ref, m2):14.2f}")
            curves[(f, pol)] = (phi, mfs, ref, m1, m2)
    return curves


def plot(curves, fname="fig_almond_benchmark"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.6})
    fig, axes = plt.subplots(2, 4, figsize=(11.0, 5.2), sharex=True, sharey=True)
    for j, f in enumerate(FREQS.values()):
        for i, pol in enumerate(("HH", "VV")):
            ax = axes[i, j]
            phi, mfs, ref, m1, m2 = curves[(f, pol)]
            # measurement: neutral context; simulations: colour plus line style / marker (never colour alone)
            ax.plot(phi, m1, "-", color=GREY, lw=1.0, label="measured (Kelley et al. 2020)", zorder=1)
            ax.plot(phi, mfs, "-", color=BLUE, lw=1.6, label="MFS (this work)", zorder=3)
            ax.plot(phi[::12], ref[::12], "o", color=ORANGE, ms=4, mfc="white", mew=1.1,
                    label="reference simulation", zorder=4)
            ax.set_title(f"{f:g} GHz, {pol}", loc="left", fontsize=9, color=INK)
            ax.set_xticks(range(0, 181, 30))
            ax.grid(True, color="#e4e3df", lw=0.6)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            if i == 1:
                ax.set_xlabel("azimuth [deg]  (0 = nose-on)")
        axes[i, 0].set_ylabel("monostatic RCS [dBsm]")
    h, lab = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, lab, frameon=False, fontsize=8.5, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.005))
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(fname + ".pdf")
    fig.savefig(fname + ".png", dpi=200)
    print("wrote", fname + ".pdf/.png")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        run()
    plot(report())
