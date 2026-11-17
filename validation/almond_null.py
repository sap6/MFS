"""How sensitive is the deepest HH null of the almond at 10.25 GHz (35.5 deg,
-66.7 dBsm in the Austin reference simulation) to the size of the almond?

The MFS and the reference simulation differ by 4.3 dB in that null. The
benchmark almond is 9.936 in long; ours is 0.2524 m (9.937 in), 0.01% longer.
Scaling the body by a factor is the same as scaling the frequency, so the
test solves at 10.25 GHz * (1 + delta) for a few delta and prints the HH
cross-section around the null, next to the reference values.

    python3 almond_null.py [n_mid:J]          # the table of sizes
    python3 almond_null.py 6:3 fit           # best-fitting size change and angular shift, on a fine grid of angles
"""
import sys

import numpy as np

import almond_rcs as r
import pec_almond_skeleton as a

FREQ = 10.25
DELTAS = (0.0, 1e-4, -1e-4, 1e-3, -1e-3)
PHI = np.arange(33.0, 38.01, 0.5)


def hh(freq, n_mid, J, s=0.5, s_col=0.3):
    k = 2 * np.pi * freq * 1e9 / r.C0
    y = a.source_rows(n_mid, J, 0.5, s=s, hmax=2 * np.pi / k / 4)
    tail = a.tip_tail(n_mid, J, 0.5, s=s, hmax=2 * np.pi / k / 4)
    d, ph, pv = r.waves(np.deg2rad(PHI))
    amp = 0
    for sy in (+1, -1):
        ee, me = a.tip_exponents(sy, +1, 2, 1)
        f, res, rows, unk = r.solve_class(k, y, sy, +1, d, ph, s_col, tail=tail, exps=ee, mexps=me)
        amp = amp + f
    return 10 * np.log10(4 * np.pi * np.abs(amp) ** 2)


def fit(n_mid, J):
    """Which change of size, or which shift of the look angle, brings the MFS
    curve onto the reference samples around the null?"""
    global PHI
    ref = np.loadtxt("austin_almond/ref_rcs.III.A.s1.fx4.H.txt")
    coarse = np.arange(33.0, 38.01, 0.5)
    refdb = np.interp(coarse, ref[:, 2], ref[:, 3])
    PHI = np.arange(32.5, 38.5001, 0.02)
    idx = [int(np.argmin(np.abs(PHI - p))) for p in coarse]
    print("size change   null at [deg]   depth [dBsm]   rms difference from the reference at its 11 angles [dB]   at 35.5 deg [dB]")
    for dl in (0.0, -1e-4, -2e-4, -3e-4, -4e-4, -5e-4):
        db = hh(FREQ * (1 + dl), n_mid, J)
        i = int(np.argmin(db))
        diff = db[idx] - refdb
        print(f"  {1e2 * dl:+.2f}%      {PHI[i]:8.2f}      {db[i]:8.2f}      {np.sqrt(np.mean(diff ** 2)):8.2f}      {diff[5]:+8.2f}", flush=True)
        if dl == 0.0:
            base = db
    print("angular shift of the nominal MFS curve   rms difference [dB]   at 35.5 deg [dB]")
    for sh in (-0.10, -0.08, -0.06, -0.04, -0.02, 0.0, 0.02):
        diff = np.interp(coarse - sh, PHI, base) - refdb
        print(f"  {sh:+.2f} deg      {np.sqrt(np.mean(diff ** 2)):8.2f}      {diff[5]:+8.2f}")


if __name__ == "__main__":
    n_mid, J = (int(v) for v in sys.argv[1].split(":")) if len(sys.argv) > 1 else (6, 3)
    if "fit" in sys.argv:
        fit(n_mid, J)
        sys.exit()
    ref = np.loadtxt("austin_almond/ref_rcs.III.A.s1.fx4.H.txt")
    refdb = np.interp(PHI, ref[:, 2], ref[:, 3])
    print("azimuth [deg]      " + " ".join(f"{p:7.1f}" for p in PHI))
    print("reference [dBsm]   " + " ".join(f"{v:7.2f}" for v in refdb), flush=True)
    base = None
    for dl in DELTAS:
        db = hh(FREQ * (1 + dl), n_mid, J)
        base = db if base is None else base
        print(f"MFS, size {1e2 * dl:+.2f}%  " + " ".join(f"{v:7.2f}" for v in db)
              + f"   largest change from nominal: {np.abs(db - base).max():.2f} dB", flush=True)
