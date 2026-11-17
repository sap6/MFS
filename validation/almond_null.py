"""How sensitive is the deepest HH null of the almond at 10.25 GHz (35.5 deg,
-66.7 dBsm in the Austin reference simulation) to the size of the almond?

The MFS and the reference simulation differ by 4.3 dB in that null. The
benchmark almond is 9.936 in long; ours is 0.2524 m (9.937 in), 0.01% longer.
Scaling the body by a factor is the same as scaling the frequency, so the
test solves at 10.25 GHz * (1 + delta) for a few delta and prints the HH
cross-section around the null, next to the reference values.

    python3 almond_null.py [n_mid:J]
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


if __name__ == "__main__":
    n_mid, J = (int(v) for v in sys.argv[1].split(":")) if len(sys.argv) > 1 else (6, 3)
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
