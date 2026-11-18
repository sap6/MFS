"""Numbers quoted in the comparison of the MFS with Bempp, from the timing
data: for each case the time of each method at given errors (interpolated
in log(time) against log(error) between neighbouring runs), their ratio,
and the error at which the two cost curves cross.

    python3 crossovers.py
"""
import numpy as np


def curve(rows):
    """Lower envelope sorted by decreasing error: arrays (error, time)."""
    out, best = [], np.inf
    for t, e in sorted(rows):
        if e < best:
            out.append((e, t))
            best = e
    return np.array(out)


def time_at(c, err):
    """Time needed for the error err (None outside the range of the curve)."""
    e, t = c[:, 0], c[:, 1]
    if err > e[0] or err < e[-1]:
        return None
    return float(np.exp(np.interp(-np.log(err), -np.log(e), np.log(t))))


def crossing(a, b):
    lo, hi = max(a[-1, 0], b[-1, 0]), min(a[0, 0], b[0, 0])
    errs = np.exp(np.linspace(np.log(hi), np.log(lo), 4000))[1:-1]
    r = np.array([time_at(a, e) / time_at(b, e) for e in errs])
    i = np.nonzero(np.diff(np.sign(np.log(r))))[0]
    return [(errs[k], time_at(a, errs[k])) for k in i]


def report(name, mfs, bem, errs):
    print(f"{name}: MFS smallest error {mfs[-1, 0]:.2e} in {mfs[-1, 1]:.1f} s; Bempp smallest error {bem[-1, 0]:.2e} in {bem[-1, 1]:.1f} s")
    for e in errs:
        tm, tb = time_at(mfs, e), time_at(bem, e)
        if tm and tb:
            print(f"   error {e:.1e}: MFS {tm:6.1f} s, Bempp {tb:6.1f} s, Bempp faster by {tm / tb:5.2f}" if tm > tb else
                  f"   error {e:.1e}: MFS {tm:6.1f} s, Bempp {tb:6.1f} s, MFS faster by {tb / tm:5.2f}")
        else:
            print(f"   error {e:.1e}: MFS {'-' if tm is None else f'{tm:.1f} s'}, Bempp {'-' if tb is None else f'{tb:.1f} s'}")
    for e, t in crossing(mfs, bem):
        print(f"   the curves cross at {e:.2e}, {t:.0f} s")


if __name__ == "__main__":
    d = np.load("fig_cost_cube.npz")
    ref = d["ref_E"]
    rel = lambda E: np.linalg.norm(E - ref) / np.linalg.norm(ref)
    jref = int(d["ref_info"][0])
    env = np.load("bem_grading_scan_k1_final.npz")["rows"]
    bem = curve([(r[5], r[6]) for r in env])
    for key, name in (("dense", "cube, standard wave"), ("tail", "cube, standard wave, one depth per octave")):
        mfs = curve([(r[1], rel(d[f"{key}_E{int(r[0])}"])) for r in d[key] if int(r[0]) >= 4 and not (key == "dense" and int(r[0]) == jref)])
        report(name, mfs, bem, (1e-3, 3.5e-4, 3e-4, 2e-4, 1.4e-4, 1.15e-4, 8.2e-5))
    for case, name, errs in (("k1o", "cube, general wave", (1e-3, 3e-4, 1.5e-4, 1.1e-4)), ("k2", "cube, two wavelengths", (1e-3, 4e-4, 2.5e-4, 2e-4)),
                             ("k3", "cube, three wavelengths", (2e-3, 1e-3, 6e-4, 5.2e-4))):
        c = np.load(f"cube_cases_{case}.npz")
        report(name, curve([(r[5], r[6]) for r in c["mfs"]]), curve([(r[5], r[6]) for r in c["env"]]), errs)
    for shape in ("sphere", "torus"):
        c = np.load(f"fig_cost_{shape}.npz")
        report(shape, curve([(r[1], r[2]) for r in c["mfs"]]), curve([(r[1], r[2]) for r in c["bem"]]), (1e-1, 1e-2, 4e-3, 1e-3, 5e-4))
    c = np.load("fig_cost_almond.npz")
    ref = c[f"mfs_E{len(c['mfs']) - 1}"]
    rel = lambda E: np.linalg.norm(E - ref) / np.linalg.norm(ref)
    report("almond", curve([(c["mfs"][i, 2], rel(c[f"mfs_E{i}"])) for i in range(len(c["mfs"]) - 1)]),
           curve([(c["bem"][i, 1], rel(c[f"bem_E{i}"])) for i in range(len(c["bem"]))]), (1e-2, 6e-3, 3e-3, 1e-3))
