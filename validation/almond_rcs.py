"""Monostatic radar cross-section of the NASA almond against azimuth in the
plane of the almond (z = 0), HH and VV, with the skeleton-source MFS.

Azimuth 0 is nose-on (looking at the tip, wave travelling in -x), 180 is
tail-on. HH: electric field in the plane; VV: electric field along z.

Any incident plane wave is split into its parts of definite parity under the
mirrors y -> -y and z -> -z. Each part is solved with a quarter of the sources
and collocation points, and the QR factorisation of a class serves all look
angles. HH needs the two classes even in z, VV the two classes odd in z.

    python3 almond_rcs.py FREQ_GHz [n_mid:J] [step_deg]
Writes almond_rcs_FREQGHz.npz (angles, sigma_HH, sigma_VV in m^2).
"""
import sys
import time

import numpy as np
from lsq import factor

import pec_almond_skeleton as a
from pec_sphere_convergence import tangent_basis

C0 = 299792458.0


def waves(phi):
    """Propagation direction and the two polarisations for look angles phi (radians)."""
    d = -np.column_stack((np.cos(phi), np.sin(phi), np.zeros_like(phi)))
    ph = np.column_stack((-np.sin(phi), np.cos(phi), np.zeros_like(phi)))
    pv = np.tile(np.array([0.0, 0.0, 1.0]), (len(phi), 1))
    return d, ph, pv


def class_part(x, k, d, p, sig):
    """Part of the plane waves (d, p) with the parity of the class, at points x: (M, n_wave, 3)."""
    out = 0
    for R, s in zip(a.REFL, sig):
        E = p[None, :, :] * np.exp(1j * k * ((x * R) @ d.T))[:, :, None]
        out = out + s * E * R[None, None, :]
    return out / 4


def solve_class(k, y, sy, sz, d, p, s_col, ratio=2.5, tail=None, exps=a.TIP_EXPS, mexps=()):
    """Coefficients for all waves (d, p) in the class (sy, sz), and its backscatter amplitudes p . F(-d)."""
    comps, sig = a.parity_class(sy, sz)
    x, w, nrm = a.collocation(a.ETA_F, s_col, a.col_hmax(k))
    if 2 * len(x) < ratio * len(comps) * len(y):
        s_col *= np.sqrt(2 * len(x) / (ratio * len(comps) * len(y)))
        x, w, nrm = a.collocation(a.ETA_F, s_col, a.col_hmax(k))
    T = tangent_basis(nrm)
    sw = np.sqrt(w / w.sum())
    A = a.matrix(x, y, k, T, comps, sig)
    if tail is not None:
        A = np.concatenate((A, a.tail_matrix(x, tail, k, T, comps, sig, exps, mexps)), axis=1)
    A *= np.repeat(sw, 2)[:, None]
    rhs = -np.einsum("mab,mnb->man", T, class_part(x, k, d, p, sig)) * sw[:, None, None]
    rhs = rhs.reshape(-1, len(d))
    c, res = factor(np.ascontiguousarray(A)).solve(rhs, residual=True)     # semi-normal equations or QR (lsq.py)
    res = np.array([res, np.linalg.norm(rhs, axis=0)])                                       # absolute, per wave
    rfar = 1.0e5
    Ffar = a.matrix(-rfar * d, y, k, None, comps, sig)                                       # (3 n_wave, unknowns)
    if tail is not None:
        Ffar = np.concatenate((Ffar, a.tail_matrix(-rfar * d, tail, k, None, comps, sig, exps, mexps)), axis=1)
    Ffar = Ffar * rfar * np.exp(-1j * k * rfar)
    F = np.einsum("nau,un->na", Ffar.reshape(len(d), 3, -1), c)                               # field of wave n in its own direction
    return np.einsum("na,na->n", p, F), res, len(x) * 2, A.shape[1]


def rcs(freq, n_mid=6, J=3, step=2.0, s_col=0.3, tip=True, s=0.5, tip_exps="cone", n_e=2, n_m=1):
    """s: spacing of the sources along a row, as a fraction of the distance to
    the surface. It is the main convergence parameter for the waves that
    drive current along the axis (0.7: 3%, 0.5: 0.2%, 0.4: 0.03% residual)."""
    k = 2 * np.pi * freq * 1e9 / C0
    y = a.source_rows(n_mid, J, 0.5, s=s, hmax=2 * np.pi / k / 4)
    tail = a.tip_tail(n_mid, J, 0.5, s=s, hmax=2 * np.pi / k / 4) if tip else None
    phi = np.deg2rad(np.arange(0.0, 180.0 + 1e-9, step))
    d, ph, pv = waves(phi)
    out = {"phi_deg": np.rad2deg(phi), "freq": freq, "N": len(y)}
    for name, p, sz in (("HH", ph, +1), ("VV", pv, -1)):
        amp, res2 = 0, 0
        for sy in (+1, -1):
            t0 = time.time()
            ee, me = a.tip_exponents(sy, sz, n_e, n_m) if tip_exps == "cone" else (a.TIP_EXPS, ())
            f, res, rows, unk = solve_class(k, y, sy, sz, d, p, s_col, tail=tail, exps=ee, mexps=me)
            amp, res2 = amp + f, res2 + res ** 2
            print(f"  {name} class (sy={sy:+d}, sz={sz:+d}): unknowns={unk} rows={rows} largest residual / largest "
                  f"right-hand side over the angles = {res[0].max() / res[1].max():.1e} ({time.time() - t0:.0f}s)", flush=True)
        out["sigma_" + name] = 4 * np.pi * np.abs(amp) ** 2
        out["res_" + name] = np.sqrt(res2[0] / res2[1])       # boundary residual per look angle, both classes
    return out


if __name__ == "__main__":
    freq = float(sys.argv[1])
    n_mid, J = (int(v) for v in sys.argv[2].split(":")) if len(sys.argv) > 2 else (6, 3)
    step = float(sys.argv[3]) if len(sys.argv) > 3 else 2.0
    tip = "--notip" not in sys.argv
    s_src = ([float(v[4:]) for v in sys.argv if v.startswith("--s=")] or [0.5])[0]
    s_col = ([float(v[7:]) for v in sys.argv if v.startswith("--scol=")] or [0.3])[0]
    print(f"almond, {freq} GHz, {a.L * freq * 1e9 / C0:.2f} wavelengths long, sources n_mid={n_mid} J={J}", flush=True)
    out = rcs(freq, n_mid, J, step, s_col=s_col, tip=tip, s=s_src)
    np.savez(f"almond_rcs_{freq}GHz_{n_mid}_{J}_s{s_src}{'' if tip else '_notip'}{'' if s_col == 0.3 else f'_c{s_col}'}.npz", **out)
    try:                                           # Bempp reference, if one was stored for this frequency
        import json
        ref = json.load(open("bem_almond_rcs_ref.json"))
        r = ref[str(freq)]
        idx = [int(np.argmin(np.abs(out["phi_deg"] - p))) for p in ref["phi_deg"]]
        if np.allclose(out["phi_deg"][idx], ref["phi_deg"]):
            for pol in ("HH", "VV"):
                db = 10 * np.log10(out["sigma_" + pol][idx] / np.array(r[pol]))
                print(f"{pol} against Bempp ({r['ndof']} dofs) at {len(idx)} angles: rms {np.sqrt(np.mean(db ** 2)):.3f} dB, "
                      f"max {np.abs(db).max():.3f} dB at {ref['phi_deg'][int(np.argmax(np.abs(db)))]} deg")
    except (OSError, KeyError):
        pass
    lam2 = (C0 / (freq * 1e9)) ** 2
    print("azimuth[deg]   HH[dBsm]   VV[dBsm]   (HH, VV in dB relative to lambda^2)   boundary residual HH, VV")
    for i in range(0, len(out["phi_deg"]), max(1, int(round(30 / step)))):
        hh, vv = out["sigma_HH"][i], out["sigma_VV"][i]
        print(f"{out['phi_deg'][i]:8.1f} {10 * np.log10(hh):11.2f} {10 * np.log10(vv):10.2f}   "
              f"({10 * np.log10(hh / lam2):.2f}, {10 * np.log10(vv / lam2):.2f})   "
              f"{out['res_HH'][i]:.1e} {out['res_VV'][i]:.1e}")
