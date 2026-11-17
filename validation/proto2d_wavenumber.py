"""2D prototype for the wavenumber question of NOTES.md item 57: which parts of
the corner layout of the unit square must follow the wavelength? Free point
sources on the corner bisectors (ratio 2^(-1/2), down to 1e-7) plus base
sources on an inner square, TE polarisation, side = S wavelengths.

    python3 proto2d_wavenumber.py

Variants: A the layout of one wavelength (10 base sources per side at depth
0.2, collocation panels 0.1); B collocation panels 0.1 / S; C also 10 S base
sources per side at depth 0.2; D also the depth 0.2 / S. The reference is D
with 1.4 times the base sources, finer panels and six more depths.
"""
import numpy as np
from scipy.linalg import lstsq

import proto2d_square as p2

Q = 2 ** -0.5


def solve(S, n_base, depth0, hmax, d_min=1e-7, extra=0, p=4, cond=1e-13):
    p2.K = 2 * np.pi * S
    d0 = depth0 * np.sqrt(2)
    J = int(np.ceil(np.log(d_min / d0) / np.log(Q))) + extra
    y = p2.sources(J, Q, (0.0,), n_base, depth0)
    lv = max(3, int(np.ceil(np.log(d0 * Q ** J / p2.H) / np.log(0.5))) + 2)
    x, w, nrm, _ = p2.boundary(lv, p, 0.5, hmax)
    rw = np.sqrt(w) * p2.corner_dist(x) ** 0.5
    A, b = p2.kernels(x, y, nrm, True) * rw[:, None], -p2.incident(x, nrm) * rw
    sc = np.linalg.norm(A, axis=0)
    c = lstsq(A / sc, b, lapack_driver="gelsy", cond=cond)[0]
    return p2.far_field(y, c / sc, True), A.shape[1], A.shape[0], np.linalg.norm(c)


if __name__ == "__main__":
    print(" S   A: layout of S = 1     B: panels / S        C: + base x S        D: + depth / S      (error of the far field, unknowns, |c|)")
    for S in (1, 2, 3, 4, 6):
        ref = solve(S, int(14 * S), 0.2 / S, 0.07 / S, extra=6)[0]
        line = f"{S:2d} "
        for n_base, depth0, hmax in ((10, 0.2, 0.1), (10, 0.2, 0.1 / S), (10 * S, 0.2, 0.1 / S), (10 * S, 0.2 / S, 0.1 / S)):
            ff, unk, rows, cn = solve(S, n_base, depth0, hmax)
            line += f"  {np.linalg.norm(ff - ref) / np.linalg.norm(ref):.1e} ({unk:4d}, {cn:7.1e})"
        print(line, flush=True)
