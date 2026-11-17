"""2D prototype (unit square, proto2d_square_tail.py): does the stall of the
3D cube family with one source depth per octave appear in 2D too? Far-field
error against a free solve with 40 depths, for depth ratios 1/2 and 2^(-1/2)
at the same deepest depths, with one tail exponent per source type (as in
3D) and with the full set.

    python3 proto2d_depth_ratio.py TE|TM
"""
import sys

import numpy as np

import proto2d_square as p2
import proto2d_square_tail as t2

if __name__ == "__main__":
    pol = sys.argv[1]
    yr, cr, _ = p2.solve(pol, J=40, q=0.6, beta=0.0 if pol == "TM" else 0.5, dip=True)   # 70 depths (down to 1e-16) break down with numpy 2.5 / scipy 1.18; 30, 40 and 50 agree to 3e-7
    ref = p2.far_field(yr, cr, True)
    S = 2 ** -0.5
    print(f"{pol}: far-field error against the free solve with 40 depths")
    print(" level  deepest depth   q=1/2, 1+1 exp.   q=1/2, 4+2 exp.   q=2^-1/2, 1+1 exp.   q=2^-1/2, 4+2 exp.   (unknowns)")
    for J in (2, 4, 6, 8, 10, 12, 14):
        line, unk = f"{J:5d}   {t2.D0 * 0.5 ** J:.1e}  ", []
        for q, JJ in ((0.5, J), (S, 2 * J)):
            for nm, nd in ((1, 1), (4, 2)):
                ff, info = t2.solve(pol, JJ, q=q, exps_m=t2.EXP_M[:nm], exps_d=t2.EXP_D[:nd], qt=0.5)
                line += f"     {np.linalg.norm(ff - ref) / np.linalg.norm(ref):.1e}     "
                unk.append(info["unk"])
        print(line + f"  {unk}", flush=True)
