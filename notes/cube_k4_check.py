"""
Four-wavelength unit cube check (k = 8π, diagonal ≈ 6.93λ; Section 8 of the paper):
verifies the out-of-the-box unified rules in mfs.cube (s = 3 -> 24 extra modes per unit
length, 5×5 / 5×4 Gauss panels, h_max = 1/16) under 48-element O_h symmetry against
both J = 9 self-convergence and O_h-symmetric graded Bempp at its 8 GB limit (N_u = 41472).

    mfs-safe-run python notes/cube_k4_check.py --part ref
    mfs-safe-run python notes/cube_k4_check.py --part mfs
    mfs-safe-run python notes/cube_k4_check.py --part bem
    python notes/cube_k4_check.py
"""
import gc
import os
import sys
import time

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO)

from mfs import core, cube

NPZ = os.path.join(REPO, "data", "cube_k4_check.npz")
K = 8 * np.pi
D_W, P_W = core.wave()
PTS = 2.0 * core.fibonacci(400)


def load_st():
    return dict(np.load(NPZ)) if os.path.exists(NPZ) else {"pts": PTS}


def run_part(part):
    st = load_st()
    if part in ("all", "ref"):
        cube.solve(2, D_W, P_W, PTS, K)
        gc.collect()
        E_ref, t_ref, info_ref = cube.solve(9, D_W, P_W, PTS, K)
        ref_info = np.array([9, max(info_ref["unknowns"]), info_ref["rows"], t_ref, info_ref["optical"], info_ref["ext"]], float)
        st["ref_E"], st["ref_info"] = E_ref, ref_info
        np.savez(NPZ, **st)
        print(f"k=8pi ref J=9 (m=2, ng=5): unk/cls={int(ref_info[1])}, rows/cls={int(ref_info[2])}, "
              f"t={t_ref:.1f}s, opt={info_ref['optical']:.2e}, ext={info_ref['ext']:.6f}", flush=True)
        if part == "ref":
            return

    E_ref = st["ref_E"]
    if part in ("all", "mfs"):
        cube.solve(2, D_W, P_W, PTS, K)
        prev = {(int(r[0]), int(r[1])): list(r) for r in st.get("mfs", np.zeros((0, 8)))}
        for m, q, j_hi in ((2, 2 ** -0.5, 9), (1, 0.5, 9)):
            for J in range(4, j_hi):
                if (m, J) in prev and f"E_m{m}_J{J}" in st:
                    continue
                gc.collect()
                E, t, info = cube.solve(J, D_W, P_W, PTS, K, q=q)
                err = core.rel(E, E_ref)
                prev[m, J] = [m, J, max(info["unknowns"]), sum(info["unknowns"]), info["rows"], t, err, info["optical"]]
                st[f"E_m{m}_J{J}"] = E
                st["mfs"] = np.array([prev[k] for k in sorted(prev)])
                np.savez(NPZ, **st)
                print(f"k=8pi MFS (m={m}, J={J}): unk/cls={max(info['unknowns'])}, t={t:.2f}s, "
                      f"err_vs_J9={err:.2e}, opt={info['optical']:.2e}", flush=True)
        if part == "mfs":
            return

    if part in ("all", "ng6"):
        cube.solve(2, D_W, P_W, PTS, K)
        orig_guv = cube.graded_uv
        try:
            cube.graded_uv = lambda levels, ng, n_max=0, hmax=cube.H / 4: orig_guv(levels, ng + 1, n_max, hmax)
            gc.collect()
            E8_ng6, t8_ng6, info8_ng6 = cube.solve(8, D_W, P_W, PTS, K)
            shift_ng6 = core.rel(E8_ng6, st["E_m2_J8"])
            err8_ng6 = core.rel(E8_ng6, E_ref)
            st["ng6_J8"] = np.array([t8_ng6, shift_ng6, err8_ng6, info8_ng6["optical"]], float)
            np.savez(NPZ, **st)
            print(f"k=8pi J=8 ng=6 check: t={t8_ng6:.2f}s, shift_vs_ng5={shift_ng6:.2e}, "
                  f"err_vs_J9={err8_ng6:.2e}, opt={info8_ng6['optical']:.2e}", flush=True)
        finally:
            cube.graded_uv = orig_guv
        if part == "ng6":
            return

    if part in ("all", "bem"):
        from mfs import bem
        bempp = bem.bempp
        bem.cube(4, 1.0, K, D_W, P_W, PTS)
        bem_rows = []
        for n_mesh, g, reg, sing in ((40, 1.5, 6, 8), (48, 1.5, 6, 8), (48, 1.5, 8, 10)):
            bempp.GLOBAL_PARAMETERS.quadrature.regular, bempp.GLOBAL_PARAMETERS.quadrature.singular = reg, sing
            gc.collect()
            t0 = time.time()
            E_b, N_b = bem.cube(n_mesh, g, K, D_W, P_W, PTS)
            tb = time.time() - t0
            eb = core.rel(E_b, E_ref)
            bem_rows.append([g, reg, sing, n_mesh, N_b, tb, eb])
            st["bem"] = np.array(bem_rows)
            np.savez(NPZ, **st)
            print(f"k=8pi Bempp (n={n_mesh}, g={g}, {reg}/{sing}): N_u={N_b}, t={tb:.1f}s, diff_vs_J9={eb:.2e}", flush=True)


def show():
    d = np.load(NPZ)
    ri = d["ref_info"]
    print("=== Unit Cube at k = 8*pi (4 wavelengths/side, O_h symmetry) ===")
    print(f"Ref J={int(ri[0])}: unk/cls={int(ri[1])}, rows/cls={int(ri[2])}, t={ri[3]:.1f}s, opt={ri[4]:.2e}, ext={ri[5]:.6f}")
    if "mfs" in d:
        for r in d["mfs"]:
            print(f"  MFS m={int(r[0])} J={int(r[1])}: unk/cls={int(r[2]):4d}, rows={int(r[4]):5d}, t={r[5]:6.2f}s, err={r[6]:.2e}, opt={r[7]:.2e}")
    if "ng6_J8" in d:
        ng6 = d["ng6_J8"]
        print(f"  J=8 ng=6 check: t={ng6[0]:.1f}s, shift_vs_ng5={ng6[1]:.2e}, err_vs_J9={ng6[2]:.2e}, opt={ng6[3]:.2e}")
    if "bem" in d:
        for b in d["bem"]:
            print(f"  Bempp n={int(b[3])} g={b[0]:g} {int(b[1])}/{int(b[2])}: N_u={int(b[4])}, t={b[5]:.1f}s, diff_vs_J9={b[6]:.2e}")


if __name__ == "__main__":
    if "--part" in sys.argv:
        run_part(sys.argv[sys.argv.index("--part") + 1])
    elif "--run" in sys.argv:
        run_part("all")
    else:
        show()
