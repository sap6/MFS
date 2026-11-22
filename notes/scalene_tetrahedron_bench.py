"""
Zero-symmetry convex tetrahedron benchmark (Section 5.6 of the paper):
compares zero-symmetry MFS (ScaleneTetrahedron, m = 1 and m = 2) against
zero-symmetry edge-graded Bempp (bem.efie on scalene_tetrahedron_grid) up to
Bempp's 8 GB RAM ceiling (n = 40, N_u = 14,400).

    python notes/scalene_tetrahedron_bench.py --ref
    python notes/scalene_tetrahedron_bench.py --chunk a|b|c|d
    python notes/scalene_tetrahedron_bench.py --summary
"""
import gc, os, resource, sys, time
import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO)

from mfs import bem, core, tetrahedron

NPZ = os.path.join(REPO, "data", "scalene_tetrahedron.npz")
K = 2 * np.pi
D, P = core.wave()                 # completely non-symmetric with respect to SCALENE_VERTS
PTS = 2.0 * core.fibonacci(400)
VERTS = tetrahedron.SCALENE_VERTS


def warmup(need_mfs=True, b_quads=()):
    if need_mfs:
        tetrahedron.solve_scalene(2, D, P, PTS[:10], K, q=0.5)
    for reg, sing in dict.fromkeys(b_quads):
        bem.bempp.GLOBAL_PARAMETERS.quadrature.regular, bem.bempp.GLOBAL_PARAMETERS.quadrature.singular = reg, sing
        bem.efie(bem.scalene_tetrahedron_grid(VERTS, 4, 1.8), K, D, P, PTS[:10])


def bem_one(n, g, reg, sing, ref):
    bem.bempp.GLOBAL_PARAMETERS.quadrature.regular, bem.bempp.GLOBAL_PARAMETERS.quadrature.singular = int(reg), int(sing)
    grid = bem.scalene_tetrahedron_grid(VERTS, int(n), float(g))
    gc.collect()
    t0 = time.time()
    E, N = bem.efie(grid, K, D, P, PTS)
    return float(n), float(reg), float(N), time.time() - t0, core.rel(E, ref)


def load_state():
    return dict(np.load(NPZ)) if os.path.exists(NPZ) else {}


def run_ref():
    warmup(need_mfs=True)
    gc.collect()
    E, t, info = tetrahedron.solve_scalene(7, D, P, PTS, K, q=2 ** -0.5)
    rss_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    print(f"Scalene ref (m=2, J=7): N={info['unknowns'][0]}, M={info['rows']}, t={t:.2f}s "
          f"(asm={info['t_assemble']:.2f}s, fac={info['t_factor']:.2f}s), opt={info['optical']:.3e}, peak_RSS={rss_mb:.0f} MiB", flush=True)
    st = load_state()
    st["pts"] = PTS
    st["verts"] = VERTS
    st["ref_E"] = E
    st["ref_info"] = np.array([7, info["unknowns"][0], info["rows"], t, info["optical"], info["ext"]])
    np.savez(NPZ, **st)


def run_chunk(cid):
    st = load_state()
    ref = st["ref_E"]
    chunks = {
        "a1": ([(1, J) for J in (2, 3, 4)], [(12, 4, 6), (16, 4, 6), (20, 4, 6), (24, 4, 6)]),
        "a2": ([(1, 5)], [(12, 6, 8), (16, 6, 8)]),
        "a3": ([(1, 6)], [(20, 6, 8)]),
        "b1": ([(2, J) for J in (2, 3, 4)], [(24, 6, 8), (32, 4, 6)]),
        "b2": ([(2, 5)], [(32, 6, 8)]),
        "c1": ([(2, 6)], []),
        "d1": ([], [(40, 4, 6)]),
        "d2": ([], [(40, 6, 8)]),
    }
    m_jobs, b_jobs = chunks[cid]
    warmup(need_mfs=bool(m_jobs), b_quads=[(reg, sing) for _, reg, sing in b_jobs])
    reps = 1 if cid in ("b1", "b2", "c1", "d1", "d2") else 2
    tm = {k: [] for k in m_jobs}
    tb = {k: [] for k in b_jobs}
    rm, rb = {}, {}
    for rep in range(reps):
        for m, J in m_jobs:
            gc.collect()
            E, t, info = tetrahedron.solve_scalene(J, D, P, PTS, K, q=2 ** -0.5 if m == 2 else 0.5)
            tm[m, J].append(t)
            rm[m, J] = (m, J, info["unknowns"][0], info["rows"], core.rel(E, ref), info["optical"])
            print(f"chunk {cid} pass {rep+1} MFS ({m},{J}): N={rm[m,J][2]}, t={t:.2f}s, err={rm[m,J][4]:.3e}, opt={rm[m,J][5]:.3e}", flush=True)
        for n, reg, sing in b_jobs:
            r = bem_one(n, 1.8, reg, sing, ref)
            tb[n, reg, sing].append(r[3])
            rb[n, reg, sing] = r
            print(f"chunk {cid} pass {rep+1} Bempp n={n} ({reg}/{sing}): N={int(r[2])}, t={r[3]:.2f}s, err={r[4]:.3e}", flush=True)

    # Merge into NPZ (keeping min time if run again!)
    m_prev = {(int(r[0]), int(r[1])): list(r) for r in st.get("mfs", np.zeros((0, 7)))}
    for k, vals in rm.items():
        best_t = min(tm[k] + ([m_prev[k][4]] if k in m_prev else []))
        m_prev[k] = [vals[0], vals[1], vals[2], vals[3], best_t, vals[4], vals[5]]
    st["mfs"] = np.array([m_prev[k] for k in sorted(m_prev)]) if m_prev else np.zeros((0, 7))

    b_prev = {(int(r[0]), int(r[1])): list(r) for r in st.get("bem", np.zeros((0, 5)))}
    for (n, reg, sing), vals in rb.items():
        key = (n, reg)
        best_t = min(tb[n, reg, sing] + ([b_prev[key][3]] if key in b_prev else []))
        b_prev[key] = [vals[0], vals[1], vals[2], best_t, vals[4]]
    st["bem"] = np.array([b_prev[k] for k in sorted(b_prev)]) if b_prev else np.zeros((0, 5))
    np.savez(NPZ, **st)


def show_summary():
    st = load_state()
    ri = st["ref_info"]
    print(f"Ref (m=2, J={int(ri[0])}): N_u={int(ri[1])}, M={int(ri[2])}, t={ri[3]:.2f}s, opt={ri[4]:.3e}")
    print("\n=== Zero-Symmetry MFS (m, J, N_u, rows, t, err_vs_J7, opt) ===")
    for r in st["mfs"]:
        print(f"  m={int(r[0])} J={int(r[1])}: N_u={int(r[2]):5d}, M={int(r[3]):5d}, t={r[4]:6.2f}s, err={r[5]:.3e}, opt={r[6]:.3e}")
    c = st["mfs"][st["mfs"][:, 0] == 2]
    p = np.polyfit(np.sqrt(c[:, 2]), np.log(c[:, 5]), 1)
    print(f"Root-exp fit (m=2, J=2..6): err = {np.exp(p[1]):.3g} exp({p[0]:.4f} sqrt(N_u))")
    print("\n=== Zero-Symmetry Edge-Graded Bempp (n, q_reg, N_u, t, err_vs_J7) ===")
    for r in st["bem"]:
        print(f"  n={int(r[0]):2d} ({int(r[1])}/{int(r[1])+2}): N_u={int(r[2]):5d}, t={r[3]:6.2f}s, err={r[4]:.3e}")


if __name__ == "__main__":
    if "--ref" in sys.argv:
        run_ref()
    elif "--chunk" in sys.argv:
        run_chunk(sys.argv[sys.argv.index("--chunk") + 1])
    else:
        show_summary()
