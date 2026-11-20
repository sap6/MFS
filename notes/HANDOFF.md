# Handoff (20 November 2026)

Read this first, then `CLAUDE.md` and the top of `notes/NOTES.md`.

## 1. State of the repository

- Branch `ccr-116c520c-fts35f` of `sap6/MFS`, everything committed and
  pushed. No pull request; not merged into `main` (merge into `main` before
  sharing Colab links). Nothing was ever pushed to `sap6/mfs-paper` (it syncs
  to Overleaf): never push there without asking the author.
- `python -m pytest tests` passes (`12` tests with `MFS_TEST_BEMPP=1`), CI on
  GitHub is green, and `python figures/tables.py` reproduces `paper/tables/`.
- Layout: `mfs/` (the package), `figures/` (one script per figure, replots
  from `data/`, `--run` recomputes), `data/`, `notebooks/` (Colab),
  `tests/`, `paper/`, `thesis/` (the author's MSci notebooks), `notes/`.
  Setup commands are in `CLAUDE.md`.

## 2. What changed in the last two sessions

All of it is recorded, with numbers, in `notes/NOTES.md` (part three, and the
experiment log). The essentials:

- **Code for readers.** The package replaced 7300 lines of research scripts.
  It reproduces every stored result: see `notes/code_for_readers_18nov.md` and
  `notes/code_size_analysis.md`.
- **Libraries, and what speed costs.**
  - Library calls replace hand-written code wherever they are as fast.
  - Every speed-up that remains was measured against its simple
    alternative: `notes/simplicity_pass_18nov.md`,
    `notes/speed_vs_simplicity_18nov.md`. Do not repeat these measurements.
- **A faster MFS on the cube, fairly** (`notes/faster_mfs_18nov.md`):
  1. The boundary panels now halve towards the edges until they are shorter
     than d_min, as Section 5 of the paper already said. The code had graded
     two levels further.
  2. The line kernel no longer computes four of the eight mirror images. A
     line is its own mirror image across the plane through its middle, so
     those images repeat the other four. The almond kernel is the same, since
     its sources lie on z = 0.
  3. The almond's optical check is no longer timed, as was already the case
     for Bempp.
  4. The comparison uses both MFS layouts (two and one depths per octave), as
     Bempp's curve uses the best of its scan.

  Measured crossovers with Bempp on the 2-core Intel Xeon Platinum 8481C (8 GB RAM) machine (the MFS is faster below):

  | case | best MFS layout (1 or 2 depths/octave) | two depths/octave alone |
  |---|---|---|
  | standard wave (`k1`) | 5.7e-4 | 2.4e-4 |
  | general wave (`k1o`) | 5.2e-4 | 2.0e-4 |
  | two wavelengths (`k2`) | 7.6e-4 | 4.8e-4 |
  | three wavelengths (`k3`) | 1.3e-3 | 1.0e-3 |

## 3. Completed on 19--20 November 2026

- **3.1 Unified cube rules & refreshed benchmarks (commit `80928a0` + 20 Nov pass):**
  - Unified `mfs/cube.py` boundary panel cap `hmax = min(1/8, pi/(2k))` (`lambda/4`), panel order `2 + max(m, 1 + ceil(s/2))`, physical line-length mode rule `round(2*s*(s+1)*2*(H-d))`, and line quadrature `12 + ceil(8*s), 8 + ceil(2*s)` with zero per-case overrides; F-contiguous per-class slices `np.moveaxis(np.zeros((rows, nmax, g), order="F"), 2, 0)` avoid implicit `zherk` and `QR` copies.
  - Eliminated 1.4 GB of matrix/tensor copies in `mfs/almond.py` (`matrix`, `tail_matrix` via chunked BLAS `zgemm`, vectorised `collocation`), speeding up `1.19 GHz` nose-on MFS by `1.5x--3.9x` (`3:0..8:5` in `0.67, 1.48, 2.38, 5.10, 12.94 s`, `10^{-15}` field agreement) so MFS is now faster than Bempp across the entire accuracy range on sphere, torus, and `1.19 GHz` almond.
  - Fixed `torus.residual` to default to `(sol.R, sol.a)` and scaled `CSNE` Gram matrix in-place.
- **3.2 Figures, tables, manuscript, notebooks & tests (20 Nov):**
  - Raised PNG export to `300 dpi`; replaced aliased `6°` dot subsampling in `fig_almond_benchmark` with full `0.5°` dashed reference curve; fixed `fig_almond_layout` `1.85x` aspect distortion and `fig_cube_layout` text/smudge overlaps; added `$2^{-9}d_{\min}$` marker in `fig_edge_residual`; unified axes and legends across `fig_cube_cases`, `fig_cost_all`, `fig_sphere`, `fig_torus_sources`, and `fig_sphere_nearfield`; expanded Table 4 (`summary.tex`) to show both `Total` and `Solved` unknowns.
  - Refined manuscript exposition, referee preemptions (open-skeleton absence of spurious Maxwell resonances, unit-cube EFIE cavity parity orthogonality, CSNE iterative-refinement sensitivity, verified `10^7` column-scaled condition bound, `0.20%`/`0.033%` almond `(1,1)` residuals, `-66.7` vs `-71.1` dBsm `10.25 GHz` null flank, DOIs and pruned `refs.bib`), compiling `main.pdf` (46 pp) and `main_3p.pdf` (27 pp) with `0` warnings and `0` Overfull/Underfull boxes.
  - Eliminated three `>7.5 GiB` RAM / swap-thrash hazards on 8 GB VMs: capped `lines._TREES` at 32 entries (saves `0.75--1.55 GiB`), made `bem.fold` F-ordered and factored one orbit class at a time while freeing `sym.A` before OpenCL potential evaluation (cuts `bem.cube` `k1o` `n=48` peak RSS by `~2.4 GiB`), and included `CSNE` Cholesky factor `nmax^2` + cgroup limits in `cube.solve` batching (`k1o J=11` `5.53 -> 2.77 GiB`).
  - Upgraded `README.md`, all 4 Colab notebooks (`1_sphere.ipynb`..`4_almond.ipynb`, re-executed cleanly), and `tests/test_mfs.py` (`12` tests with `MFS_TEST_BEMPP=1`).

## 4. Decisions only the author can make (asked, not yet answered)

- Merge into `main` without squashing. Tag the research scripts:
  `git tag research-scripts 9f5a721` existed only in the old container.
- The LICENSE holder ("Sap"); a Zenodo release and DOI; whether `CLAUDE.md`
  and `notes/` stay public.
- The JCP declarations (`notes/jcp_checklist.md`).
- The presentation ideas in `notes/humanising_pass_18nov.md` (not done).

## 5. Facts that are easy to get wrong

- **Fairness rule (the author's):** every symmetry or speed-up goes into both
  codes or into neither. Bempp gets its scan over mesh, grading and
  quadrature, the symmetry, the compiled fold and LU.
- **Bempp:** import `mfs.bem` before `bempp` (it sets `HWLOC_SYNTHETIC` for
  the OpenCL memory). Its default quadrature in bempp-cl 0.4.2 is 4/4.
- **What failed for speed:**
  - interior panels 0.15 fail at three wavelengths and at J = 12;
  - grading four levels short of the paper's rule gives errors near 1 (three
    is fine but close to the cliff);
  - three nodes per panel with two depths per octave give garbage;
  - fewer modes do not pay;
  - SciPy's `lstsq` is 2 to 4 times slower than our CSNE/QR;
  - NumPy kernels are 20 to 50 times slower than the compiled ones.
- **QR fallback:** CSNE falls back to in-place QR only when the Cholesky diagonal of $R$ spans more than four orders of magnitude (`d.min() <= 1e-4 * d.max()`); at `1--3` wavelengths with the unified mode rule, CSNE succeeds at every level `J = 4..12` (the fallback triggers e.g. at `4` wavelengths, `J = 6`).
- **The author's voice:** `notes/style_profile.md` (from `thesis/`). The
  title is settled. The abstract has at most 250 words.
