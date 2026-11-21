# Handoff (21 November 2026)

Read this first, then `CLAUDE.md` and the top of `notes/NOTES.md`.

## 1. State of the repository

- Branch `ccr-116c520c-fts35f` of `sap6/MFS`, everything committed. No pull
  request; not merged into `main` (merge into `main` before sharing Colab
  links). Nothing was ever pushed to `sap6/mfs-paper` (it syncs to Overleaf):
  never push there without asking the author.
- `python -m pytest tests` passes (`16` tests with `MFS_TEST_BEMPP=1`), CI on
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
- **A faster MFS across all five shapes, fairly** (`notes/faster_mfs_18nov.md` + 21 Nov campaign):
  1. The boundary panels halve towards the edges until they are shorter
     than `d_min`, as Section 5 of the paper states.
  2. The line and almond kernels skip midpoint/mid-plane duplicate mirror images,
     accumulate across active symmetry images in registers before writing to F-ordered
     class matrices, and evaluate near and exact `r -> inf` far fields directly from
     contracted mode/multipole moments without re-forming `A`.
  3. The sphere uses fast Numba tangential block assembly + `CSNE` with automatic
     `QR` fallback (`cond > 10^8`), the torus vectorises block-circulant assembly and
     batch `CSNE`/SVD solves across azimuthal Fourier modes, and the almond uses a single
     fused symmetric dipole kernel (`almond.rcs` sweeps 361 azimuths in 3.0 min/frequency).
  4. Both MFS and Bempp were re-timed across **two alternating passes** on the same idle
     2-core Intel Xeon @ 2.3 GHz (8 GB RAM) VM, keeping the smaller time of each configuration.

  Measured crossovers with symmetric graded Bempp on the 2-core Intel Xeon @ 2.3 GHz (8 GB RAM) VM (the MFS is faster below):

  | case | best MFS layout (1 or 2 depths/octave) | two depths/octave alone |
  |---|---|---|
  | cube, standard wave (`k1`) | 6.0e-4 | 2.6e-4 |
  | cube, general wave (`k1o`) | 5.6e-4 | 2.6e-4 |
  | cube, two wavelengths (`k2`) | 8.8e-4 | 5.0e-4 |
  | cube, three wavelengths (`k3`) | 1.3e-3 | 6.2e-4 |
  | regular tetrahedron (`k = 2*pi`) | 8.5e-4 | 8.5e-4 |

## 3. Completed on 19--21 November 2026

- **3.1 Unified cube rules & refreshed benchmarks (commit `80928a0` + 20 Nov pass):**
  - Unified `mfs/cube.py` boundary panel cap `hmax = min(1/8, pi/(2k))` (`lambda/4`), panel order `2 + max(m, 1 + ceil(s/2))`, physical line-length mode rule `round(2*s*(s+1)*2*(H-d))`, and line quadrature `12 + ceil(8*s), 8 + ceil(2*s)` with zero per-case overrides; F-contiguous per-class slices `np.moveaxis(np.zeros((rows, nmax, g), order="F"), 2, 0)` avoid implicit `zherk` and `QR` copies.
  - Eliminated 1.4 GB of matrix/tensor copies in `mfs/almond.py`, speeding up `1.19 GHz` nose-on MFS and fixing `torus.residual` to default to `(sol.R, sol.a)`.
- **3.2 Anti-stall & referee audit pass (20 Nov, commit `8a17127`):**
  - Capped `lines._TREES` at 32 entries, made `bem.fold` F-ordered and single-LU-at-a-time, and added `CSNE` `nmax^2` + cgroup RAM awareness to `cube.solve` and `tetrahedron.solve`.
- **3.3 Fifth benchmark geometry + 20-hour speed, accuracy & figure campaign (20--21 Nov):**
  - **Regular tetrahedron (`mfs/tetrahedron.py`, `figures/tetrahedron.py`, `data/tetrahedron_convergence.npz`, Section 5.6):** Eliminated the 1-depth/octave (`m = 1`) error uptick (`3.77e-2 -> 1.74e-2 -> 7.89e-3 -> 3.70e-3 -> 2.12e-3 -> 1.84e-3`, strictly monotonic across `J = 3..8` via `max(4, ng)` kite collocation, `0.56`-scaled `n_base = 12` inner base meeting `d_0 = 0.14`, and `1.2 d_min` corner rays), cut `m = 2` optical mismatch 13x (`2.25e-5` at `J = 10`, fit `2.05 exp(-0.0695 sqrt(N_u))` within `12%`, `9.8e-5` at `J = 9` in `60 s`, oblique wave `2.4e-2 -> 6.4e-4`), and improved symmetric graded Bempp confirmation 2x (`2.21e-4` vs `J = 10`, `2.03e-4` vs `J = 9` at `N_u = 17424` in `257 s` using `8/10` quadrature on `n = 16, 24, 32, 44`). Also verified untuned multi-wavelength convergence on the regular tetrahedron for referee readiness: at `k = 4*pi` (`J = 3..9` vs `J = 10`: `2.8e-2 -> 1.7e-4` in `79 s`, `opt = 2.6e-5`, confirmed by symmetric graded Bempp `n = 44` to `3.2e-4` in `257 s`) and at `k = 6*pi` (`J = 3..7` vs `J = 8`: `2.5e-2 -> 1.1e-3`, `opt = 7.2e-5`).
  - **Fair MFS speedups & complete 2-pass alternating re-timing across all 5 shapes (`commits 6315a9b..27c62b9`):** Accelerated `sphere` (`0.56 s` at `N = 400`, `10 s` at `N = 1200` for `3.5e-14`), `torus` (`0.31 s` at `21120` unknowns for `7.3e-12`), `almond` (`8.2 s` at `6:3` `1.19 GHz` for `5.6e-5`, `8.4 s` at `3.5 GHz` for `4.6e-5`, `179 s = 3.0 min/freq` for 361-azimuth RCS up to `10.25 GHz`), `cube` (`k1`, `k1o`, `k2`, `k3`, `cube_convergence`, `cube_bempp_uniform` — e.g. `k1 J = 9` `5.5e-5` in `76 s` vs Bempp `6.2e-5` in `636 s`), and `tetrahedron` (`8.5e-4` crossover at `24 s`, `3.2e-4` in `39 s`). Exact `r -> inf` far-field evaluation (`Cube.far`, `Tetrahedron.far`, `almond.far_field`) replaces finite-radius spheres.
  - **Super-clear placement & benchmark figures + 24-page JCP audit (`commits 4a3d309..2eb3442`):** Redesigned all 13 figures (`fig_cube_layout`, `fig_almond_layout`, `fig_tetrahedron` 4-panel `(a)-(d)`, `fig_wedge`, `fig_edge_residual`, `fig_sphere`, `fig_torus_sources`, `fig_graphical_abstract`, and convergence/cost figures) so geometric/lightning clustering, stretched `t(u)` nodes, singular tails, focal sets, and dyadic boundary panels are unmistakable at journal print width (`[final,3p,times]`, visually inspected across all 24 rendered pages). Qualified sphere interior TE resonances (`j_n(k R_s) != 0`), verified `cond(A_scaled) < 10^7` across all shapes (`3e3..2.7e6`), updated and re-executed all 4 notebooks (`notebooks/1..4`), and expanded `tests/test_mfs.py` to `16` passing checks (`main_3p.pdf`: **24 pages**, `main.pdf`: **40 pages**, `0` warnings/bad boxes).

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
