# Handoff (22 November 2026)

Read this first, then `CLAUDE.md` and the top of `notes/NOTES.md`.

## 1. State of the repository

- Branch `ccr-116c520c-fts35f` of `sap6/MFS`, everything committed. No pull
  request; not merged into `main` (merge into `main` before sharing Colab
  links). Nothing was ever pushed to `sap6/mfs-paper` (it syncs to Overleaf):
  never push there without asking the author.
- `python -m pytest tests` passes (`19` tests with `MFS_TEST_BEMPP=1`), all
  `figures/*.py` scripts and `git diff --exit-code paper/tables` pass locally,
  and `python figures/tables.py` reproduces `paper/tables/`.
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
- **A faster MFS across all five shapes, fairly** (`notes/faster_mfs_18nov.md` + 21--22 Nov campaigns):
  1. The boundary panels halve towards the edges until they are shorter
     than `d_min`, as Section 5 of the paper states, with universal 2-branch
     anisotropic Gauss--Legendre panels (`ng x ng` for `|p-q| <= 1`,
     `ng x max(3, ng-1)` for `|p-q| >= 2`, cutting rows by `20%` with zero loss
     of conditioning).
  2. `lines.py` vectorises `build_tree` moment panels (`~30x` faster), centers
     `x - y(t) = p_perp - (t - p_t) t_hat` orthogonally at each target's
     projection (`rel |A_new - A_orig| = 3e-16`), integrates 6 scalar radial
     weights per panel, flushes L2-tiled (`BM = 64`) columnwise target blocks into
     F-ordered class matrices (`2.3x` faster line assembly), contracts root-panel
     moments for distant exterior targets (`~100x` faster `field()`), exact
     `r -> inf` far fields, and `@njit` exact transpose panel/moment adjoints
     (`eval_points_adj`, `eval_line_adj`).
  3. Under the full `48`-element octahedral group `O_h` (matched in both `mfs/cube.py`
     and `mfs/bem.py`), Green's kernels are evaluated once on the `y <= z` triangle of
     face `x = 1/2` (`1/48` of the surface), cyclic `(x, y, z)` orbits map onto at most
     four representative classes, and invariant/mixed classes split into three `120°`-rotation
     (`C_3`) or two `y <-> z` reflection blocks.
  4. Keeping two free vertex dipoles per octave (`mv = max(2, m)`) when `q = 1/2`
     restores `CSNE` conditioning across all levels, taking `1 depth/octave` down to
     `7.12e-5` in `3.8 s` (`5.1 s` in the alternating pass) on the cube (`J = 9`) and
     `1.45e-3` in `4.40 s` on the regular tetrahedron (`J = 7`), matching the
     Trefethen--Weideman strip-analyticity floor `exp(-2*pi*theta_e / |ln q|)`.
  5. Both MFS and Bempp were re-timed across **two alternating passes** on the same idle
     2-core Intel Xeon @ 2.3 GHz (8 GB RAM) VM, keeping the smaller time of each configuration.

  Measured crossovers with symmetric graded Bempp on the 2-core Intel Xeon @ 2.3 GHz (8 GB RAM) VM (the MFS is faster below):

  | case | best MFS layout (1 or 2 depths/octave) | two depths/octave alone |
  |---|---|---|
  | cube, standard wave (`k1`) | faster throughout (`1.5 s` vs `2.1 s` at `10^-3`, `15 s` vs `148 s` at `6.3e-5`) | 4.2e-4 |
  | cube, general wave (`k1o`) | faster throughout (`1.9 s` vs `2.5 s` at `10^-3`, `25 s` vs `165 s` at `7.6e-5`) | 4.9e-4 |
  | cube, two wavelengths (`k2`) | faster throughout (`2.2 s` vs `5.3 s` at `10^-3`, `46 s` vs `148 s` at `1.3e-4`) | 8.2e-4 |
  | cube, three wavelengths (`k3`) | faster throughout (`4.8 s` vs `12 s` at `10^-3`, `11 s` vs `147 s` at `3.1e-4`) | 5.6e-4 |
  | regular tetrahedron (`k = 2*pi`) | faster throughout (`4.4 s` vs `6.8 s` at `1.5e-3`, `36 s` vs `267 s` at `10^-4`) | 1.8e-3 |

## 3. Completed on 19--22 November 2026

- **3.1 Unified cube rules & refreshed benchmarks (commit `80928a0` + 20 Nov pass):**
  - Unified `mfs/cube.py` boundary panel cap `hmax = min(1/8, pi/(2k))` (`lambda/4`), panel order `2 + max(m, 1 + ceil(s/2))`, physical line-length mode rule `round(2*s*(s+1)*2*(H-d))`, and line quadrature `12 + ceil(8*s), 8 + ceil(2*s)` with zero per-case overrides; F-contiguous per-class slices `np.moveaxis(np.zeros((rows, nmax, g), order="F"), 2, 0)` avoid implicit `zherk` and `QR` copies.
  - Eliminated 1.4 GB of matrix/tensor copies in `mfs/almond.py`, speeding up `1.19 GHz` nose-on MFS and fixing `torus.residual` to default to `(sol.R, sol.a)`.
- **3.2 Anti-stall & referee audit pass (20 Nov, commit `8a17127`):**
  - Capped `lines._TREES` at 32 entries, made `bem.fold` F-ordered and single-LU-at-a-time, and added `CSNE` `nmax^2` + cgroup RAM awareness to `cube.solve` and `tetrahedron.solve`.
- **3.3 Fifth benchmark geometry + 20-hour speed, accuracy & figure campaign (20--21 Nov):**
  - **Regular tetrahedron (`mfs/tetrahedron.py`, `figures/tetrahedron.py`, `data/tetrahedron_convergence.npz`, Section 5.6):** Eliminated the 1-depth/octave (`m = 1`) error uptick, cut `m = 2` optical mismatch 13x, and improved symmetric graded Bempp confirmation 2x (`2.19e-4` vs `J = 10` at `N_u = 17424`). Also verified untuned multi-wavelength convergence on the regular tetrahedron at `k = 4*pi` (`1.7e-4` in `79 s`, confirmed by Bempp to `3.2e-4`) and `k = 6*pi` (`1.1e-3`).
  - **Austin RCS exact specification & competitor evaluation (`Section 6`, `Section 7`):** Traced official Austin Problem III-A target length `L = 9.936 in = 0.2523744 m` (exact SI `c_0 = 299,792,458 m/s`), reducing HH mean RCS difference to `0.038, 0.049, 0.065, 0.121 dB` and VV to `0.007, 0.010, 0.011, 0.018 dB` across `3.5--10.25 GHz`. Documented why symmetry-reduced edge-graded OpenCL dense-LU Galerkin EFIE (`Bempp-cl`) is the strongest homogeneous-medium PEC sharp-edge benchmark against volumetric FEM/FDTD+PML, MFIE/smooth Nyström (`burbano2026`), and iterative/Calderón EFIE.
- **3.4 8-hour polyhedral MFS speed & accuracy campaign (21--22 Nov):**
  - Upgraded `mfs/lines.py` (vectorised `build_tree`, projection-centered L2-tiled `_line_kernel_single` / `_line_kernel_multi`, column-outer `_points_kernel`), `mfs/cube.py` and `mfs/tetrahedron.py` (`mv = max(2, m)` vertex-ray spacing fix, universal 2-branch anisotropic boundary quadrature `graded_uv`, root-moment far-target batching in `field()`, class-by-class `CSNE` factorisation with `max_gb = min(6.0, 0.60 * ram)`).
- **3.5 10-hour JCP strengthening campaign: full `O_h` symmetry, zero-symmetry scalene tetrahedron & matrix-free block-QR `LSMR` (22 Nov):**
  - **Full `48`-element `O_h` symmetry in both `mfs/cube.py` and `mfs/bem.py`:** Evaluated kernels on `1/48` of the cube (`y <= z` of face `x = 1/2`), mapped cyclic orbits onto at most 4 representative classes, and factored `C_3` / `Z_2(y <-> z)` blocks in both MFS and Bempp (`O_h` folding shifts Bempp by at most `3e-5` at `4/6` and `7e-6` at `8/10`). Re-timed `cube_convergence` and `cube_cases_{k1,k1o,k2,k3}` across 2 alternating passes: MFS now beats Bempp in every cell of Table 5 (`1.5 vs 2.1 s`, `1.9 vs 2.5 s`, `2.2 vs 5.3 s`, `4.8 vs 12 s` at `10^-3`; `15 vs 148 s`, `25 vs 165 s`, `46 vs 148 s`, `11 vs 147 s` at Bempp's `8 GB` limit), and the combined MFS envelope lies strictly below Bempp throughout (`crossover = []`) in all four cube cases.
  - **Zero-symmetry convex polyhedron (`ScaleneTetrahedron` / `solve_scalene`, `notes/scalene_tetrahedron_bench.py`, `data/scalene_tetrahedron.npz`):** Generalized `cone.CubeVertex` to arbitrary trihedral corners and tested an asymmetric tetrahedron with six distinct edge lengths (`[0.93, 1.09]`), dihedral angles `61.3°--78.3°` (`gamma_e in [0.603, 0.639]`), and four corner cones (`nu_1 in [0.309, 0.342]`) with zero symmetry in either MFS or Bempp. Against `J = 7` (`N_u = 12654`, `506 s`, `opt = 3.7e-4`), `m = 2` error falls monotonically from `5.6e-2` (`J = 2`, `N_u = 3444`, `21 s`) to `1.6e-3` (`J = 6`, `N_u = 9894`, `302 s`), following `8.7 exp(-0.087 sqrt(N_u))`, while non-symmetric graded Bempp (`g = 1.8`, `4/6` and `6/8`) hits its `8 GB` limit at `N_u = 14400` (`n = 40`, `337--556 s`) and agrees with `J = 7` to `1.24e-3`.
  - **Matrix-free coarse-bordered block-QR `LSMR` (`notes/iterative_mfs_bench.py`, `data/iterative_check.npz`, Section 8):** Added `@njit` exact adjoints `eval_points_adj` and `eval_line_adj` (`3.3e-14` transpose exactness) and coarse-bordered edge-wall block-QR preconditioning (`build_bordered_qr`, storing `11--14x` less memory than dense `A` and cutting `cond(A P^-1)` by `4300--7000x` down to `180--400`). Converges to `10^-4` relative normal-equation residual in **`17--22` `LSMR` iterations** on the cube (`J = 4, 6, 8` at `k = 2*pi` and `J = 6` at `k = 4*pi`; `7--53 MiB` vs `88--625 MiB`) and **`33--35` iterations** on the zero-symmetry scalene tetrahedron (`J = 4, 5`; `158--259 MiB` vs `2.10--3.51 GiB`), matching direct solves to `1.3--5.5e-4`.
  - Regenerated all figures and tables, verified `notes/paper_numbers.py` and `19` passing `pytest` tests (`MFS_TEST_BEMPP=1`), re-executed notebooks `1..4`, and compiled `paper/main_3p.pdf` (**strictly 24 journal pages**) and `paper/main.pdf` (**41 preprint pages**) + `paper/mfs-paper-overleaf.zip` with `0` warnings and `0` bad boxes.

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
