# Handoff (19 November 2026)

Read this first, then `CLAUDE.md` and the top of `notes/NOTES.md`.

## 1. State of the repository

- Branch `ccr-116c520c-fts35f` of `sap6/mfs`, everything committed and
  pushed (last commit before this file: `9c7b3d4`). No pull request; not
  merged into `main`. Nothing was ever pushed to `sap6/mfs-paper` (it syncs
  to Overleaf): never push there without asking the author.
- `python -m pytest tests` passes (8 tests with `MFS_TEST_BEMPP=1`), CI on
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
  | standard wave (`k1`) | 4.5e-4 | 1.8e-4 |
  | general wave (`k1o`) | 3.9e-4 | 1.5e-4 |
  | two wavelengths (`k2`) | 8.6e-4 | 4.0e-4 |
  | three wavelengths (`k3`) | 2.2e-3 | 1.0e-3 |

## 3. Completed on 19 November 2026

- **3.1 Re-run of Figure 7 and Table 5 (`data/cube_cases_{k1,k1o,k2,k3}.npz`):** Done in 2 alternating passes on 2 cores (8 GB RAM), including `k3` after fixing peak-memory copies in `core.QR` / `core.CSNE` and setting `cube.solve` to assemble one class at a time when `k > 5 pi` and `J >= 8`. Tables (`paper/tables/*.tex`) and all 11 figures regenerated and visually checked.
- **3.2 Numbers in the text (`paper/sections/*.tex`):** Updated throughout from `notes/paper_numbers.py` and verified by an independent review pass.
- **3.3 Standalone claim re-checks with the new boundary grading:**
  - One depth per octave at `J = 10, 11, 12`: `1.27e-4`, `1.44e-4`, `1.78e-4` (stalls between `1.3e-4` and `1.8e-4`, unchanged).
  - Refinements at `J = 9` vs `dense_E9`: `+4` modes `9.9e-7`, `12x12` base `1.1e-7`, `14` tail steps `7.3e-7`, `+2` boundary edge levels `3.43e-6` (`+2` levels & `hmax = 0.07`: `3.43e-6`) -> updated `< 3e-6` to `< 4e-6`.
  - Tail exponents at `J = 6, 8`: adding next edge term changes error by `-1.6%` (`J = 6`) and `-3.1%` (`J = 8`, `< 5%`); lowering bisector exponents by 1 increases error by `3.37x` at `J = 6` and `7.60x` at `J = 8` (updated in `cube.tex`); adding 2nd vertex exponent at `J = 6` changes error by `-0.15%` (`< 1%`).
  - Two wavelengths (`k = 4 pi`), `J = 8` (`+8` vs `+4` modes): `9.5e-6` (`1e-5`).
  - Three wavelengths (`k = 6 pi`): `+8` vs `+12` modes changes `J = 7` by `3.18e-4` (`3e-4`) and `J = 8` by `1.92e-4` (`2e-4`); `+16` vs `+12` modes changes `J = 7` by `7.36e-5` (`7e-5`).
  - Almond at `3.5 GHz` (`s = 0.5`, 2 alternating passes): MFS `5:2` reaches `4.2e-5` in `5.7 s` (`6:3`: `8.3e-6` in `9.3 s`), while Bempp reaches `3.1e-3` in `8.4 s` and `1.1e-3` in `64 s`.
- **3.4 Notebooks (`notebooks/1_sphere.ipynb` .. `4_almond.ipynb`):** Cell 1 updated for both Colab and local execution; all four re-executed cleanly.
- **3.5 Manuscript & Overleaf bundle:** `paper/main.pdf` (45 pages, preprint), `main_3p.pdf` (26 pages, JCP `3p`), and `paper/mfs-paper-overleaf.zip` built with zero warnings or overfull boxes; all appendix hyperlinks fixed (`\theHsection`).

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
- **QR fallback:** CSNE refuses the matrix only at three wavelengths with
  J >= 9; QR then takes over, and needs one more copy of the matrix.
- **The author's voice:** `notes/style_profile.md` (from `thesis/`). The
  title is settled. The abstract has at most 250 words.
