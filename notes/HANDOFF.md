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

  Measured crossovers with Bempp (the MFS is faster below):

  | case | now | before |
  |---|---|---|
  | standard wave | 4.1e-4 | 1.5e-4 |
  | general wave | 5.1e-4 | 1.4e-4 |
  | two wavelengths | 7.3e-4 | 3.0e-4 |
  | three wavelengths | 1.6e-3 | 5e-4 |

  Bempp's stored times were re-checked today and agree within 5%.

## 3. Unfinished: do this first

**3.1 Re-run the data of Figure 7 and Table 5** (both codes, alternately,
both MFS layouts). The script is ready. The stored rows are still the old
two-depth runs, relabelled to the new format; the figure is unchanged until
the re-run.

    python figures/cube_cases.py --run k1     # then k1o, k2, k3: about 1 to 1.5 hours each
    python figures/cube_cases.py && python figures/cost.py && python figures/tables.py

- Run nothing else at the same time: the times are the result.
- Background tasks in this environment are stopped after 2 hours, and the
  container restarted once. So start the runs detached, one case at a time
  (`setsid nohup ... &`), check on them, and commit each case when it ends.
- Each run first writes `data/cube_cases_<case>_ref.npz` (the reference),
  then `data/cube_cases_<case>.npz` at the end.
- Already re-run and committed: `cube_convergence` (Figure 6, Table 2,
  Figure 10d), `edge_residual` (Figure 5), the almond timings (Figure 10c,
  Table 4).

**3.2 Update the numbers in the text.** The tables are generated; the text is
not. `python notes/paper_numbers.py` prints every number of Sections 5 and 7
that depends on the cube runs, from `data/`. Change only numbers and leave the prose to the author's style.
List every change for the author.

- **Section 5.4** (known now, from `data/cube_convergence.npz`):
  - "from 5.1e-3 at J = 4 to 2.1e-5 at J = 10" becomes 5.3e-3 and 2.0e-5;
    the factor of 2.4 to 2.7 per level is unchanged.
  - Equation (rootexp): 0.46 and 0.0474 become **0.48 and 0.0480** (within
    5%).
  - J = 11: 5.8e-6 ("6e-6" stays).
  - One depth per octave: 56 to 68% of the unknowns, so "about 60%" stays;
    "between a third and a half of the time" becomes **between a quarter and
    a half**; "within 20% down to J = 7" stays.
  - Without tails, 4.9e-2, 2.2e-2 and 9.5e-3: unchanged.
  - Reference: J = 11 and 12 differ by 5.8e-6, and its optical mismatch is
    8.6e-7 ("6e-6" and "9e-7" stay; the error of the reference is still
    about 4e-6).
- **Section 5.4, after 3.1:** "Bempp's finest mesh differs from the
  reference by 6e-5" and the extinction cross-sections (2.92786 against
  2.92798).
- **Section 5.5, after 3.1:** the general-wave, two- and three-wavelength
  errors per level, the cost ratios, and the references' own errors.
- **Section 7, after 3.1:**
  - the paragraph "In between, the two methods are level near 1.5e-4 ... 3e-4
    at two and 5e-4 at three ... 24 s ... 20 s ... 30 s ... 2e-5 in less than
    a minute";
  - the paragraph on one depth per octave ("near 3e-4 ... 8 s ... 11 s ...
    11 and 30 s");
  - "At an error of about 1e-3, the boundary element method is the faster ...
    ten times at one wavelength, six times at two and four times at three"
    (at three wavelengths the MFS may now be faster);
  - "its most accurate run ... the MFS reaches a smaller error in half of
    that time or less, with J = 9 ... J = 8 at three".

  The comparison should now be described as "the best MFS layout against the
  best Bempp configuration at each accuracy".
- **Figure 7 caption:** mention the two MFS curves (two depths per octave,
  and one depth per octave).
- **Conclusion:** check "confirms the field to 6e-5" and "about 4e-6" after
  3.1.

**3.3 Re-check the claims from separate runs that used the old grading**
(Section 5.4), and change the text only if the new number differs:

- one depth per octave stalls between 1.3e-4 and 1.8e-4 at J = 10 to 12;
- at J = 9 the field changes by less than 3e-6 when the modes, the tails, the
  base grid or the boundary grid are refined;
- adding the next term of each edge expansion to the tail changes the error
  by less than 5%, and lowering the bisector exponents by one makes it 4
  times worse at J = 6 and 10 times at J = 8 (`cube.TAIL_E`);
- at two wavelengths, J = 8 changes by 8e-6 with four more modes;
- at three wavelengths, 8 or 16 modes in place of 12 change the field by
  3e-4 or 2e-4, and 6e-5.

**3.4 Re-run the notebooks** `notebooks/3_cube.ipynb` and `4_almond.ipynb`
(their times changed):

    jupyter nbconvert --to notebook --execute --inplace notebooks/3_cube.ipynb notebooks/4_almond.ipynb

**3.5 Rebuild the paper:**

    cd paper && latexmk -pdf main.tex

Also build `main_3p.tex` (the journal layout), read the result and send the
PDF to the author.

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
