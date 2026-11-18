# The repository for readers of the paper (18 November)

## What a visitor sees

```
README.md      the landing page: results figure, Colab table, install,
               which script makes which figure
mfs/           the method, 1862 lines, one module per shape
notebooks/     one notebook per shape; opens in Colab with one click
figures/       one script per figure of the paper, and the tables
data/          the results behind every figure and table (1.6 MB)
tests/         checks against exact solutions and the paper's results
paper/         the manuscript
thesis/        your MSci notebooks, as submitted
notes/         our research log
```

The old `validation/` folder had 67 scripts and about 7300 lines. It is now
1862 lines in the package and 846 in the figure scripts. Every module
reproduces the stored results of the paper:

- the sphere errors are the same;
- the torus fields are identical;
- the cube fields agree to 1e-14, and the cube runs at the same speed;
- the almond fields agree to 1e-15, and its radar cross-section to 1e-13 dB;
- the exponents of Table A.1 are the same;
- the Bempp runs with symmetry are identical to the old ones.

Every table the scripts generate is identical to the old one. Most figures
are pixel-identical. Figures 2, 3 and 9 changed slightly: Figure 3 now draws
the real boundary points of the J = 6 layout.

The comments follow your notebooks: a docstring that starts with the formula,
then the arguments with their shapes; short block comments such as
"# Form matrix A"; and aligned comments on constants and shapes. A fresh
reviewer compared the code with your thesis notebooks. I applied its
findings, including two real bugs:

- the Colab cell could not import the package;
- `cube.solve` did not free memory between groups.

## Are Jupyter notebooks the best way to show this?

On their own, no. They are the best entry point, though, so the repository
uses three layers:

- **Notebooks, for trying things.** A reader clicks the Colab badge and sees
  the sphere, torus, cube and almond results in a few minutes in the browser,
  with nothing to install. GitHub also shows each notebook with its pictures.
- **One plain script per figure, for checking a result.** A reader who asks
  "how was Figure 6 made?" finds a short `.py` file. It is easier to read
  and run than a notebook, has no hidden state, and is tested on every push.
- **A small package, for the method itself.** Both the notebooks and the
  scripts import it, so the method is written once.

I looked at newer formats such as marimo, which runs notebooks in the browser
without a server. It cannot run numba there, and the cube needs numba, so
Colab is the better choice for now.

## What only you can do

1. **Merge this branch into `main` without squashing.** The Colab links, the
   CI badge and `pip install git+https://github.com/sap6/MFS` all point to
   `main`. The old research scripts live in the git history at commit
   9f5a721. I made a local tag, `research-scripts`; should I push it?
2. **Licence.** `LICENSE` says "Copyright (c) 2025 Sap". You may want your
   full names and 2025–2026 there.
3. **Zenodo.** Turn on GitHub's Zenodo integration and make a release. That
   gives a DOI to put in the paper's code statement and in `CITATION.cff`.
4. **Public notes.** `CLAUDE.md` and `notes/` contain your messages to me and
   the status sheets. Keep them public, or move them to a private branch
   before the merge.

## Found while porting

The library's default quadrature in Bempp 0.4.2 is 4/4, not 4/6 as
Appendix C said. The cube runs set 4/6, 6/8 and 8/10 themselves, while the
sphere, torus and almond used 4/4. I tested 4/6 and 6/8 on all three shapes:
Bempp's error changes by less than 1%, so the comparison stands. I corrected
Appendix C.

## The paper

- **Abstract:** 242 words. The JCP guide for authors sets a limit of 250.
- **Page cuts:** the only part I found to be truly redundant was the almond's
  "In summary" paragraph, which repeated two sentences from just above it. I
  removed it. The other candidates each carry evidence that the text relies
  on, so they stay:
  - the 2D square shows why the depth ratio is 2^(-1/2);
  - the torus arguments are the support for the conjecture;
  - the edge-residual figure shows what the tails do.

  The paper is 41 preprint pages and 24 in the journal layout.
- **Code statement:** it now describes the repository.
- **Declarations:** the JCP declarations are still to be written by you
  (`notes/jcp_checklist.md`).
