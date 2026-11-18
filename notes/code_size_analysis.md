# How big is the code, and can it be simpler? (18 November)

## Your original

The seven thesis notebooks hold **1851 lines** of code cells:

- 1312 lines of code: 801 for the method and 511 for plotting;
- the remaining 539 lines are comments, docstrings and blank lines.

The four Maxwell notebooks alone hold **1299 lines** (924 of code: 579 for the
method, 345 for plotting).

Much of that method code is repeated. `fs()` (the dyadic Green's function)
and `E()` (the plane wave) are written out in every Maxwell notebook. Each
shape was an interior problem with sources on a scaled copy of the surface,
solved with `lstsq`.

Before today the paper's code was `validation/`: 67 scripts and **7301 lines**.

## Now

| | total lines | code | docstrings | comments | blank |
|---|---|---|---|---|---|
| package `mfs/` | 1872 | 1170 | 402 | 53 | 247 |
| figure scripts `figures/` | 859 | 605 | 109 | 38 | 107 |
| tests | 85 | 52 | 5 | 8 | 20 |
| notebooks (code cells) | 147 | 131 | 0 | 7 | 9 |

The plotting that sat inside your notebooks is now in `figures/`, so the
package itself holds only computation.

## Where the 1170 lines of package code go

1. **The method, as the paper describes it: about 640 lines.** The thesis
   method code was 579 lines, so this is about the same size, for a harder
   problem: exterior scattering, the edge lines and tails, the tip tail, and
   the symmetry.

   | | lines |
   |---|---|
   | kernels and helpers | 62 |
   | Mie series | 34 |
   | sphere | 11 |
   | torus | 69 |
   | cube: layout, tails, symmetry classes, solve | 194 |
   | almond: geometry, source rows, tip tail, boundary points, RCS | 165 |
   | cone exponents (Appendix A) | 103 |

   The symmetry is not optional. Without it the cube at J = 10 would need
   a dense matrix of about 176 GB, and the finest torus runs one of about
   85 GB.

2. **Speed: about 270 lines.** These make the cube 25 times faster with
   identical fields (notes item 41).

   | | lines |
   |---|---|
   | the fast line quadrature, `lines.py` | 191 |
   | the compiled kernel of the almond | 43 |
   | QR without forming Q, and the semi-normal equations | 42 |

   Section 7 compares wall times with Bempp, and these timings depend on
   this code. Without it the MFS would be slower than Bempp on the cube
   at every accuracy.

3. **The Bempp comparison: 261 lines** (`bem.py`). This is the boundary
   element side, with the same symmetries as the MFS. It is not part of the
   method.

## Options

**A. Keep it as it is.** 1872 lines.

**B. Same results, same speed: about 150 lines fewer (to about 1700).**

- Shorten the module docstrings to two or three lines with a pointer to the
  paper section (about -60). They now repeat the paper.
- Merge the two symmetric Bempp solvers (mirrors and rotation) into one
  (about -40).
- Use one symmetry helper for the cube, the almond and Bempp (about -15).
- Smaller merges elsewhere (about -30).

The tests would show that nothing changes.

**C. Simpler but slower: about 450 lines fewer (to about 1400).**

- Do everything in B.
- Replace the fast line quadrature with plain NumPy Gauss panels.
- Replace the compiled almond kernel with NumPy.
- Replace QR and the semi-normal equations with `scipy.linalg.lstsq`.

The results stay the same, but the cube becomes 5 to 25 times slower and
the almond 2 to 5 times. The repository would then no longer reproduce the
timings of Section 7, and the cube comparison would read worse for the MFS
than the paper says.

**D. Move `bem.py` into `figures/`.** This adds or removes no lines, but
the package would then hold only the method (about 1600 lines, or about 1450
with B). A reader would see that the comparison code is separate.

## My recommendation

B and D. Going below about 1400 lines would mean losing either the speed the
paper's comparison rests on or the symmetry that makes the large runs fit in
memory. The method itself is already about the size of your thesis's method
code.
