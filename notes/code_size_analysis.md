# How big is the code, and can it be simpler? (18 November)

## Outcome

**B is done. D is not done, on purpose (see below).** A later pass that used
library functions where they are as fast (`notes/simplicity_pass_18nov.md`)
brought the package to 1767 lines.

| | before B | after B |
|---|---|---|
| package `mfs/`, total lines | 1872 | **1788** (-84) |
| of which code | 1170 | **1111** (-59) |
| of which docstrings | 402 | 383 |

B saved 84 lines, not the 150 I estimated. The difference is the docstrings:
I kept the argument lists with shapes, as in your notebooks, instead of cutting
each docstring to two lines. Merging the two Bempp solvers into one base class
also saved less than I expected, because the mirror and rotation solves really
are different. Last, I put 8 lines back: three expressions that had become
too dense to read at a glance are now short loops with their shapes commented.

What B changed:

- one symmetry helper (`core.class_part`, `core.classes_of`, `core.sign`) for
  the cube, the almond and Bempp, in place of three copies;
- one base class `bem.Symmetric` for the two symmetric Bempp solvers;
- a shorter QR, without the rank-revealing branch;
- SciPy's `BarycentricInterpolator` in place of a hand-written Lagrange
  interpolation in `lines.py`;
- the almond's source rows and tip tail made by one function, `source_rows`;
- tighter docstrings.

The results do not change:

- the tests pass, including Bempp with symmetry (8 tests, 38 s);
- every table is identical;
- every figure redraws, and the recomputed edge residual agrees with the
  stored one to 4e-13;
- the Bempp fields agree with those before the change to 1e-15;
- the almond's radar cross-section agrees to 1e-13 dB;
- the four notebooks run with no errors.

One known difference: the stored sphere runs used the SVD, and the package
now uses QR. The two give the same errors except where rounding dominates
(R_s <= 0.2 with N = 400), where they differ by up to a factor 2. This is
noted in `figures/sphere.py`.

### Why not D (moving `bem.py` out of the package)

I don't think it is a good idea, and it would not be intuitive:

- **Where readers look.** A reader who wants to know how we ran Bempp will
  look in the package. `pip install "mfs[bempp]"` with
  `from mfs import bem` is the usual pattern for an optional part of a
  package. Code hidden in `figures/` looks like plotting.
- **Who uses it.** The figure scripts and the tests both import it. If it
  sat in `figures/`, the tests would import from a folder of plotting
  scripts, and the scripts would import from each other.
- **The comparison is part of the paper.** The paper's claim is that the MFS
  is faster than Bempp with the same symmetries. A referee who doubts it
  should find the Bempp side as easy to run as the MFS side.
- **It is already marked.** The module docstring and the README both say it
  is the boundary element side of the comparison and not part of the method.
  It is a separate module, and it needs a separate install.

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

## Now (after B)

| | total lines | code | docstrings | comments | blank |
|---|---|---|---|---|---|
| package `mfs/` | 1788 | 1111 | 383 | 55 | 239 |
| figure scripts `figures/` | 848 | 604 | 97 | 38 | 109 |
| tests | 84 | 52 | 5 | 8 | 19 |
| notebooks (code cells) | 150 | 134 | 0 | 7 | 9 |

The plotting that sat inside your notebooks is now in `figures/`, so the
package itself holds only computation.

Lines of code in each module:

| module | what it does | code |
|---|---|---|
| `core` | kernels, plane wave, symmetry helpers, QR and semi-normal equations, cross-sections | 112 |
| `mie` | the Mie series | 34 |
| `sphere` | the sphere | 11 |
| `torus` | the torus, mode by mode | 69 |
| `cube` | layout, tails, symmetry classes, solve | 186 |
| `lines` | the fast line quadrature of the cube (speed only) | 182 |
| `almond` | geometry, source rows and tip tail, compiled kernel, radar cross-section | 189 |
| `cone` | the exponents of Appendix A | 103 |
| `bem` | the Bempp side of the comparison, with the same symmetries | 223 |

## The options as I first wrote them

**A. Keep it as it is.** 1872 lines.

**B. Same results, same speed: about 150 lines fewer.** (Done: 84 fewer.)

**C. Simpler but slower: about 450 lines fewer (to about 1400).** Replace
the fast line quadrature with plain NumPy Gauss panels, the compiled almond
kernel with NumPy, and QR and the semi-normal equations with
`scipy.linalg.lstsq`. The results stay the same, but the cube becomes 5 to
25 times slower and the almond 2 to 5 times. The repository would then no
longer reproduce the timings of Section 7, and the cube comparison would
read worse for the MFS than the paper says. Not recommended.

**D. Move `bem.py` into `figures/`.** Not recommended (see above).

Going below about 1400 lines would mean losing either the speed that the
paper's comparison rests on or the symmetry that makes the large runs fit in
memory. Without the symmetry the cube at J = 10 would need a dense matrix of
about 176 GB. The method itself is already about the size of your thesis's
method code.
