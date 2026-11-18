# Simplicity pass: where a library could do the work (18 November)

## Summary

The package already uses a library wherever one does the same job at the same
speed. The code still written by hand falls into two groups:

- code that is faster than any library route I found (the least-squares solve
  and the compiled kernels);
- code that is specific to the paper (source layouts, boundary points, the
  symmetry, the meshes).

I made six replacements. The package went from 1788 to 1767 lines, and every
result is the same as before.

## The least-squares solve

You asked whether SciPy's own least squares, perhaps with some parameter,
would do as well. I timed every route on a matrix the size of a cube class
(9000 × 3000 complex, four right-hand sides, four cores):

| method | time |
|---|---|
| ours: semi-normal equations (CSNE) | **2.1 s** |
| ours: QR without forming Q (LAPACK geqrf, unmqr) | 4.1 s |
| `scipy.linalg.qr_multiply` | 4.3 s |
| `scipy.linalg.qr(mode="economic")`, which forms Q | 8.2 s |
| `numpy.linalg.lstsq` | 8.5 s |
| `scipy.linalg.lstsq`, default driver `gelsd` | 8.6 s |
| `scipy.linalg.lstsq`, driver `gelsy` | 10.9 s |
| `scipy.linalg.lstsq`, driver `gelss` | 182 s |

All three `lstsq` drivers are built for matrices that may not have full rank.
`gelsd` and `gelss` compute an SVD. `gelsy` uses QR with column pivoting,
whose search for pivots runs at matrix-vector speed. No parameter turns this
off. Our matrices have full rank once their columns are scaled, so plain QR is
enough.

`qr_multiply` makes the same LAPACK calls as our QR, at the same speed, but it
factorises the matrix again for every right-hand side. The cube factorises each
class once and then solves it for up to three waves (the wave turned about the
body diagonal), so it needs to keep the factorisation.

SciPy has no semi-normal equations, and they are where the time goes. At J = 9
the cube takes 39 s, of which 18 s is the factorisation, and CSNE costs half as
much as QR.

**Both solvers are needed.** I logged which one the paper's runs use:

- CSNE solves every cube and almond run;
- the one exception is the cube at three wavelengths with J >= 9 (Figure 7,
  Table 5), where CSNE refuses the matrix and QR takes over.

**What changed.** The two solvers now use SciPy's own functions in place of
bare LAPACK calls:

- `scipy.linalg.qr(mode="raw")` for QR without forming Q;
- `scipy.linalg.cho_factor` and `cho_solve` for the Cholesky factor.

Applying Q^H still needs LAPACK's `unmqr`, because SciPy offers no other way to
reach it. The Gram matrix A^H A stays a BLAS call (`zherk`). Plain NumPy
(`A.conj().T @ A`) takes 5.4 s instead of 3.0 s on a 12000 × 4000 matrix, and
it makes a copy of A, which is several GB at J = 12. The coefficients and
residuals are identical to the last bit.

## The rest of the package

### Replaced

| where | before | now | check |
|---|---|---|---|
| `core.magnetic` | six assignments build the matrix of β → e × β | `np.cross(np.eye(3), e)` | identical |
| `bem.efie` | LU factorisation, then LU solve | `scipy.linalg.solve` (the same LAPACK routine) | identical |
| `cone.smallest_sigma` | SVD, then cut off the small singular values | `scipy.linalg.orth` | singular values to 8e-12; Table A.1 unchanged |
| `almond.matrix` | its own copy of the dipole formulas | the dipole function the cube's kernels use | fields to 3e-16, RCS to 6e-14 dB |

The last one uses no library, but the idea is the same: write the formula once.
It is 25 lines shorter, and the kernel itself is 1.4 to 1.7 times faster. The
total time of the almond does not change (3.7 s for the 6:3 layout and 8.5 s
for 8:5, before and after), so the paper's timings stand.

### Already library calls

- NumPy: the torus FFT, Gauss-Legendre nodes, Legendre polynomials.
- SciPy: the spherical Bessel functions of the Mie series, the Legendre
  functions of the cone, the bounded minimisation, barycentric interpolation,
  `lstsq` for the torus, LU for Bempp, sparse matrices.
- Bempp: assembly, potentials, the sphere mesh.
- numba for the kernels, matplotlib for the figures.

### Considered, not done

| what | library | why not |
|---|---|---|
| the cube's compiled point-dipole kernel (30 lines) | NumPy (`core.dyadic`) | 20 to 30 times slower (0.2 s becomes 3.7 to 6.7 s per assembly); the cube would be at least 17% slower |
| the line quadrature (`lines.py`) | `scipy.integrate.quad_vec` | it calls Python at every node; the tree, with the solver and the third mirror plane, made the cube 25 times faster (notes, item 41) |
| the Gram matrix | NumPy | 1.8 times slower, and it copies A |
| the quadrature on the sphere (7 lines) | `scipy.integrate.lebedev_rule` | the cross-sections agree to 1e-13, but the same rule places the 200 test points of the torus behind Table 4, so the torus numbers would change in the last digit; it also needs SciPy 1.15 |
| merging the vertices in `bem.mirror_grid` | `np.unique` | saves about 3 lines, but it sorts and rounds the vertices, so the Bempp meshes and their numbering change |
| the meshes of the cube, torus and almond | `bempp.shapes` (gmsh) | the symmetry needs structured meshes made of one piece and its images; `bempp.shapes` gives unstructured meshes (the sphere does use it) |
| the symmetry in Bempp | none | Bempp has no support for symmetry |
| the Mie series (34 lines) | none | it is the independent reference for every sphere result; SciPy already supplies the Bessel functions |

## Checks

- tests, including Bempp with symmetry (8 tests);
- every figure redraws and every table is identical;
- almond nose-on field, extinction and scattering identical; RCS at 3.5 GHz
  (tip tail, all four classes) to 6e-14 dB;
- sphere coefficients and torus optical theorem identical;
- the Bempp sphere EFIE identical.

The paper's text on the solver (Section 2) is still accurate: QR for the
sphere, the SVD for the torus, CSNE for the cube and the almond with QR as
the fallback.
