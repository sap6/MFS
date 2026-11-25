# Cover Letter — *Journal of Computational Physics*

**Date:** [Submission Date]

**To:** The Editors, *Journal of Computational Physics*

**Manuscript Title:** *The Method of Fundamental Solutions for Maxwell's Equations*
**Authors:** Saptarshi Bandyopadhyay and Timo Betcke (Department of Mathematics, University College London)

---

Dear Editors,

Please consider our manuscript, *"The Method of Fundamental Solutions for Maxwell's Equations"*, for publication in the *Journal of Computational Physics*.

For smooth three-dimensional scatterers, the Method of Fundamental Solutions (MFS, also known in electromagnetics as the method of auxiliary sources) is attractive because it requires neither a surface mesh nor singular surface integrals and converges rapidly when its interior source surface encloses the singularities of the analytically continued field. On bodies with sharp edges, corners or conical tips, however, those singularities reach the boundary itself, so placing sources on a scaled inner copy of the boundary fails. In two dimensions, corner singularities have recently been resolved with root-exponential convergence by lightning solvers (Gopal and Trefethen, *SIAM J. Numer. Anal.*, 2019), which closed with the remark that the extension to three dimensions remained an open problem. To our knowledge, no MFS formulation for the three-dimensional time-harmonic Maxwell equations has previously achieved convergence on conductors with unrounded edges, polyhedral vertices and sharp conical tips.

This paper develops a unified singularity-guided source-placement framework for three-dimensional perfectly conducting Maxwell scattering across five benchmark geometries (sphere, torus, cube, regular tetrahedron and the NASA almond) and makes four principal contributions of direct interest to the JCP readership:

1. **Medial-surface representation for convex polyhedra (Theorem 1).** Using facewise image continuation and the Stratton–Chu formula, we prove that the field scattered by any convex conducting polyhedron equals the field radiated by integrable electric and magnetic currents on its interior medial surface (skeleton). This explains why scaled inner surfaces stall, shows where interior sources can sit without introducing spurious cavity resonances, and establishes why both electric and magnetic dipoles are necessary.
2. **Stretched Legendre line sources and parameter-free singular tails.** Along polyhedral edges, discrete point sources at depth $d$ would require $O(1/d)$ points per edge; replacing them by continuous Legendre dipole densities in a log-stretched coordinate needs only $O(\ln(1/d))$ modes per line. Below the smallest free depth (and near polyhedral vertices and conical tips), we continue the sources with *singular tails* whose weights follow without fitting from the analytical wedge exponent $\gamma = \pi/\Phi$ and from spherical Laplace–Beltrami conical eigenvalues $(\nu, \tau)$ computed via the subspace-angle method of particular solutions.
3. **Root-exponential convergence on polyhedra and benchmark validation on the NASA almond.** With these layouts, the scattered-field error decays like $\exp(-c\sqrt{N_u})$ on the unit cube (down to $2\times 10^{-5}$, confirmed by graded `Bempp-cl` to $6\times 10^{-5}$, and up to four wavelengths per side and oblique incidence), on the regular tetrahedron (down to $10^{-4}$, confirmed by graded `Bempp-cl` to $2\times 10^{-4}$, testing non-orthogonal $70.53^\circ$ triangular faces, an irrational wedge exponent $\gamma \approx 0.621823$, and a sharper corner $\nu = 0.32221$), and without any symmetry on an asymmetric convex tetrahedron, whereas scaled inner surfaces stall above $20\%$ and $50\%$. On smooth non-convex tori, we conjecture and verify numerically a focal-ring lower bound $a_s > R - \sqrt{R^2 - a^2}$ on the minor radius of the source torus. On the NASA almond, the monostatic radar cross-section agrees with the $397{,}674$-unknown Austin RCS Benchmark reference to $0.007\text{--}0.12$\,dB on average up to $10.25$\,GHz ($8.6$ wavelengths, three minutes per frequency on two cores).
4. **Equal-footing comparison with an edge-graded boundary element solver and matrix-free iterative solves.** Across all five shapes, we compare accuracy and wall time against `Bempp-cl`, equipping the boundary element solver with edge-graded meshes, raised singular quadrature, and the exact same octahedral ($O_h$), mirror, half-turn and rotational symmetry reductions used by the MFS. In addition, column-scaled skeleton and mid-plane matrices remain moderately conditioned ($\operatorname{cond}(A) < 3\times 10^7$), permitting fast Cholesky semi-normal equations with iterative refinement as well as matrix-free coarse-bordered block-QR `LSMR` solves beyond $8$\,GB of memory on zero-symmetry polyhedra.

All source code, precomputed benchmark data, figure-generation scripts and browser notebooks reproducing the results of this paper are publicly available at `https://github.com/sap6/MFS`.

This manuscript is original, has not been published previously (early preliminary results on three of the shapes appeared in the first author's unpublished MSci project report at UCL, which is cited in the text), and is not under consideration for publication elsewhere. Both authors have approved the manuscript and declare no competing interests.

Thank you for your time and consideration.

Yours sincerely,

Saptarshi Bandyopadhyay (`saptarshi.bandyopadhyay.21@ucl.ac.uk`)
Timo Betcke (`t.betcke@ucl.ac.uk`)
Department of Mathematics, University College London, Gower Street, London WC1E 6BT, United Kingdom

---

### Optional Suggested Reviewers (for Editorial Manager)
*(Check for any recent co-authorship or institutional conflicts before entering in Editorial Manager)*
1. **Alex H. Barnett** (Flatiron Institute, Simons Foundation) — Fundamental contributions to MFS analytic-continuation theory, toroidal/axisymmetric boundary methods, and high-accuracy wave scattering.
2. **Lloyd N. Trefethen** (Harvard University / University of Oxford) or **Pablo D. Brubeck** (University of Oxford) — Pioneers of lightning solvers and root-exponential rational/fundamental-solution approximations for corner and edge singularities.
3. **Oscar P. Bruno** (California Institute of Technology) — High-order computational electromagnetics, singular surface/edge scattering, and NASA almond RCS benchmarking.
4. **Shidong Jiang** (Flatiron Institute) or **Zydrunas Gimbutas** (NIST) — Integral equations, dipole/multipole representations, and singular asymptotics for Maxwell's equations.
5. **Simon N. Chandler-Wilde** (University of Reading) or **Jeffrey Ovall** (Portland State University) — Numerical analysis of polyhedral/wedge wave scattering and conical/corner eigenvalue problems.
