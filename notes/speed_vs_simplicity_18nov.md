# Is the speed worth the code? (18 November)

## Short answer

Yes, for each piece I measured. Each one is 10 to 30 lines. Dropping any one
of them would slow down the method it serves by a factor between 1.2 and 50.
The two largest would also change what the paper can claim against Bempp:

- **The semi-normal equations (CSNE).** No library has them. The library
  routes are QR (the same speed as our QR) or SciPy's `lstsq` (slower still).
- **The precomputed moments in the line quadrature.** These are our own; no
  library does this quadrature.

Bempp keeps all of its own speed-ups (LU, symmetry, the compiled fold, tuned
gradings and quadrature orders), so removing ours would no longer compare the
best of both codes.

## Can a library skip the extra work?

- **SciPy `lstsq`: no.** Its `lapack_driver` offers only `gelsd` and `gelss`
  (an SVD) and `gelsy` (QR with column pivoting). `cond` only changes which
  singular values are dropped, not the work done. No other option avoids the
  extra work.
- **LAPACK's plain QR driver `gels`: yes, but only from the low-level
  interface.** SciPy has it as `scipy.linalg.lapack.zgels`. With the optimal
  workspace (from `zgels_lwork`) it is as fast as our QR: 3.95 s against
  3.83 s. With SciPy's default workspace it takes 23 s, because LAPACK then
  falls back to its slow code. It would not shorten our code: a wrapper is as
  long as our QR class, and it factorises again for every wave the cube solves.
- **NumPy `lstsq`:** an SVD only.
- **PyTorch** offers `gels` through `torch.linalg.lstsq`, but it would add a
  large dependency for one call, and it is still QR.
- **No library has the semi-normal equations**, so nothing reaches CSNE's
  speed: 1.8 s against 3.9 s for QR on a 9000 × 3000 matrix.
- **No library computes these matrices.** Fast-multipole libraries give
  matrix-vector products, but least squares needs the matrix entries.

## What each speed-up costs and buys

Cube times are totals for the standard wave, measured now on this machine.

| our code | lines | without it | effect on the paper |
|---|---|---|---|
| **CSNE** (core) | about 25 | QR: cube 1.4 to 1.9 times slower (J = 4 to 9); at J = 10 QR ran out of memory here, because it needs one more copy of the matrix. SciPy `lstsq`: 1.5 to 2.6 times slower | standard wave, at Bempp's best accuracy (6e-5 in 131 s): the MFS takes 69 s (QR) or 93 s (`lstsq`) instead of 36 s, and is faster below 8e-5 or 7e-5 instead of below 1.5e-4. General wave with `lstsq`: the MFS is slower than Bempp at every accuracy that both reach |
| **moments on the coarsest panels** of the line quadrature | about 25 | Gauss panels only (the fields agree to 3e-9): assembly twice as long, cube 1.6 times slower (J = 8: 38 s instead of 24 s) | about as large as dropping CSNE |
| **compiled point-dipole kernel** (cube) | about 30 | NumPy: 3.7 to 6.7 s per assembly instead of 0.2 s | cube at least 17% slower |
| **compiled almond kernel** | about 20 | NumPy: 16.8 s per class instead of 0.31 s | the 6:3 layout takes about 20 s instead of 3.7 s; Bempp's best takes 22 s |
| **compiled fold** (Bempp side) | about 10 | NumPy | Bempp 33.8 s instead of 20.6 s (recorded in item 72) |

**Both cube speed-ups dropped** (QR, Gauss panels only): the cube is 2.1 to
2.4 times slower.

- Standard wave: at Bempp's best accuracy the MFS takes 87 s instead of 36 s,
  and it is faster only below 7e-5.
- General wave: the MFS is slower than Bempp at every accuracy that both reach
  (143 s against 133 s at Bempp's best).

The MFS would still reach accuracies that Bempp cannot reach (6e-6 on the
standard wave). But the paper's statement that the MFS is cheaper on the cube
at high accuracy would no longer hold for the general wave.

## Recommendation

Keep all of them. They add about 110 lines in total. They are what lets the
paper compare the best MFS with the best Bempp, as you asked. The paper already
describes CSNE (Section 2) and the moments (Appendix B) in a few sentences,
and the tests check the fast code against the paper's stored results.

One small change I did not make: the QR class could call LAPACK's `gels`
(`scipy.linalg.lapack.zgels`, with the optimal workspace) instead of the two
calls it makes now. This would gain no speed and no lines, and the cube could no longer reuse
a factorisation, so I left it.
