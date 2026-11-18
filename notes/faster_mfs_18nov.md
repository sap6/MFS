# Making the MFS faster than Bempp at larger errors, fairly (18 November)

## Summary

On the cube the MFS can overtake Bempp at about twice the error it does now:

- standard wave: below **3.2e-4** instead of 1.5e-4;
- general wave: below about **3e-4** instead of 1.4e-4.

None of this needs new code or anything taken away from Bempp. It needs two
changes:

1. **Compare the best MFS layout at each accuracy, as we already do for
   Bempp.** Bempp's curve is the lower envelope of a scan over 7 gradings, 3
   pairs of quadrature orders and 7 or 8 meshes. The MFS curve is a single
   layout (two source depths per octave). The code already has a second
   layout, one depth per octave (`q=0.5`), and Figure 10(d) already plots it.
   It is three times faster at moderate accuracy, but it stalls near 1.3e-4.
   The fair MFS curve is the lower envelope of the two layouts.
2. **Optional: a coarser interior collocation grid**, (n_base, hmax) = (6, 0.15)
   instead of (8, 0.1). It gives the same errors and is 20 to 45% faster up to
   J = 10. The research log (item 65) shows that J = 11 and 12 need (8, 0.1),
   so it would be a second rule.

## Fairness

- Bempp keeps everything it has: the symmetry, the compiled fold, LU, and the
  best of its scan.
- I re-timed Bempp's best runs between 1.5e-4 and 1e-3 today, alternating
  with the MFS, as the paper's protocol does. They match the stored times to
  within 5% (for example 10.9 s for 2.1e-4, and 19.9 s for 1.5e-4).
- Choosing the best of the MFS's two layouts is the counterpart of choosing
  the best of Bempp's 150 configurations.

## Numbers (standard wave; errors against the J = 12 reference)

Measured today, each run twice, keeping the smaller time:

| | one depth per octave (8, 0.1) | one depth (6, 0.15) | two depths (8, 0.1), the paper's | two depths (6, 0.15) | Bempp's best at that error |
|---|---|---|---|---|---|
| J = 6 | 3.7 s, 9.4e-4 | 2.7 s, 9.6e-4 | 10.1 s, 8.3e-4 | 7.6 s, 8.3e-4 | about 1.2 s |
| J = 7 | 5.3 s, 4.1e-4 | 4.2 s, 3.8e-4 | 16.0 s, 3.4e-4 | 12.2 s, 3.4e-4 | about 3 s |
| J = 8 | 7.2 s, 2.1e-4 | 5.8 s, 1.8e-4 | 23.1 s, 1.4e-4 | 17.5 s, 1.4e-4 | 10.9 s at 2.1e-4 |
| J = 9 | 10.1 s, 1.35e-4 | 8.1 s, 1.5e-4 | 37.5 s, 5.6e-5 (stored) | | 19.9 s at 1.5e-4 |

Where the MFS becomes faster than Bempp:

| MFS curve | standard wave | general wave |
|---|---|---|
| the paper's (two depths, (8, 0.1)) | below 1.5e-4 | below 1.4e-4 |
| best of both layouts, (8, 0.1) | below 3.2e-4 | level at 3e-4 (11.3 s against 11.2 s), faster below, except 1.4e-4 to 1.7e-4 |
| best of both layouts, (6, 0.15) | below 3.3e-4 | the band 1.4e-4 to 1.7e-4 closes: two depths at J = 8 take 29.2 s instead of 40.8 s (Bempp about 34 s) |

At 2e-4 on the standard wave the MFS would need 5.5 to 7.3 s, against 12.6 s
for Bempp. Above about 3e-4 Bempp stays faster (it needs 1 to 3 s there).

## What did not work

- (n_base, hmax) = (4, 0.25): the error grows from J = 6 to 7 on (8e-4 at J = 7
  with two depths; 2e-2 at J = 8 with one depth).
- Fewer nodes per panel: already known to fail (item 47).
- Fewer modes: already known not to pay (item 56).

## What it would mean for the paper (your decision)

- Figure 7 and Table 5 would show the one-depth layout as well, and the text
  would give the new crossovers. The standard and general waves are measured
  above. The two- and three-wavelength cases (k2, k3) have not been tried with
  one depth per octave; those runs take about an hour each.
- The paper would describe the comparison as "the best MFS layout against the
  best Bempp configuration at each accuracy". Section 5 already describes both
  layouts.
- The grid (6, 0.15) adds a rule (it holds up to J = 10). It gains another 20
  to 30%. I would add it only if you want every percent; the envelope alone
  doubles the crossover.

Not tried within the hour: rewriting the assembly kernel. Accumulating the
modes costs about 300 operations per quadrature node, and it might run faster
as small matrix products. That would be our own code again, so it is the kind
of trade you asked about earlier.
