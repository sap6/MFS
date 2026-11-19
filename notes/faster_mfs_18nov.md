# Making the MFS faster than Bempp at larger errors, fairly (18–19 November)

## Where it stands

The MFS on the cube now becomes faster than Bempp at errors 2.5 to 3 times
larger than in the paper. These are measured today and placed against Bempp's
stored best runs, which I re-timed and found within 5%:

| case | the paper: MFS faster below | now: MFS faster below | MFS / Bempp time at 3e-4 |
|---|---|---|---|
| standard wave | 1.5e-4 | **4.1e-4** | 3.3 s / 6.7 s |
| general wave | 1.4e-4 | **5.1e-4** | 6.3 s / 11.2 s |
| two wavelengths | 3.0e-4 | **7.3e-4** | 5.4 s / 27.7 s |
| three wavelengths | about 5e-4 | **1.6e-3** | 9.4 s / 35 s (at 5e-4) |

At the accurate end the MFS also got faster. On the standard wave, J = 11
(5.5e-6) takes 48 s instead of 108 s. At three wavelengths, J = 9 takes 120 s
instead of 258 s.

Status (19 November): the changes are committed. Re-run and committed so far:
Figure 6, Table 2, Figure 10(d), Figure 5, and the almond timings. Figure 7 and
Table 5 (`figures/cube_cases.py --run`, about 5 hours) are still to be re-run:
see `notes/HANDOFF.md`.

## What changed (about 10 lines of code)

1. **The boundary points are graded as the paper says.** Section 5 says the
   panels "halve in length towards both edges until they are shorter than
   d_min". The code graded two levels further, down to about d_min/6. With the
   rule as written, the errors are the same or slightly smaller and the cube is
   about 1.3 times faster. At J = 12 the result agrees with the stored reference
   to 3e-7. One line.
2. **No duplicate work in the kernels.** A line source is its own mirror
   image across the plane through its middle. So four of the eight mirror
   images of a boundary point repeat the other four, up to a sign (−1)^n for
   mode n. The parity class keeps exactly the columns where the two copies add.
   The kernel now computes four images and doubles them: the fields are the same
   to 1e-13, and the cube is 1.25 to 1.5 times faster. The almond's sources lie
   on its mid-plane z = 0, so the same applies there. Two lines in each kernel.
3. **The best MFS layout at each accuracy, as for Bempp.** Bempp's curve is
   the best of about 150 runs (mesh, grading, quadrature). The MFS now
   contributes both of its layouts, two depths per octave and one depth per
   octave. The paper already describes both and already remarks that the
   second one moves the crossover. Figure 7 draws both MFS curves, as
   Figure 10(d) already does. Table 5 takes the faster one.
4. **The almond timing is now like Bempp's.** The timed almond run included the
   far field for the optical-theorem check (about 10% of its time). Bempp's
   timed run does not compute it. The check now runs outside the timer
   (`almond.optical`), as it already did for the cube and the torus.

## Fairness

- Bempp keeps everything it had: the symmetry, the compiled fold, LU, and the
  best of its scan.
- The duplicate images exist only in the MFS. Bempp has no element that is its
  own mirror image, and its symmetric assembly already computes each distinct
  interaction once.
- I re-timed Bempp's best runs between 1.5e-4 and 1e-3 today: within 5% of
  the stored times. The re-run times both codes again, alternately.

## Tried and not adopted

| idea | result | why not |
|---|---|---|
| interior panels 0.15 instead of 0.1 | standard wave 1.3 times faster (crossover 5.1e-4) | J = 12 then misses the reference (2.3e-5); at three wavelengths it fails (errors 1e-2 to 1.5); it would need a rule that depends on the wavelength |
| 10/6 quadrature nodes instead of 12/8 | 10 to 15% faster, field changes 4e-7 | changes Appendix B for little gain |
| 5 tail lines instead of 9 | 5 to 15% faster, errors within 3% | changes the description of the tails for little gain |
| depth ratio 1/4 in the tails | errors 2 to 3 times larger | |
| grading 3 levels short of the current code | 15% faster still, errors as good | one level from the cliff (4 levels short gives errors near 1); the rule as written is safer |
| plain NumPy for the kernels, or SciPy `lstsq` | 20 to 50 times, and 2 to 4 times, slower | see `notes/speed_vs_simplicity_18nov.md` |

The one-depth layout stops converging above about J = 8 at two and three
wavelengths. Figure 7 shows its curve as it is. The comparison uses whichever
MFS run is fastest at each accuracy, just as Bempp's curve uses its fastest runs.

## What the paper will need once the re-run is done

- Figures 3, 5, 6, 7 and 10(d), and Tables 2, 5 and C.1, come from the
  scripts.
- Section 7 quotes the crossovers ("level near 1.5e-4 ... 3e-4 at two and
  5e-4 at three") and some times (24 s, 20 s, 30 s, "2e-5 in less than a
  minute"). Those numbers will change. I will list each one with the new
  value.
- Section 5's description of the boundary points stays as it is: the code
  now follows it.
