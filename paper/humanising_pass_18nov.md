# Humanising pass, 18 November

You asked for a humanising pass in your own style, not a rewrite. So the
structure, every number, every citation and every claim are unchanged.
A script compared the old and new versions of each section: no number or
citation was added or removed. Only the sentences changed.

## Second round, after an independent reading

A fresh reader compared the text with your thesis and the humanizer rules.
It found about forty sentences that still read as machine-written, and I
applied nearly all of the fixes. Seven of them undid mistakes from my own
first pass, for example:

- a "whereas" where the numbers actually agree;
- an "As a result" with no cause behind it;
- "Its largest value", which pointed to the wrong thing;
- a clause on the cube ("its error then decreases slowly") that added a
  claim.

Three small wording corrections change what a sentence says, so please
check them:

- **Abstract.** It said "successive layouts of the torus agree to 7e-12".
  The number is actually the difference from the finest MFS run, and the
  sentence now says so.
- **Conclusion.** It stated the almond's singular curve as a fact. It now
  says "we expect", as Section 4.2 does.
- **Introduction.** "The bound" had nothing to refer back to. It now reads
  "a lower bound for its radius", as in the abstract.

The paper is now 42 preprint pages (41 before), because the figures fall
on different pages.

## Where the style comes from

Your thesis. The profile I used is in `paper/style_profile.md`:

- the guiding "we" ("We begin by focusing on the familiar domain of the unit
  sphere", "We now present our results for the cube");
- "Note that", "Recall that", "Here," and "Specifically,";
- reasons and contrasts stated outright ("because", "since", "so that",
  "As a result", "However", "In contrast", "whereas");
- medium-to-long sentences instead of runs of short ones.

These thesis sentences are back in the text, where they are accurate:

- "The general idea of the MFS is to approximate the solution by a linear
  combination of fundamental solutions of the governing differential
  equation"
- "A well-documented disadvantage of the MFS is that this system of
  equations is often ill-conditioned"
- "Here, the conductor occupies a bounded domain ... the physical
  requirement is that the tangential component of the total electric
  field vanishes"
- "Note that p and d are orthogonal, which is essential for electromagnetic
  waves in free space"
- "indicating that the field oscillates in the y-direction"
- "the resulting system is not a simple scalar one: it is composed of 3x3
  blocks"
- "Its flat faces, combined with sharp edges and abrupt corners, produce
  discontinuities in the surface normals"
- "The cube is defined with side length 1, so that its faces lie on the
  planes x = ±1/2 ..."
- "a familiar test case in computational electromagnetics"
- "In conclusion, ..."

Your thesis's promotional phrases stay out ("highly sought after",
"daunting", "truly meshless", "exceptional capability"), as do the usual
AI words.

## Two examples

Before:
> The rule that the sources must enclose the singularities of the continued
> field, or reach them where they lie on the surface, organises all four
> layouts of this paper. For the sphere the singularity is the centre, and
> every source sphere works. ... Before we conclude, we state the limits of
> what we have shown.

After:
> One rule organises all four layouts of this paper: the sources must
> enclose the singularities of the continued field, or reach them where
> they lie on the surface. For the sphere the singularity is the centre,
> so every source sphere works. ... Before concluding, we note the limits
> of this work.

Before:
> At an error of about 1e-3 the boundary element method is the faster of
> the two in every case, about ten times at one wavelength, six times at
> two and four times at three. At the other end, the most accurate
> boundary element run ... takes a little more than two minutes in every
> case. The MFS reaches a smaller error in half of that time or less.

After:
> At an error of about 1e-3, the boundary element method is the faster of
> the two in every case: about ten times at one wavelength, six times at
> two and four times at three. However, its error then decreases slowly.
> Its most accurate run ... takes a little more than two minutes in every
> case, and the MFS reaches a smaller error in half of that time or less.

## Ideas for presenting the results (not done; your call)

You asked earlier whether the MFS clearly comes across as effective,
especially on the cube. I have not acted on this, because you asked for a
humanising pass only. The facts are already in the paper. These changes
would only make them easier to see:

1. **A headline number for every shape.** At the accuracy of Bempp's most
   accurate run, the MFS is faster on every shape and in every cube case:
   - sphere: 0.10 s against 86 s;
   - torus: 0.16 s against 71 s;
   - almond: 2.5 s against 22 s;
   - cube, standard wave: 37 s against 131 s;
   - cube, general wave: 66 s against 133 s;
   - cube, two wavelengths: 66 s against 137 s;
   - cube, three wavelengths: 53 s against 132 s.

   This could be one small table in the introduction and one sentence in
   the abstract. The numbers come from the existing data, so no new runs
   are needed.
2. **The cube table.** Table 5 now opens with "time for 1e-3", the one
   column where Bempp wins. Swapping its two halves would put Bempp's
   most accurate run against the MFS first and the 1e-3 times second.
   No number would change.
3. **The cube figures.** Bempp's best runs are drawn as about twenty
   crosses, so they dominate the panels visually. Drawing them as a plain
   line with small markers would be just as accurate and less heavy.
4. **Shortening, if wanted.** These parts could move to an appendix or be
   cut:
   - the torus image argument and the Millar/Aberra paragraph (one
     paragraph would remain);
   - the paragraph on the two-dimensional square;
   - the edge-residual figure.

   The near-field figure could also become a panel of the cost figure.
   Together this saves about three pages.
5. **The abstract is 351 words.** Elsevier abstracts are often expected to
   be about 250 words. I could not find JCP's own limit in its guide for
   authors, so please check it before submission. If it needs cutting,
   the torus conjecture and the exponent sentence could go.
