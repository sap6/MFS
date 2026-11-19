You are continuing work on the repository `sap6/mfs`, branch `ccr-116c520c-fts35f`. It holds the code and the manuscript (`paper/`) of the JCP paper "The Method of Fundamental Solutions for Maxwell's Equations", by me and my supervisor Timo Betcke; the paper grew out of my MSci thesis. Read `notes/HANDOFF.md` first, then `CLAUDE.md` and the top of `notes/NOTES.md`. They record what has been done, what was tried and failed (do not repeat it), and the rules below in more detail.

Your goal is to bring the paper and its code to a polished state, ready for submission. Work through the following in order. Do not stop until all of it is done. If you need a decision that only I can make, ask, and carry on with the rest meanwhile.

1. **Finish the open work** in `notes/HANDOFF.md`, section 3:
   - the re-run of the cube comparison (Figure 7, Table 5);
   - the numbers in the text;
   - the re-checks of earlier claims;
   - the notebooks and the PDF.

2. **Speed, fairly.** Make the MFS implementation as fast as a simple code can be, and do the same for the Bempp side.
   - "Simple" means: prefer changes that remove work or code. Examples are duplicate computations, over-resolved discretisations, and work counted in the timer that the other code does not do. Avoid changes that add machinery. A change that adds code must buy a large, measured speed-up.
   - The fairness rule is absolute. Every symmetry or speed-up goes into both codes or into neither. Each code is timed on exactly the same task (the solve plus the field at the same points), on the same machine, alternately.
   - Check that Bempp is used at its best: quadrature orders, gradings, mesh sizes, OpenCL settings, LU, and the symmetry.
   - Record every experiment as one line in the log, with its number.

3. **Code that readers of the paper will visit.** Make the repository organised, impressive and easy to understand. It should have:
   - a clear README landing page;
   - a small package with one module per shape;
   - one script per figure;
   - notebooks that run on Colab;
   - tests against exact solutions and against the stored results;
   - continuous integration (CI).

   Keep the code short and readable. Comment it in my style, as in the thesis notebooks in `thesis/`: a docstring that starts with the formula or the purpose and lists the arguments with their shapes, and short block comments. Every module must keep reproducing the stored results (`python -m pytest tests`).

4. **The paper.**
   - Make it readable and natural: humanly worded and in my voice (`notes/style_profile.md`, taken from my thesis), in a tone that conforms to academic writing.
   - It should be accessible to a Master's student. Define every concept before you use it (the MFS, parity classes, the skeleton, singular tails, the semi-normal equations, root-exponential convergence, and so on). Give the intuition behind each step. Make every assumption and every choice of parameter visible. It must still stay at the level of a JCP paper.
   - Do not change a claim without the data to support it. Every number in the text must come from `data/`, through the scripts in `figures/` or `notes/paper_numbers.py`.
   - Keep the title. Keep the abstract to at most 250 words. Follow the JCP guide for authors (`notes/jcp_checklist.md`).

5. **Figures.** Every figure must be clear and easily readable at the size it has in the paper. Check the journal layout (`paper/main_3p.tex`) as well as the preprint. That means:
   - fonts that can be read at that size;
   - the same colours and markers for the same method in every figure;
   - no overlapping legends or labels;
   - axes with units;
   - a caption and a text that explain everything the figure shows.

   Look at every page of the built PDF as an image to check.

How to work:
- Follow `CLAUDE.md`:
  - Run one heavy job at a time (4 cores, 15 GB).
  - Start long runs detached: background tasks stop after 2 hours, and the container can restart.
  - Commit and push after every step, to the branch above.
  - Never push to `sap6/mfs-paper`, which syncs to Overleaf.
  - Do not open a pull request unless I ask.
- Verify everything: the tests, `python figures/tables.py` (CI checks the tables), every figure script, and the PDF build.
- Before you finish, have a fresh subagent review the prose and, separately, the code. Apply the findings that hold up.
- I read on my phone. Send the PDF, the figures and any lists as files; only your last message is shown in full.
- At the end, give me a short list:
  - what changed;
  - what you tried and rejected, with numbers;
  - the decisions left to me (`notes/HANDOFF.md`, section 4).
