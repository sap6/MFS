# MFS for Maxwell: working rules

Goal: turn the MSci thesis "Using the Method of Fundamental Solutions to Solve
Maxwell's Equations" into a JCP paper (authors: the repo owner and their
supervisor Timo). Code lives here; the manuscript lives in `sap6/mfs-paper`
(synced with Overleaf).

**Read `notes/HANDOFF.md` (21 November) and then `notes/NOTES.md` before doing anything.** It holds what works, what
failed, the current best parameters and the task list. Do not re-derive or
re-test anything recorded there.

**Layout (since 18 November 2026).** The repository is the one readers of the
paper visit: `mfs/` (the method, one module per shape), `figures/` (one
script per figure of the paper, and `tables.py`; each replots from `data/`,
`--run` recomputes), `data/`, `notebooks/` (one per shape, run on Colab),
`tests/` (`python -m pytest tests`, 20 s), `paper/`, `thesis/` (the MSci
notebooks), `notes/` (this research log). The old research scripts
(`validation/`, named throughout `notes/NOTES.md`) are no longer in the tree:
`git checkout 9f5a721 -- validation` restores them. New experiments go in a
scratch folder or `notes/`, not in the package; the package stays short.

## Rules for every session

1. One task per session, taken from the task list in `notes/NOTES.md`.
2. Prototype in 2D first (the `proto2d_*.py` scripts at commit 9f5a721 run in
   seconds); use 3D runs only to confirm.
3. Run one heavy job at a time. Antigravity VMs have **2 cores and 7.8 GiB RAM**
   (with 8 GiB swap; older containers had 4 cores / 15 GB); a 3D solve can take
   3-6 GB, and swap thrashing freezes the VM. Wrap heavy runs with `mfs-safe-run`.
4. Pipe long output through `tail`; never print full matrices or logs.
5. Before ending: add one line per experiment to the log in `notes/NOTES.md`
   (what was run, the number, the conclusion), update the task list, commit
   and push.
6. Report failures and stalls as plainly as successes. Every accuracy claim
   needs a reference (Mie, Bempp, or a known-solution check) and a number.
7. Do not write or rewrite manuscript prose unless asked. The author writes;
   Claude critiques, checks maths and references, runs experiments, formats.
8. Never push to `sap6/mfs-paper` without asking: it syncs to Overleaf.

## Standing instructions from the author (18 November 2026)

These come from the author's own messages and take precedence over rule 1
above for as long as the paper is being finished.

- "Don't stop unless I explicitly say, as I'm working on other projects
  this week." "You can always keep improving and polishing the
  paper/code/results/etc. or simplifying things that are too complex."
  Work through the work list at the top of `notes/NOTES.md`, keep it current,
  commit and push often. **On Antigravity:** never use bash `while pgrep ... sleep`
  polling loops (which stall turns); use reactive background `run_command`
  notifications, delegate multi-file audits/PDF checks to subagents so context
  stays small, and use `/goal` or `/teamwork` for multi-turn autonomous runs.
- 18 November, after a turn that ended with a report once its plan was
  finished: "Again you stopped even though I said multiple times not to
  stop." A finished plan is not a reason to end the turn: write the next
  work list at the top of `notes/NOTES.md` and continue. Send files as they
  become ready.
- "Make sure it's not something that you're not doing for BEM. The idea
  is that we're comparing the best for both." Every symmetry or speed-up
  goes into both codes or into neither.
- The title is settled ("keep the title"); do not offer alternatives.
- "I would like to keep it clear and not have tons of optimisations and
  extra lines to explain." The paper should be "simple and appropriate
  for a first year PhD student publishing their master's thesis".
- The author reads on a phone: only the last message of a turn is shown
  in full, and repository files cannot be opened there. Send PDFs,
  figures and lists as files.
- arXiv is reachable with `curl` from the shell even when the page-fetch
  tool fails.
- Writing (18 November): "make the paper more natural and use my wording from
  the thesis where it's accurate"; then "Don't rewrite the paper, just do a
  humanising pass. Use my thesis to get my style of writing." The profile is
  `notes/style_profile.md`.
- Code (18 November): "make the codebase simple and easy to read, with as few
  lines as possible (still readable and clear)", "the best way to share
  Python code in a GitHub repository accompanying a paper", "a repo that
  people visit and interact with from the paper", and "use comments in a very
  similar way to the way I used comments in my original repo" (the thesis
  notebooks in `thesis/`): docstrings that start with the formula or purpose
  and list the arguments with their shapes, indented two spaces; short block
  comments such as `# Form matrix A`, `# Solve for alpha`; aligned inline
  comments on constants and shapes. Every module must keep reproducing the
  stored results (`tests/`).

## Setup (fresh container)

```
pip install -e ".[test,bempp]"         # the package, pytest, and Bempp with its fast (OpenCL) backend
python -m pytest tests                 # 20 s warm; MFS_TEST_BEMPP=1 also checks Bempp with symmetry (16 tests)
for f in figures/*.py; do case "$f" in *style.py|*cone_exponents.py) ;; *) python "$f";; esac; done   # every figure and table from data/
apt-get install -y --no-install-recommends texlive-latex-base texlive-latex-recommended texlive-latex-extra texlive-fonts-recommended texlive-science latexmk cm-super   # only to build paper/main.pdf
```
