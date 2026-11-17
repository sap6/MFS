# MFS for Maxwell: working rules

Goal: turn the MSci thesis "Using the Method of Fundamental Solutions to Solve
Maxwell's Equations" into a JCP paper (authors: the repo owner and their
supervisor Timo). Code lives here; the manuscript lives in `sap6/mfs-paper`
(synced with Overleaf).

**Read `NOTES.md` before doing anything.** It holds what works, what failed,
the current best parameters and the task list. Do not re-derive or re-test
anything recorded there.

## Rules for every session

1. One task per session, taken from the task list in `NOTES.md`.
2. Prototype in 2D first (`validation/proto2d_*.py` run in seconds); use 3D
   runs only to confirm.
3. Run one heavy job at a time. The container has 4 cores and 15 GB shared by
   all processes; a 3D solve can take 4-7 GB.
4. Pipe long output through `tail`; never print full matrices or logs.
5. Before ending: add one line per experiment to the log in `NOTES.md`
   (what was run, the number, the conclusion), update the task list, commit
   and push.
6. Report failures and stalls as plainly as successes. Every accuracy claim
   needs a reference (Mie, Bempp, or a known-solution check) and a number.
7. Do not write or rewrite manuscript prose unless asked. The author writes;
   Claude critiques, checks maths and references, runs experiments, formats.
8. Never push to `sap6/mfs-paper` without asking: it syncs to Overleaf.

## Setup (fresh container)

```
pip install numpy scipy matplotlib numba bempp-cl gmsh meshio pyopencl pocl-binary-distribution   # the last two give Bempp its fast (OpenCL) backend
apt-get install -y libglu1-mesa libxcursor1 libxinerama1 libxft2 libxrender1 libgl1   # for gmsh
cd validation && python3 mie.py        # self-check, should print errors ~1e-14
apt-get install -y --no-install-recommends texlive-latex-base texlive-latex-recommended texlive-latex-extra texlive-fonts-recommended texlive-science latexmk cm-super   # only to build paper/main.pdf
```
