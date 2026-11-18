# Draft of the paper

Third draft (17 November 2026, timings and comparison of 18 November) of the journal version of the MSci thesis.
It lives here, next to the code that produces its numbers, and is not yet in
`sap6/mfs-paper` (the Overleaf copy).

    python3 make_tables.py                       # table bodies from ../validation/*.npz into tables/
    latexmk -pdf main.tex                        # or: pdflatex main && bibtex main && pdflatex main && pdflatex main

- `main.tex`, `sections/*.tex`: the text. One paper, no supplement: the
  derivations are in the text and three appendices (`sections/appendix.tex`)
  hold the cone exponents, the line quadrature and the Bempp settings.
  Lines starting with `%% TODO` are open points for the authors. The title
  is settled (the author, 18 November).
- Length: the target is 22-24 pages in the journal's layout (at present 23, and 38 in the
  12 pt preprint layout). To count them,
  build `main.tex` with `\documentclass[final,3p,times]{elsarticle}`
  (`sed 's/\[preprint,12pt\]/[final,3p,times]/' main.tex > main_3p.tex`).
- `tables/*.tex`: generated; do not edit by hand. `make_tables.py` also
  prints the numbers that are quoted in the running text.
- Timings: the text quotes few seconds; they are in the tables. To measure
  all of them again on one machine: `sh ../validation/time_paper.sh`
  (about three hours), then the figure scripts and `make_tables.py`.
- `refs.bib`: every entry with a DOI was compared with its Crossref record
  (`check_bib.py`, report in `check_bib_report.txt`); arXiv entries with
  `check_arxiv.py`.
- Figures are read from `../validation/` (`fig_*.pdf`); each has a script of
  the same name there.
- `elsarticle.cls`, `elsarticle-num.bst`: copies of the journal class.
- `highlights.txt`: the highlights.

Earlier versions: the first draft (43 pages) is in `long/` and builds there
with `TEXINPUTS=..: pdflatex main`; the second draft (25 pages plus an
11-page supplement) is in the git history (commit 485b981).

LaTeX packages: see the `apt-get` line in `../CLAUDE.md`.

## After a new timing run

`sh ../validation/time_paper.sh` (three hours, nothing else on the machine),
then `python3 ../validation/fig_cost_all.py`, `fig_cube_cases.py`,
`fig_cube_convergence.py`, `python3 make_tables.py`, `latexmk -pdf main.tex`.
Tables and figures follow by themselves. These sentences quote times or
ratios and have to be compared with the output of `make_tables.py` (lines
starting with `cube_ends`) and of `python3 ../validation/crossovers.py`:

- `sections/cost.tex`, "Smooth bodies and the almond": "in about five
  seconds or less", "after about twenty seconds", the sentence on 3.5 GHz
  (`../validation/almond_cost_freq_3.5GHz.npz`), "three to ten times
  faster" and "ten to twenty times faster" (symmetry gains, measured by
  hand: NOTES.md item 83).
- `sections/cost.tex`, "The cube": "after three quarters of a minute";
  "about ten times at one wavelength, six times at two and four times at
  three"; the levels J = 9, 8, 7; "within 17%"; the three runs of the
  standard wave (24, 20, 30 seconds); "level near 1.5e-4", "1.5e-4, 3e-4,
  5e-4"; "in less than a minute"; the one-depth sentence (8, 11, 11 and 30
  seconds).
- `sections/cost.tex`, "What the comparison shows": "below about 5e-3"
  (almond).
- `sections/cube.tex`: "between a third and a half of the time" (one depth
  per octave); "about one and a half times ... about twice" (general wave).
- `sections/almond.tex`: "takes about two minutes".
- `sections/appendix.tex`, Appendix C: "five to ten times", "three to seven
  times", "ten to twenty times".

## A copy for Overleaf

`sh export_overleaf.sh` writes `overleaf/` and `mfs-paper-overleaf.zip`: a
self-contained copy (text, tables, figures, bibliography, journal class)
that compiles on its own. Upload the zip to Overleaf, or copy the folder
into `sap6/mfs-paper`. Nothing is pushed there by the scripts.
