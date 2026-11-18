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
