# Draft of the paper

First draft (17 November 2026) of the journal version of the MSci thesis.
It lives here, next to the code that produces its numbers, and is not yet in
`sap6/mfs-paper` (the Overleaf copy).

    python3 make_tables.py                       # table bodies from ../validation/*.npz into tables/
    pdflatex main && bibtex main && pdflatex main && pdflatex main

- `main.tex`, `sections/*.tex`: the text. Lines starting with `%% TODO` are
  open points for the authors.
- `tables/*.tex`: generated; do not edit by hand. `make_tables.py` also
  prints the numbers that are quoted in the running text.
- `refs.bib`: every entry with a DOI was compared with its Crossref record
  (`check_bib.py`, report in `check_bib_report.txt`); arXiv entries with
  `check_arxiv.py`. Entries marked "unverified" in a note still need a
  bibliographic record.
- Figures are read from `../validation/` (`fig_*.pdf`); each has a script of
  the same name there.
- `elsarticle.cls`, `elsarticle-num.bst`: copies of the journal class.

`main.pdf` is the preprint layout (12 pt, 42 pages). `main_3p.pdf` is the same
text in the journal's layout (26 pages): change the first line of `main.tex`
to `\documentclass[final,3p,times]{elsarticle}` to rebuild it.

LaTeX packages: see the `apt-get` line in `../CLAUDE.md`.
