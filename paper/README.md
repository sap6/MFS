# Draft of the paper

First draft (17 November 2026) of the journal version of the MSci thesis.
It lives here, next to the code that produces its numbers, and is not yet in
`sap6/mfs-paper` (the Overleaf copy).

    python3 make_tables.py                       # table bodies from ../validation/*.npz into tables/
    pdflatex main && bibtex main && pdflatex main && pdflatex main
    pdflatex supplement && bibtex supplement && pdflatex supplement && pdflatex supplement

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

`main.pdf` is the short paper (second draft, 12 pt preprint layout, 25 pages)
and `supplement.pdf` its supplementary material (11 pages): derivations,
settings and the tables that the text only summarises. `highlights.txt` holds
the highlights. The first draft (43 pages) is kept in `long/`; it is
superseded and builds there with `TEXINPUTS=..: pdflatex main`.

LaTeX packages: see the `apt-get` line in `../CLAUDE.md`.
