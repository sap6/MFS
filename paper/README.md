# The paper

The manuscript, in LaTeX (Elsevier's `elsarticle`, preprint layout). In `paper/`:

    python ../figures/tables.py      # table bodies from ../data into tables/ (do not edit those by hand)
    latexmk -pdf main.tex            # or: pdflatex main && bibtex main && pdflatex main && pdflatex main
    sh export_overleaf.sh            # a self-contained copy for Overleaf: overleaf/ and mfs-paper-overleaf.zip

- `main.tex`, `sections/*.tex`: the text; `refs.bib`: the references; `highlights.txt`.
- The figures are read from `../figures/`, where each has a script.
- `main_3p.tex` is the same paper in the journal's single-column `3p` layout, to count pages
  (`sed 's/\[preprint,12pt\]/[final,3p,times]/' main.tex > main_3p.tex && latexmk -pdf main_3p.tex`).
- Before submission the authors add the declarations that JCP requires (templates are
  commented in `main.tex`; see `../notes/jcp_checklist.md`).
