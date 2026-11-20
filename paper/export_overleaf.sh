#!/bin/sh
# A self-contained copy of the paper for Overleaf (or for the manuscript repository): main.tex, sections,
# generated tables, the figures as PDF, bibliography and the journal class, in overleaf/ and as
# mfs-paper-overleaf.zip. The figures are read from figures/ there (main.tex already lists that folder).
#     sh paper/export_overleaf.sh          (after figures/tables.py and the figure scripts)
set -e
cd "$(dirname "$0")"
if [ -f main.log ]; then
    ! grep -E "Warning:.*(undefined|multiply-defined)|Rerun to get" main.log || exit 1
fi
rm -rf overleaf mfs-paper-overleaf.zip
mkdir -p overleaf/sections overleaf/tables overleaf/figures
cp main.tex refs.bib elsarticle.cls elsarticle-num.bst highlights.txt overleaf/
if [ -f main.bbl ]; then cp main.bbl overleaf/; fi
cp sections/*.tex overleaf/sections/
cp tables/*.tex overleaf/tables/
for f in $(grep -hoE '\{fig_[a-z0-9_]+(\.pdf)?\}' main.tex sections/*.tex | tr -d '{}' | sed 's/\.pdf$//' | sort -u); do
    cp "../figures/${f}.pdf" overleaf/figures/
done
if [ -f ../figures/fig_graphical_abstract.pdf ]; then cp ../figures/fig_graphical_abstract.pdf overleaf/figures/; fi
if command -v pdflatex >/dev/null 2>&1 && command -v bibtex >/dev/null 2>&1; then
    (cd overleaf
     pdflatex -interaction=nonstopmode main.tex >/dev/null
     bibtex main >/dev/null
     pdflatex -interaction=nonstopmode main.tex >/dev/null
     pdflatex -interaction=nonstopmode main.tex >/dev/null
     if [ -f main.log ]; then ! grep -E "Warning:.*(undefined|multiply-defined)|Rerun to get" main.log || exit 1; fi
     rm -f main.aux main.log main.out main.blg main.spl main.pdf)
fi
(cd overleaf && python3 -c "import shutil; shutil.make_archive('../mfs-paper-overleaf', 'zip', '.')")
ls -lh mfs-paper-overleaf.zip
