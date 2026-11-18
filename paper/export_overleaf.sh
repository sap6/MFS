#!/bin/sh
# A self-contained copy of the paper for Overleaf (or for the manuscript repository): main.tex, sections,
# generated tables, the figures as PDF, bibliography and the journal class, in overleaf/ and as
# mfs-paper-overleaf.zip. The figures are read from figures/ there (main.tex already lists that folder).
#     sh export_overleaf.sh          (run in paper/ after make_tables.py and the figure scripts)
set -e
rm -rf overleaf mfs-paper-overleaf.zip
mkdir -p overleaf/sections overleaf/tables overleaf/figures
cp main.tex refs.bib elsarticle.cls elsarticle-num.bst highlights.txt overleaf/
cp sections/*.tex overleaf/sections/
cp tables/*.tex overleaf/tables/
for f in $(grep -ho "{fig_[a-z_]*\.pdf}" sections/*.tex | tr -d '{}' | sort -u); do cp ../validation/$f overleaf/figures/; done
(cd overleaf && pdflatex -interaction=nonstopmode main.tex > /dev/null && bibtex main > /dev/null && pdflatex -interaction=nonstopmode main.tex > /dev/null && pdflatex -interaction=nonstopmode main.tex > /dev/null \
  && grep -c "undefined" main.log | sed 's/^/undefined references in the exported copy: /' && rm -f main.aux main.log main.out main.blg main.spl main.pdf)
(cd overleaf && python3 -c "import shutil; shutil.make_archive('../mfs-paper-overleaf', 'zip', '.')")
ls -la mfs-paper-overleaf.zip | cut -c20-
