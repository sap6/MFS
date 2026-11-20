# JCP submission checklist (20 November 2026)

Checked against the Elsevier guide for authors of the Journal of
Computational Physics. The guide itself could not be downloaded from this
machine because sciencedirect.com is blocked, so these rules come from
search results that quote it. Please check them once against the live
guide.

| Rule | Ours | Status |
|---|---|---|
| Abstract at most 250 words | 237 words | done |
| Highlights: 3 to 5 bullets, at most 85 characters each, including spaces (encouraged, not required) | 5 bullets, 71 to 76 characters (`73` to `78` with `- `, `paper/highlights.txt`) | done |
| Keywords: 1 to 7, in English; avoid multi-word terms joined by "and" or "of" | 6 keywords | "Method of Fundamental Solutions" contains "of", but it is the name of the method, so I would keep it |
| Declaration of competing interest | missing | **authors**: template in `main.tex` (commented) |
| CRediT authorship statement (required from the corresponding author) | missing | **authors**: template with suggested roles in `main.tex` |
| Declaration of generative AI use | missing | **authors**: Elsevier template in `main.tex`; the work used AI tools for code, experiments and editing |
| Corresponding author marked on the title page | not marked | **authors**: add `\corref{cor1}` and `\cortext` |
| Acknowledgements, funding | missing | **authors**: template in `main.tex` (includes MSci project report citation) |
| Data and code availability | "Code and data availability" section with the GitHub link | done; add the Zenodo DOI once a release is archived |
| Graphical abstract (encouraged, optional) | `figures/fig_graphical_abstract.pdf` / `.png` (`figures/graphical_abstract.py`) | done |
| Cover letter | `notes/cover_letter_jcp.md` | done (fill `[Submission Date]`) |
| Reference style: numbered (elsarticle-num) | elsarticle-num | done |
| Length | no page limit known; 28 pages in the journal `3p,times` layout (49 preprint pages), in line with the JCP papers we cite | fine |

The JCP home page also says the journal "strongly encourages authors to
make available all software or data that would allow published results to
be reproduced". The new repository layout is built for that.
