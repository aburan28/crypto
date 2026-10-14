# Index-calculus programme paper: figures (no new measurement)

Browsing copies of the figures of
[`papers/ecdlp-index-calculus-program/`](../../papers/ecdlp-index-calculus-program/)
([paper PDF](paper.pdf)), the programme record of the repository's ECDLP
index-calculus work, 2026-07-01 to 2026-10-14. The paper directory is the
source of truth: every plotted value and its source artefact are in
[`figures/data.json`](../../papers/ecdlp-index-calculus-program/figures/data.json),
rendered by `figures/make_figures.py` (SVG plus a 150 dpi PNG from the same
source). No experiment was run for this entry.

| figure | what it shows |
|:--|:--|
| [`fig-timeline`](fig-timeline.svg) | commits per ISO week on `main`, with keyword-matched themes |
| [`fig-rootfind-oracle`](fig-rootfind-oracle.svg) | root finding vs Gröbner solvers: oracle stage and complete solves |
| [`fig-k13-progression`](fig-k13-progression.svg) | Koblitz n = 13 root-finding IC arms across four ecbench loops, against rho and F4 |
| [`fig-cold-gap`](fig-cold-gap.svg) | cold one-target IC/rho ratios; every arm above parity |
| [`fig-koblitz-online`](fig-koblitz-online.svg) | one-target online rho/IC ratios, exploratory vs admitted |
| [`fig-prime`](fig-prime.svg) | FGHR D4 vs Gaudry S4 at k = 4; PKM tower IC/rho by size |
| [`fig-engine`](fig-engine.svg) | selected Gröbner-engine ratios |
