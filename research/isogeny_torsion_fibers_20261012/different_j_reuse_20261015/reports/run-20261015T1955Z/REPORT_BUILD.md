# Report rendering receipt

The first `typst compile REPORT.typ REPORT.pdf` attempt failed because
`../../GRAPH_V3.svg` was outside Typst's default project root:
`path "../../GRAPH_V3.svg" would escape the project root`.
No source or graph was replaced. The successful command was:

```sh
typst compile --root /Volumes/SSD990-2/.codex/worktrees/d92d/crypto/research/isogeny_torsion_fibers_20261012/different_j_reuse_20261015 REPORT.typ REPORT.pdf
```

`pdfinfo` reports two A4 pages. Both were rendered to PNG and
visually inspected; the table, two quantitative graphs and limits are
fully visible. The SVGs and their 150 dpi PNG renders were generated
from the same native graph outputs and Typst render wrappers in the
study directory.
