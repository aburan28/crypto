# Historical factor-base evaluation

This is a browsing copy of the [study](../../research/factor_base_evaluation_20261015/README.md).
The [combined run/results report](REPORT.pdf) is the unchanged historical PDF;
its source HTML, Sage/C source snapshots, raw evidence and original commands are
linked from the study. It is not a repository-native performance claim.

The SVG and 150 dpi PNG pairs are rendered from the same committed summary by
[render_figures.py](../../research/factor_base_evaluation_20261015/reports/render_figures.py).
Run that rendering-only script, then copy the four rendered files here. Preserve
the original PDF and its frozen figures. No new statistical analysis is performed.

Historical stage CPU rates at 24 bits: low-x 1,352, coset 1,167, random 1,144,
low digit weight 859 independent rows/CPU second; coset/random 1.020
(95% field-bootstrap interval 0.959–1.083). No full-DLP or rho ratio is available.
The host was Xeon Platinum 8573C, pinned CPU 0; no repository isolation receipt
or A/A measurement exists, so timings remain descriptive and exploratory.
