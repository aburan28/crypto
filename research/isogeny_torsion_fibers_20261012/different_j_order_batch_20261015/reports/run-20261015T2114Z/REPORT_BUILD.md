# Visual render receipt

The native `different_j_order_batch_summary` binary read the frozen
run, meters, independent replay and calibration, then generated
`SUMMARY_V1.json` and `GRAPH_V1.svg`. The SVG was embedded directly in
both Typst report sources. The 150 dpi graph PNG was rendered from that
same SVG through `GRAPH_RENDER_V1.typ`:

```sh
cd research/isogeny_torsion_fibers_20261012/different_j_order_batch_20261015
typst compile --format png --ppi 150 GRAPH_RENDER_V1.typ GRAPH_V1.png
typst compile --root . reports/run-20261015T2114Z/REPORT.typ reports/run-20261015T2114Z/REPORT.pdf
typst compile RESULTS.typ RESULTS.pdf
```

The `cd` starts from the repository root. Both PDFs
were checked with `pdfinfo` (two A4 pages each), rendered with
`pdftoppm`, and visually inspected page by page. The graph labels,
four target counts, min–max whiskers, table values, source IDs and
limitations are visible. The `benchmarks/` PDF, SVG and PNG copies
were compared byte-for-byte with their study originals.
