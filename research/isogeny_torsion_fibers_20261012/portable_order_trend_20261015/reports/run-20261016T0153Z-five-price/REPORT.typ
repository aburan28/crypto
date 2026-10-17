#set page(paper: "a4", margin: 16mm)
#set text(size: 10pt)
= Five-order missing-price sensitivity run
*2026-10-16 01:53 UTC · Darwin arm64 shared host.* The source and
binary were frozen before execution; their SHA-256 hashes and the
five input run hashes are in `sensitivity_five/RUN_SPEC_V1.json`.
The native derivation checked 420 child records and 8,925 scalar
outputs against independent replay receipts, and all group totals
against the five-order trend. Its first 48 rows and four older N64
summaries were unchanged from V1. Sixty matched rows were emitted.

The graph plots N64 median extra omitted rho group-addition equivalents
per observed rho step required to tie the full-image counted lower
bound, setting all omitted full-image work to zero. At order 1,201 the
median was 1,390.865, with three-round range 1,266.482–1,575.602.
This is a price threshold, not a measured primitive price or complete
portable cost. The source-linked `REPORT.md` and `sensitivity_five/RUN.md`
contain commands, hashes, values, and limits.

#image("GRAPH.svg", width: 100%)
