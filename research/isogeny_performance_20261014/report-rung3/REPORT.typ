#set page(paper: "a4", margin: 16mm)
#set text(font: "Arial", size: 9pt)
#set heading(numbering: none)
= BMSS constructor optimization
2026-10-14 · Native exact-output evidence

All eight BMSS methods at 192, 224, 256 and 511 bits, degrees 7, 31 and 101: 96 exact baseline/candidate output matches. Kernel, codomain and rational map coefficients agree with the Kohel reference. Four public point transports are checked per constructor case.

== Arithmetic counts at 224 bits, degree 101
#table(columns: (1.8fr, 1fr, 1fr, 0.8fr, 0.8fr, 0.8fr), inset: 4pt,
  [Algorithm], [Baseline mul], [Candidate mul], [Quotient], [Baseline inv], [Candidate inv],
  [linear algebra], [273659435], [4640748], [58.969], [611], [209],
[stark1972], [10853630], [10856762], [1.000], [1151], [109],
[atkin1992], [547470], [548379], [0.998], [409], [109],
[atkin modcomp], [627414], [628032], [0.999], [417], [213],
[elkies1992], [20109], [20565], [0.978], [255], [106],
[elkies1998], [11127], [11277], [0.987], [249], [200],
[fast elkies], [46958], [47108], [0.997], [108], [59],
[fast elkies prime], [1073345], [906695], [1.184], [110], [110],

)
#image("COUNTS.svg", width: 100%)

The quotient counts generic Field multiplication calls. Inversions are opaque and counted separately. Setup, preparing sigma, independent reference construction and verification are excluded. Counts are deterministic for these fixed inputs; their model is not a CPU instruction count or a priced total cost.

== Changes and limits
Row-independent series products are computed once per column. The second rung batches fixed series and Newton denominators, with direct division retained at small sizes. The third uses forward elimination followed by back substitution. Raw receipts distinguish the rungs and retain any additional multiplications incurred by batching.

This shared Apple M4 Pro host has 14 logical cores and 48 GiB RAM. Wall time is unqualified; no wall samples are promoted. Qualified multi-CPU A/A and A/B comparisons, the other inventory families, named NIST whole-search/replay totals, and actual release publication remain open. The 2,000-additional-map search remains separate and unfinished.

Comparison SHA256:
#text(size: 7pt)[678d9a8974c3324fec62e74d1adc8011b63efb69379e3804ee19c6136ac5e906]
Baseline and candidate run digests and complete case inputs are in COMPARISON.json and the source-linked raw JSONL records. Read PROTOCOL.md, INVENTORY.md and ATTEMPTS.md for boundaries and censored attempts.
