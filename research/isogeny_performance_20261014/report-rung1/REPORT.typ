#set page(paper: "a4", margin: 16mm)
#set text(font: "Arial", size: 9pt)
#set heading(numbering: none)
= BMSS constructor optimization
2026-10-14 · Native exact-output evidence

All eight BMSS methods at 192, 224, 256 and 511 bits, degrees 7, 31 and 101: 96 exact baseline/candidate output matches. Kernel, codomain and rational map coefficients agree with the Kohel reference. Four public point transports are checked per constructor case.

== Arithmetic counts at 224 bits, degree 101
#table(columns: (1.8fr, 1fr, 1fr, 0.8fr, 0.8fr, 0.8fr), inset: 4pt,
  [Algorithm], [Baseline mul], [Candidate mul], [Quotient], [Baseline inv], [Candidate inv],
  [linear algebra], [273659435], [5639435], [48.526], [611], [611],
[stark1972], [10853630], [10853630], [1.000], [1151], [1151],
[atkin1992], [547470], [547470], [1.000], [409], [409],
[atkin modcomp], [627414], [627414], [1.000], [417], [417],
[elkies1992], [20109], [20109], [1.000], [255], [255],
[elkies1998], [11127], [11127], [1.000], [249], [249],
[fast elkies], [46958], [46958], [1.000], [108], [108],
[fast elkies prime], [1073345], [1073345], [1.000], [110], [110],

)
#image("COUNTS.svg", width: 100%)

The quotient counts generic Field multiplication calls. Inversions are opaque and counted separately. Setup, preparing sigma, independent reference construction and verification are excluded. Counts have no sampling error for these fixed inputs; their model is not a CPU instruction count or a priced total cost.

== Changes and limits
Row-independent series products are computed once per column. The later rung batches fixed series and Newton denominators, with direct division retained at small sizes. Raw receipts distinguish the rungs and retain any additional multiplications incurred by batching.

This shared Apple M4 Pro host has 14 logical cores and 48 GiB RAM. Wall time is unqualified; no wall samples are promoted. Qualified multi-CPU A/A and A/B comparisons, the other inventory families, named NIST whole-search/replay totals, and actual release publication remain open. The 2,000-additional-map search remains separate and unfinished.

Comparison SHA256:
#text(size: 7pt)[7f5cadad4bc607b456d15fcf60c4363a9243ab28732e44546ad309d883ba42bf]
Baseline and candidate run digests and complete case inputs are in COMPARISON.json and the source-linked raw JSONL records. Read PROTOCOL.md, INVENTORY.md and ATTEMPTS.md for boundaries and censored attempts.
