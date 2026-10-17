#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
= A17: the presence filter on huge pages from its first touch

*Exploration, no round.* v8 (`ic-r12-v7-23138d20`) against v8 with the filter
copied into a fresh 2 MiB-aligned anonymous mapping, advised `MADV_HUGEPAGE`
before it is written (`ic-a17-v8-fb4284d8`, `filter-huge-pages.patch`). A
third arm, `heap`, is the same binary with `KIC_FILTER_HUGE=0` (the copy alone).

*The timed runs:* 2026-10-17 05:08–05:35 UTC; 108 isolated processes, all completed. Four were contended and are excluded. Every logarithm is the same in all three arms.

*The stage probes:* R04's probes on both arms, interleaved; 48 processes.

*Outcome:* a gain at 2^47.2 only, so no round. A round needs every target size above 1.03.

#figure(image("cold-ratio.svg", width: 100%), caption: [Cold-time ratios per pair (open: a contended run) and the clean pairs' geometric mean with its 95% t interval. Dashed: 1.03.])

#table(columns: 6, stroke: 0.4pt, inset: 3pt,
  [*curve*], [*log₂ r*], [*v8 / a17 [95%]*], [*v8 / heap*], [*heap / a17*], [*filter stage, v8 → a17 (ns)*],
  [`icv1-f2m53-tm56619371-dac20a85`], [44.3], [1.047 [0.999, 1.098]], [1.014], [1.023], [4.75 → 4.78],
  [`icv1-f2m59-tm943548413-98844ecc`], [44.5], [0.976 [0.944, 1.009]], [0.981], [0.994], [4.92 → 4.56],
  [`icv1-f2m61-t158598901-ab42b6c5`], [47.2], [1.051 [1.029, 1.075]], [0.991], [1.054 [1.026, 1.082]], [6.43 → 5.32],
)

*Reading:*
- *The copy costs nothing measurable,* and the gain at 2^47.2 is the pages'.
- *The gain follows the filter's size.* The 8 MiB filter spans 2,048 base pages, as many as this core's second-level TLB holds. The 4 MiB filters span 1,024 base pages and show no consistent change.
- *This is a hypothesis.* The VM exposes no performance counters, so TLB misses were not counted.
- *Huge pages were confirmed* through `AnonHugePages` on unisolated diagnostic runs: 4,096 kB and 8,192 kB.

Scope: one host, one core, three sizes, two rows a size. THP is `madvise` here.
