#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
= A16: the folded table's scan arrays collapsed onto 2 MiB pages

*Exploration, no round.* One binary (`ic-a16-v8-e7136eb8`, `huge-pages.patch` on v8) with three arms set by `KIC_HUGE_PAGES`:
- `off`: v8's pages;
- `filter`: the filter alone collapsed;
- `all`: the filter, the bucket starts and the stored keys collapsed.

Collapsing is `MADV_HUGEPAGE`, then `MADV_COLLAPSE`, after the build. The runs: 2026-10-17 04:24–04:54 UTC on the Emerald Rapids host R12 ran on; `M1`'s two rows at the three target sizes; six rounds with the order rotating; 108 isolated processes, all completed.
- *Contended runs:* one (`off`, 2^44.3, round 1), excluded.
- *Logarithms:* every pair recovered the same logarithm in all three arms.

*Outcome:*
- *The filter alone* gains 4.6% cold at 2^47.2 ([1.026, 1.067]) and is inconclusive at 2^44.3 and 2^44.5.
- *Every array* gains nothing: the build pays for collapsing the stored keys.
- *Why the filter's gain is partial:* only the whole 2 MiB pages inside the filter's block collapse, and the block starts off a boundary.

#figure(image("cold-ratio.svg", width: 100%), caption: [Cold time, `off` over each arm, per pair (open: a contended run), with the clean pairs' geometric mean and 95% t interval. Dashed: 1.03.])

#table(columns: 5, stroke: 0.4pt, inset: 3pt,
  [*curve*], [*log₂ r*], [*off / filter [95%]*], [*off / all [95%]*], [*collection off / filter*],
  [`icv1-f2m53-tm56619371-dac20a85`], [44.3], [1.017 [0.968, 1.068]], [0.970 [0.910, 1.034]], [1.027],
  [`icv1-f2m59-tm943548413-98844ecc`], [44.5], [0.981 [0.936, 1.029]], [0.979 [0.941, 1.020]], [1.018],
  [`icv1-f2m61-t158598901-ab42b6c5`], [47.2], [1.046 [1.026, 1.067]], [1.012 [0.976, 1.049]], [1.062],
)

*Smoke records* (single processes at 2^47.2, before the timed runs; observations, not measurements):
- glibc's `malloc.hugetlb=1` on every allocation made collection faster (1395 → 1177 ms) but the build slower (213 → 565 ms).
- A16's full collapse gave cold 1626 → 1558 ms.

Scope: one host, one core, three sizes, two rows a size. THP is `madvise` here.
