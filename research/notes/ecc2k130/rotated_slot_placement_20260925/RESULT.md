# n19 full-rank rotated slot placement: saturated support at every position

**Decision: FAIL the preregistered positional-support follow-up gate.** All
six positions support every one of the same `q=130873` projected subgroup
targets. Four position pairs have exactly matched labelled tuple counts `N`
and per-slot signed-column budgets `B`; each has zero support difference and
zero symmetric difference, below the frozen thresholds of 1% and 2% of q.
The experiment therefore gives no support reason to pursue position tuning
on this saturated n19, m6, full-rank `3+3+3+3+3+4` toy rung.

The [protocol](PROTOCOL.md), [input](INPUT.json), source and workflow were
frozen at `5b2b577` in draft [PR #793](https://github.com/aburan28/crypto/pull/793)
before any new six-sum result. `FROZEN.json` SHA-256 is
`ff4bd870a0c21e823fc7efa117630876118e0685da92079232c48cb2f07c5fcf`.
Both exact-head hash-only checks passed first. The subsequent single run
completed all 13 children without a timeout, nonzero exit, RSS breach or
replacement. Archive-only independent replay passed and compared every
one of the six complete subgroup multiplicity arrays, every shared target,
rank and pair decision. The raw archive SHA-256 is
`495e522064f440f0317640d28d812c5ef012f3dd12580a07e3b941c806b92d68`.

Each ordinary slot has the three normal conjugates at indices `i,i+6,i+12`;
the designated slot receives index 18. The six slots span all 19 normal
coordinates at every position. Here `B` counts per-slot distinct nonzero
signed `[4]P` pairs and `C` counts their global union; both are observed,
not assumed. `N/q` is a tuple-capacity ratio, and the Cauchy value
`N²/(M2 q)` is a proven support lower-bound ratio, with
`M2=sum_k c_k²`. Neither is an observed PDP success probability.

| Fourth-bit slot | Physical factor sizes | B / C | N / q (capacity ratio) | Exact support / q | Cauchy bound / q | M2 | Archived first-witness rank / C, same 64 Q | Sparse integer updates | Exact replay | Full S / rho |
|:--|:--|--:|--:|--:|--:|--:|--:|--:|:--|:--|
| 0 | 23×13⁵ | 41 / 41 | 8,539,739 / 130,873 = 65.252 | 130,873 / 130,873 | 98.2265% | 567,297,007 | 31 / 41 | 2,279,404 | PASS | unset / unset |
| 1 | 19×13⁵ | 39 / 39 | 7,054,567 / 130,873 = 53.904 | 130,873 / 130,873 | 97.9329% | 388,295,399 | 34 / 39 | 2,172,586 | PASS | unset / unset |
| 2 | 19×13⁵ | 39 / 39 | 7,054,567 / 130,873 = 53.904 | 130,873 / 130,873 | 97.5873% | 389,670,213 | 34 / 39 | 2,162,914 | PASS | unset / unset |
| 3 | 19×13⁵ | 39 / 39 | 7,054,567 / 130,873 = 53.904 | 130,873 / 130,873 | 97.9085% | 388,392,101 | 33 / 39 | 2,174,526 | PASS | unset / unset |
| 4 | 21×13⁵ | 40 / 40 | 7,797,153 / 130,873 = 59.578 | 130,873 / 130,873 | 97.9842% | 474,095,485 | 34 / 40 | 2,242,026 | PASS | unset / unset |
| 5 | 23×13⁵ | 41 / 41 | 8,539,739 / 130,873 = 65.252 | 130,873 / 130,873 | 98.2265% | 567,297,007 | 34 / 41 | 3,170,428 | PASS | unset / unset |

The equal-budget pairs are `(0,5)` at `N=8,539,739, B=41`, and
`(1,2)`, `(1,3)`, `(2,3)` at `N=7,054,567, B=39`. Every pair has exact
support intersection 130,873, symmetric difference zero, and zero of the
64 frozen target labels flipping membership. Position 4 has no equal-budget
partner and cannot isolate placement from factor-size changes. All 15 pairs,
including unmatched pairs, have identical complete *support sets* (the
whole subgroup); their multiplicity arrays are generally different. Even
positions 0 and 5 have different array hashes although their second moments
are equal. The archived first-witness rank differs within two equal-budget
pairs, but the chosen 64 rows depend on a deterministic prior witness
selection in a complete oracle. This is only a toy archive-row diagnostic,
not measured independent relation acquisition or n131 matrix rank.

The exact second moments prove strong lower bounds here (97.59–98.23% of q),
but the observed support is already q in all arms. Their small differences
therefore do not change the target-membership decision. The earlier
rank-two beta-3 corpus [#767](https://github.com/aburan28/crypto/pull/767)
had 117,649 labelled tuples and 62,389 supported targets; its 327,325
second moment was already reported in [#769](https://github.com/aburan28/crypto/pull/769).
That is a different tuple/column budget and is shown only as the historical
unsaturated control, never as an equal-budget position arm.

The sequential charged run took 72.947 s wall from the first child to the
final receipt, including group setup, all q-point tables, factor lifts,
six exact convolutions, 64 selected-target row checks per arm, serialization,
independent verification and pair analysis. Producer children used
133,465–133,481 point additions, 217–221 scalar calls and
2,162,914–3,170,428 sparse integer convolution updates per arm. Independent
verifiers used 131,895–131,911 point additions per arm and separately
recomputed all q multiplicities. Child walls were 2.297–2.877 s producer
and 8.645–10.805 s verifier; the maximum recorded child RSS was 75,628,544
bytes, below 512 MiB. Field-operation types and individual child costs are
preserved in the raw result/verify JSON and [receipt](evidence/receipt.json).
The integer update counts are versioned as explicit *post-outcome cost
reconstruction* in [accounting.json](evidence/accounting.json); they do not
alter the frozen decision.

The rational next positional test, if pursued, is an unsaturated n19 frame:
retain the six rank-two slots of #767 and assign one predeclared unused
normal coordinate to each slot in turn, comparing only equal-`N`, equal-`B`
pairs. That would test whether placement affects actual support while
avoiding the present all-q ceiling. It requires a separate pre-outcome PR.
The present experiment measured no solver, negative PDP difficulty,
relation-collection cost, challenge logarithm, n131 support transfer,
matched rho or end-to-end index-calculus `S`.
