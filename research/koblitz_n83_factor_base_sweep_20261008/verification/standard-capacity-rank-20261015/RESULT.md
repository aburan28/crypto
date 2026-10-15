# Conditional N83 standard-subspace capacity and rank floors

The [screen](screen.json) covers every archived standard polynomial-subspace dimension 2–18 at arities 2–6: **85 cells**. It binds six prior construction, write-once S3 replay and independent-audit receipts by SHA-256. Every input point count and sign-only column count is cross-checked across those receipts. Dimension 1 has an empty cofactor image and no usable archive. This stage makes no S3 write and launches no N83 solver.

For a base with `N` distinct nonidentity subgroup points and `K` sign-only columns, let `M = C(N+m-1,m)` be the exact number of unordered `m`-point multisets with repetition. Every multiset has one subgroup sum. If each query target has a uniform marginal over the `r-1` nonidentity primary-subgroup elements, then `E[decompositions per query] <= M/(r-1)`; multisets summing to the identity can only reduce the expectation. The union/Markov hit ceiling is `min(1,M/(r-1))`. Across `A` such queries, even if an oracle emits every valid multiset and every row is independent,

`Pr(at least K relation rows) <= min(1, A M / (K (r-1)))`.

Thus `ceil(K(r-1)/(2M))` is the first query count at which this **upper bound permits**, but does not guarantee, a 50% probability of `K` rows. The saved integer count and exact clipped rational hit ceiling are authoritative; rounded decimal means are only for display. The independent [audit](audit.json) recomputed every multiset count by a product/factorial route, checked both sides of every tight floor inequality, and verified the complete 85-cell grid. Four focused tests include a one-byte parent-receipt mutation and a rehashed semantic multiset-count corruption.

| Standard dimension | Arity | Projected points | Sign-only columns `K` | `M/(r-1)` display | First `A` permitting 50% ceiling for `K` rows | Nominal all-affine chain variables |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 15 | 6 | 32,654 | 16,327 | 0.69671582 | 11,718 | 422 |
| 16 | 5 | 64,904 | 32,452 | 0.0039702087 | 4,086,939 | 329 |
| 16 | 6 | 64,904 | 32,452 | 42.950379 | 378 | 428 |
| 17 | 5 | 130,604 | 65,302 | 0.13097993 | 249,283 | 334 |
| 17 | 6 | 130,604 | 65,302 | 2,851.1929 | 12 | 434 |
| 18 | 5 | 261,444 | 130,722 | 4.2101650 | 15,525 | 339 |
| 18 | 6 | 261,444 | 130,722 | 183,457.24 | 1 | 440 |

The largest archived base at arity four has `M/(r-1) = 0.000080516297`, so even its idealized one-hit Markov ceiling remains below 50% until **6,210** fresh uniform queries. Smaller dimensions have still lower multiset counts at the same arity. The m5/d18 and m6/d15 cases have comparable conditional row-count floors but different algebraic systems: the respective nominal all-affine chained-S3 counts are 339 variables/332 equations and 422 variables/415 equations. These model counts come from the existing source-coordinate and intermediate-x layout; they are not native build or search receipts for every cell.

The `K`-row goal here is one row per **sign-only** column. The signed-Frobenius projection map separately reduces `K` at d12, d16 and d18, so these sign-only floors must not be presented as necessary for that quotient pipeline. A query can yield zero, one or multiple valid multisets; the number found within a solver cap, their independence, and the cost of finding them are all unmeasured. The uniform-marginal premise excludes the frozen public target and any nonuniform or adaptive query schedule. These bounds cannot select the fastest factor base or establish total index-calculus runtime.

Reproduce the saved result and independent audit with:

```sh
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/standard-capacity-rank-20261015/screen.py --check
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/standard-capacity-rank-20261015/audit.py --check
python3 -m unittest discover -s research/koblitz_n83_factor_base_sweep_20261008/verification/standard-capacity-rank-20261015 -p 'test_*.py' -v
```
