# Rank-query necessary condition for the retained N83 bases

This screen applies to a **one-row full-smooth oracle**: each query asks for a decomposition of one subgroup point into `m` points of a fixed base, and returns at most one relation row. Its query targets must each be uniform over the whole prime-order subgroup. Dependence between queries is allowed. The current capped cold driver has not demonstrated this sampling model for the 81-bit primary arm, so these numbers are a guide for a future full-width sampler, not measurements of that driver.

The generic Koblitz IC sampler now draws a uniform full-width additive scalar for subgroup orders above 64 bits, which supplies this target distribution for its `[a]G+[b]Q` probes when `b` is nonzero. The retained-base cold driver is still diagnostic `a=1` only. Applying this query screen to an actual primary cold run additionally requires a verified adapter from the stored base and a one-row full-smooth oracle with the stated semantics.

For a base of `B` distinct points and subgroup order `r`, exactly `M = binomial(B+m-1,m)` unordered multisets with repetition exist. Each has one group sum. A uniform query hits a full-smooth sum with probability at most `p = min(1,M/r)`, including possible identity sums. If `S` of `q` queries return a row, then `E[S] <= qp` without any independence assumption. Matrix rank is at most `S`, so for `K` required independent rows:

`Pr(rank >= K) <= 0` for `q < K`; otherwise `Pr(rank >= K) <= min(1, qp/K)`.

Consequently, a necessary number of queries before this upper bound can reach one half is `max(K, ceil(Kr/(2 min(M,r))))`. This is a **lower bound on a query count needed even to permit a 50% completion probability under the bound**. It does not assert that the probability actually reaches 50%, that rows are independent, or that a solver finds an existing decomposition. It is not a lower bound on every index-calculus algorithm: a solver returning multiple independent rows for one query, a large-prime graph, or a different accepted relation event needs a separate analysis. Biased or deliberately chosen query targets also fall outside the stated uniform-target model.

The table gives exact necessary query counts for the retained orbit-column sizes. Each column corresponds to `166K` distinct signed-Frobenius points. The 81-bit `a=0` arm is the primary workload; the 53-bit `a=1` arm is diagnostic.

| Arm | K required rank rows | m=4 | m=5 | m=6 |
| --- | ---: | ---: | ---: | ---: |
| primary a=0 | 64 | 145,677,800,576 | 68,534,909 | 38,688 |
| primary a=0 | 256 | 2,277,179,841 | 267,904 | 256 |
| primary a=0 | 600 | 176,888,106 | 8,880 | 600 |
| diagnostic a=1 | 64 | 517 | 64 | 64 |
| diagnostic a=1 | 256 | 256 | 256 | 256 |
| diagnostic a=1 | 600 | 600 | 600 | 600 |

At K=600 on the primary arm, the four-summand bound alone requires more than 176 million uniform one-row queries before a half-probability full-rank outcome is even compatible with this first-moment bound. The six-summand support bound is vacuous at this size; the query condition then reduces to the obvious `q >= K` and says nothing about the cost of solving a six-summand equation. These facts prioritize implementing and measuring higher-arity or partial-relation routes, while preserving four-summand rows as diagnostics. They cannot select a total-runtime winner.

`pilot-01/rank-query-screen.json` contains the exact integers and fractions for all 30 `(arm,K,m)` cases. `rank_query_screen.py` reads the immutable support-moments receipt, verifies its panel hash and replay binding, recomputes the multiset counts, and writes this new immutable receipt. `test_rank_query_screen.py` checks the threshold boundary in all cases and exhaustively counts a small-group example, including correlated queries. Reproduce with:

```sh
python3 research/koblitz_n83_factor_base_sweep_20261008/rank_query_screen.py
python3 -m unittest discover -s research/koblitz_n83_factor_base_sweep_20261008 -p test_rank_query_screen.py -v
```
