# Conditional relation-query floors for retained N83 pilot pairs

The [saved aggregate receipt](rank-floor.json) derives a necessary
query-count floor from every one of the **2,232** admissible pair/shape
support counts in the frozen [archive-backed support screen](../pilot-pair-support-screen-20261014/RESULT.md).
It hashes the complete per-pair floor table without copying the 1.4-MiB
parent receipt. This is exact integer arithmetic conditional on the parent
archive evidence, not a new N83 solver experiment.

Let `K` be a *specified goal* of `K` independent relation rows, here set to
the small base's signed-Frobenius orbit-column count (64, 256 or 600).
For one query, the parent screen bounds the expected number of point-multiset
decompositions by `M/(r-1)` when the query is uniformly distributed over the
nonidentity primary subgroup. Across `A` queries, each with that uniform
marginal, linearity of expectation gives at most `A M/(r-1)` decompositions;
independence between queries is unnecessary. Each decomposition can supply
at most one independent folded row, and any large-prime graph closure can
combine rows but cannot create rank beyond the input-row count. Markov then
gives

`Pr(rank from this relation stream >= K) <= min(1, A M / ((r-1) K))`.

Thus `ceil(K (r-1)/(2 M))` queries are **necessary before this upper bound
even permits** a 50% chance of `K` rows. The floor is not sufficient for a
hit, graph closure, rank, or an individual logarithm. It assumes fixed
archived domains, one specified shape and uniform query marginals; a
nonuniform or target-adaptive schedule, another relation source, or a
different one-target method is outside its scope. The selected public target
and the worker's deterministic first trial are not covered.

For the five-summand **3+2** shape, the floors over all current-worker
admissible retained pairs are:

| Desired rows `K` | Ordered pairs | Lowest floor | Median floor | Highest floor |
| ---: | ---: | ---: | ---: | ---: |
| 64 | 206 | 78,029 | 428,616 | 7,001,510,487 |
| 256 | 186 | 4,878 | 26,795 | 437,687,082 |
| 600 | 166 | **889** | 4,878 | 79,681,727 |

The same-base K600 sequential-seed-`2026100802` contrast from the parent
screen has a **33,975-query** floor with the Gray-prefix pool's 97 residual
orbits, versus **889 queries** with the hash pool's 600 residual orbits. Pool
policy differs, so this is a support comparison rather than a causal timing
experiment. In the five-summand 4+1 shape, even the most favorable K600
pair has a 1,776-query floor. The six-summand shapes can have floors of one
for some K600 pairs, but their SAT search time and natural relation yield
have not been established. These contrasts prioritize future complete-cost
experiments; they do not select a factor base.

The receipt also records all twelve shape/`K` summary panels, medians,
extremes, counts above 100 and 1,000 queries, two exact selected-pair rows,
and a SHA-256 over the full sorted 2,232-row floor table. Its [replay
program](rank_floor.py) pins the parent support JSON by SHA-256, independently
recomputes every multiset count, checks the tight integer ceiling boundary,
and fails closed on a changed parent. Reproduce with:

```sh
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/pilot-pair-rank-floor-20261014/rank_floor.py
python3 -m unittest discover -s research/koblitz_n83_factor_base_sweep_20261008/verification/pilot-pair-rank-floor-20261014 -p 'test_*.py' -v
```

The focused tests check a tight small-order Markov threshold, all archived
shape populations and a modified-parent rejection. No S3 object was read or
written, no N83 solver call occurred, and the conservative search ledger
remains **22,200 seconds**. In particular, no per-query time was measured
here, so these query floors must not be multiplied by a guessed solver wall
time to claim a cold-runtime bound or a comparison to rho.

Validation: the frozen receipt replayed byte-for-byte; three new focused
tests, four parent-support tests, the 63-test N83 study suite and the 43-test
boundary autolab suite passed. This extension changes Python analysis and
documentation only. The current release-library rerun remains [recorded as
interrupted](../known-witness-stage-preflight-20261014/RESULT.md); these
source-only checks do not represent a fresh full-library pass.
