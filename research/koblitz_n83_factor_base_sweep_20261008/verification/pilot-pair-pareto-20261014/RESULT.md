# Source-derived frontier of retained N83 ordered pairs

The [receipt](receipt.json) replays all 648 ordered pairs in the pinned pilot
support census against the retained object manifest and the current worker's
3+2 capacity formula. Ninety pairs have fewer than two residual signed
orbits and are unavailable to this worker. All 558 admitted pairs pass its
source variable and domain-clause caps. These are construction estimates;
none is a measured solver wall or memory peak.

For each admitted pair, the audit independently recomputes the number of
three-base/two-residual point multisets with repetition. Under a fresh uniform
nonidentity query model, it derives the necessary query count before the
first-moment upper bound permits a 50% chance of obtaining `K` rows, where
`K` is that pair's base-orbit count. It then retains every pair whose triple
`(conditional query floor, source domain clauses, source total variables)` is
not strictly dominated componentwise. Ties remain distinct configurations.
The result is **311 pairs in 24 unique metric shapes**. This is a
source-derived screening frontier; it does not establish a runtime Pareto
frontier or a factor-base winner.

| Selected pair or source alternative | Base columns | Residual orbits | Conditional query floor for its `K` rows | Domain clauses | Variables |
| --- | ---: | ---: | ---: | ---: | ---: |
| Frozen sequential-K600 / Gray-K600 arm | 600 | 97 | 33,975 | 183,314 | 95,068 |
| One strict source dominator: hash-K256 / sequential-K600 | 256 | 600 | 4,878 | 181,078 | 95,016 |
| Frozen sequential-K600 / hash-K600 arm | 600 | 600 | 889 | 269,830 | 97,080 |

The Gray arm has **114 strict source dominators** among the retained pairs;
the example above is lowest in query floor among them. The hash arm has none
and is tied with other K600/disjoint-pool configurations on these three
metrics. The frozen paired gate remains a controlled comparison because it
holds its K600 base fixed while changing the pool. The source dominance of
the Gray arm is not a runtime verdict: changing the base or pool also changes
SAT structure, relation dependencies, import work and linear algebra.

The audit pins the support screen, worker capacity source and manifest by
SHA256, rejects changed inputs, checks every pair's stored identity, and
recomputes exact multiset counts. It launches no N83 solver, adds no S3
objects, and leaves `selected_best_total_runtime` null. Reproduce with:

```sh
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/pilot-pair-pareto-20261014/pareto.py
python3 -m unittest discover -s research/koblitz_n83_factor_base_sweep_20261008/verification/pilot-pair-pareto-20261014 -p 'test_*.py' -v
```
