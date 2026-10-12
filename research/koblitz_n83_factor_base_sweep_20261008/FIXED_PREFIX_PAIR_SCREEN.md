# Fixed-prefix pair splitting: exact N83 coverage screen

This screen evaluates a specific five-summand split: choose three factor-base
points before observing a query, subtract their sum from the query, and ask a
complete two-point oracle whether the residual is a pair sum. The oracle is
assigned zero cost and never misses a pair. The model thus gives this fixed
schedule more help than an implementation can expect. Prefixes may repeat,
and each query contributes at most one full-smooth matrix row.

For a base of `B` distinct subgroup points, at most `P=C(B+1,2)` unordered
two-point multisets exist. A query uniform over the whole order-`r` subgroup
has a uniform residual after any predetermined prefix. One prefix therefore
hits with probability at most `min(1,P/r)`. A union bound gives at most
`min(1,A P/r)` for `A` fixed prefix checks on one query. Across uniform
queries, even if their outcomes are dependent, `E[successful queries] <= A
P/r` for `A` total prefix checks whose schedules are fixed before their
queries. Since matrix rank cannot exceed the number of successful queries,
Markov's inequality gives `Pr(rank >= K) <= min(1,A P/(r K))`. A 50% hit is
therefore compatible with the bound only after `ceil(r/(2P))` checks; 50% rank
`K` requires at least `max(K,ceil(Kr/(2P)))` total checks. These are necessary
counts, not predicted search work or achieved probabilities.

The [machine receipt](pilot-01/fixed-prefix-pair-screen.json) covers all 21
primary signed-Frobenius column sizes addressable in the v1 and v2 designs.
`B=166K` is conditional on complete construction except at K=2,048, where
the stored 339,968-point manifest and successful S3 download-hash receipt are
bound by SHA-256. Selected exact thresholds are:

| Signed columns K | Point count B | Pair multisets P | Checks before 50% single-hit bound | Checks before 50% K-rank bound |
| ---: | ---: | ---: | ---: | ---: |
| 600 | 99,600 | 4,960,129,800 | 243,728,666,055,352 | 146,237,199,633,210,983 |
| 1,182 | 196,212 | 19,249,672,578 | 62,802,409,480,827 | 74,232,448,006,337,231 |
| 2,048 retained object | 339,968 | 57,789,290,496 | 20,919,547,709,261 | 42,843,233,708,565,920 |
| 8,192 | 1,359,872 | 924,626,608,128 | 1,307,474,616,227 | 10,710,832,056,131,306 |
| 16,627 | 2,760,082 | 3,809,027,703,403 | 317,384,359,934 | 5,277,149,752,622,457 |

This rejects fixed-prefix enumeration as a short-pilot route to a useful
K=2,048 five-summand rank: even an exact, free pair oracle needs the stated
count before the upper bound permits 50% completion. It does **not** reject
splitting with target-adaptive prefixes, algebra shared across prefixes,
multiple independent rows per query, or partial large-prime cycles. It does
not compare wall time with rho or identify a factor-base winner. The existing
full-five-sum support bound is a different screen: it permits all five points
to vary simultaneously, while this one fixes the first three for each test.

Recompute and compare the committed receipt without overwriting it:

```sh
python3 research/koblitz_n83_factor_base_sweep_20261008/fixed_prefix_pair_screen.py
python3 -m unittest discover -s research/koblitz_n83_factor_base_sweep_20261008 -p test_fixed_prefix_pair_screen.py -v
```

The source reads only the frozen v1/v2 designs, K=2,048 manifest, and S3
round-trip receipt. Its small-group test enumerates the actual pair-sum image
to check the union-bound premise. The calculation spends no solver budget and
provides no new relation, rank or runtime measurement.
