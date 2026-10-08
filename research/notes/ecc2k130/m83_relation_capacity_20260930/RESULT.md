# m=83 natural m=3 controls had negligible support capacity

**Decision: `M83_M3_FIXED_BASES_UNDERPOWERED`.** The archived full-orbit
factor bases contained 332 and 498 physical subgroup points. Even if every
unordered three-point sum were distinct, their one-target support could cover
at most `2.54535220e-18` and `8.56483486e-18` of the m=83 prime-order
subgroup. For the four independently sampled natural targets per seed, the
union bounds are `1.01814088e-17` and `3.42593394e-17`; jointly across all
eight targets the bound is `4.44407482e-17`. The recorded zero natural hits
and rank zero were therefore the expected outcome under that frozen uniform
target law, even with a perfect m=3 solver. They cannot discriminate FES,
SAT, F4/F5, or factor-base quality on natural yield.

This is exact integer **counting**, not a measured relation probability. For
`B` distinct physical subgroup points and unordered `m`-tuples with
repetition, at most `C(B+m-1,m)` tuples exist. Their group sums occupy no
more than `min(r,C(B+m-1,m))` targets in a subgroup of order `r`. Collisions,
projection, killed points and solver failures can only reduce that number.
For each uniformly distributed target the support probability is bounded by
this count divided by `r`; the four-target figure uses the union bound and
does not require the target draws to be independent. A fixed or adversarially
chosen target is outside this probability statement.

| Archived seed | Physical B | m | Tuple ceiling | One-target support ceiling | Four-target union ceiling | Observed natural m=3 hits/rank |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| 260938 | 332 | 3 | 6,154,284 | 2.54535220e-18 | 1.01814088e-17 | 0 / 0 |
| 260939 | 498 | 3 | 20,708,500 | 8.56483486e-18 | 3.42593394e-17 | 0 / 0 |
| 260938 | 332 | 4 | 515,421,285 | 2.13173247e-16 | 8.52692988e-16 | untested |
| 260939 | 498 | 4 | 2,593,739,625 | 1.07274557e-15 | 4.29098226e-15 | untested |
| 260938 | 332 | 10 | 5,127,780,552,883,038,908 | 2.12080033e-6 | 8.48320132e-6 | untested |
| 260939 | 498 | 10 | 282,830,784,825,729,136,425 | 1.16976071e-4 | 4.67904284e-4 | untested |

The following **necessary** physical-base sizes make the tuple-counting
ceiling reach 1% of the group for one uniform target. This 1% line is a
planning threshold, not a rank or hit prediction. The last column assumes a
literal 16-byte entry for every unordered pair; it illustrates one explicit
pair-table layout and is not a memory lower bound for an implicit solver.

| Field degree | m | Smallest B for 1% ceiling | Unordered pair entries | 16-byte pair-table scenario |
| ---: | ---: | ---: | ---: | ---: |
| 83 | 3 | 52,544,464 | 1,380,460,374,795,880 | 20,570,462.6 GiB |
| 83 | 4 | 872,790 | 380,881,628,445 | 5,675.6 GiB |
| 83 | 7 | 5,325 | 14,180,475 | 0.211 GiB |
| 83 | 10 | 780 | 304,590 | 0.00454 GiB |
| 131 | 3 | 3,443,553,994,135 | 5,929,032,055,263,277,584,196,180 | 8.835e16 GiB |
| 131 | 4 | 3,574,951,633 | 6,390,139,590,932,159,161 | 9.522e10 GiB |
| 131 | 7 | 617,668 | 190,757,187,946 | 2,842.5 GiB |
| 131 | 10 | 21,837 | 238,438,203 | 3.553 GiB |

This removes tiny-base m=3 natural-yield comparisons at m=83 from the
promising path. A useful next m=83 confidence test must first pick a physical
base and arity with adequate **counting capacity**, then measure actual
distinct support, a natural positive/negative PDP corpus, rank progress,
full-log recovery and matched rho with every phase charged. At m=83,
`B=780,m=10` merely clears the chosen 1% *ceiling*; it does not guarantee any
relations or a practical Semaev solver. At n=131 the corresponding m=10
threshold is 21,837 physical points, and even a 16-byte pair materialization
would use 3.553 GiB before the rest of an attack. An implicit high-arity
representation may avoid that pair table, so the separately review-gated
[m10 capacity attempt #937](https://github.com/aburan28/crypto/pull/937)
remains an independent prerequisite rather than a failed consequence of this
audit.

The [protocol](PROTOCOL.md) was committed on [PR #1094](https://github.com/aburan28/crypto/pull/1094)
before the machine result. [BOUNDS.json](BOUNDS.json) records the exact tuple,
pair and probability-bound numerators and denominators plus frozen input and
derivation-source hashes. The independent [REPLAY.json](REPLAY.json) checks
all six archived rows, all eight minimal thresholds against their predecessor,
the original pair counts and input hashes, and exhaustive multiset support in
three small cyclic groups. Its derived JSON SHA-256 is
`da364105c77e8ef04c11d57430f076c23b660e75560e479c90ab9b7995226eac`.
Reproduce from the repository root with fresh paths:

```sh
python3 research/notes/ecc2k130/m83_relation_capacity_20260930/derive.py \
  --out /tmp/m83-capacity-bounds.json
python3 research/notes/ecc2k130/m83_relation_capacity_20260930/verify.py \
  --input /tmp/m83-capacity-bounds.json --out /tmp/m83-capacity-replay.json
```

The old m=83 solver-stage results and their failed/timeout cases remain
unchanged. No discrete logarithm, useful natural relation, m=83 full-rank
matrix, end-to-end cost `S`, matched rho ratio, or n=131 speedup is established
by this counting audit. The ECC2K-130 feasibility goal remains open.
