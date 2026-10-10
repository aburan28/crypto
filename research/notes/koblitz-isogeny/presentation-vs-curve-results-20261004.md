# Results: no curve effect detected in index-calculus cost (n = 16–23)

Design: [`presentation-vs-curve-design-20261003.md`](presentation-vs-curve-design-20261003.md).
Pilot: [`presentation-vs-curve-pilot-20261003.md`](presentation-vs-curve-pilot-20261003.md).
Data: `experiments/koblitz_presentation/`.
Analysis: `presentation_icc.py`, a two-way crossed random-effects model,
95 % CI from 2 000 bootstraps over curves.

## Amendments to the design (recorded before analysis)

1. **n = 31 dropped.** With a random V at l = 16 the pilot measured
   257 s/call and 802 reductions/call, which projects to 7.2 core-h per
   100 relations against the rule of ≤ 2. n = 31 is therefore out of
   budget for this solver.
2. **m ≥ 3 dropped.** m = 3 costs about 295 s/call even at n = 17, and
   m = 4 is not built by the solver.
3. **Samples.** n = 19: 128 curves. n = 23: 32 curves. The full classes
   did not fit the budget once random V turned out about 8× costlier per
   call than the monomial V.
4. **Cell order.** n = 16, 17 and 19 ran curves in order, so each curve's
   four subspaces ran back to back. n = 23 and an n = 19 re-timing ran in a
   seeded random order. Wall-clock metrics from runs in curve order are
   confounded with machine drift (see "Validity" below). Algebraic metrics
   (relations, reductions, first fall degree, |F|) are deterministic and
   unaffected.
5. **The n = 16 a₂ = 1 run duplicates a₂ = 0.** For even n the two
   families are isomorphic by y ↦ y + sx with s² + s = 1. The map fixes
   every abscissa, so relations, reductions and |F| are identical in all
   660 cells. The duplicate counts once and serves only as a timing-noise
   replicate.

Every cell recovered the transported d (BSGS-checked). There were 0
inconsistent relations across all 3 180 cells.

## Curve ICC, the share of cost variance attributable to the curve

Setting m = 2, l ≈ n/2, full-group probes, 4 random subspaces V per curve,
400 probes per cell.

| class | curves | yield / \|F\|² | reductions / call | log µs / call | log projected solve cost |
|---|---|---|---|---|---|
| n=16, l=8 (165 = full class) | 165 | 0.000 [0, 0.021] | 0.000 [0, 0.052] | 0.001 [0, 0.061] | 0.000 [0, 0.061] |
| n=17, l=8 (full class) | 273 | 0.009 [0, 0.059] | 0.014 [0, 0.062] | 0.030 [0, 0.101]† | 0.000 [0, 0.041] |
| n=19, l=9 | 128 | 0.000 [0, 0.045] | 0.000 [0, 0.039] | 0.920 [0.82, 0.95]† **invalid** | 0.175 [0.05, 0.29]† **invalid** |
| n=19, l=9, shuffled re-timing | 32 | 0.000 [0, 0.000] | 0.000 [0, 0.078] | 0.000 [0, 0.038] | 0.000 [0, 0.000] |
| n=23, l=11, shuffled | 32 | 0.086 [0, 0.260] | 0.031 [0, 0.185] | 0.000 [0, 0.075] | 0.030 [0, 0.186] |

† Wall-clock metric from a run in curve order.

Across random subspaces the share of variance from V itself is at most
0.005 in every class. Once V is random, which particular random V it is
does not matter. The pilot's large gap was monomial V against random V,
not one random V against another.

## Validity: the n = 19 timing signal is an artefact

In the curve-order n = 19 run, time per call showed a curve ICC of 0.92,
while reductions per call showed 0.000 on the same cells. The completion
log shows slow drift in machine speed of about ±1.5 %, plus a block of the
last 32 cells at +22 %. Each curve's four cells ran back to back, so the
drift landed on curves. Dropping the last block still left 0.27.

Re-timing 32 curves in a seeded random order gives **0.000 [0, 0.038]**.
The curve-order timing ICCs are therefore excluded as confounded. The n = 23
run was restarted in random order for the same reason; the 4 cells it had
completed in curve order were discarded and kept on file. The n = 16 timing
replicate (a₂ = 1, identical algebra, independent timings) bounds pure
timing noise at ICC 0.069 [0, 0.14].

## Verdict against the pre-registered criteria

- **H1 (curve structure): not supported.** No valid metric in any class has
  ICC ≥ 0.10 with a lower bound above 0. The only signals that met the
  threshold were the n = 19 curve-order timings, and they are shown above
  to be run-order artefacts.
- **H0 (folklore, upper bound < 0.05 on every metric): not established by
  the strict rule.** The contention-free metrics (yield, reductions) have
  upper bounds of 0.02–0.06 at n ≤ 19, so they sit at or just above the
  0.05 line. At n = 23, 32 curves give bounds of 0.19–0.26. The design
  calls this **inconclusive**, and it is reported as such.
- **Plainly stated:** no curve effect was detected. Every point estimate is
  ≤ 0.09, and at n ≤ 19 the curve contributes at most about 5–6 % of cost
  variance at 95 % confidence. What dominates instead is the choice of
  presentation: monomial V against random V changes per-call cost by
  2.7× at n = 19 and about 150× at n = 31.

## Scope

This covers m = 2, l ≈ n/2, the S₃ Gröbner solver in this repository,
400 probes per cell, full-group probes and random F₂-subspaces, at
n = 16, 17, 19 (|E| ≈ 2^16–2^19) and n = 23 (32-curve sample). It does
**not** cover n ≥ 31, m ≥ 3, other decomposition strategies or full solves.
The projected solve cost is (unknowns + 4)/yield × per-call time, not a
measured solve.

To tighten H0 at n = 23 and reach n = 31 would need a faster S_{m+1}
solver, or more compute than this session's 4 cores.
