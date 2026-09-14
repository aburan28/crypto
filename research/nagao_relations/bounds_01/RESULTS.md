# Bounds identify coefficient branching as the next obstacle

The existing hybrids cannot beat the measured direct S3 table in field-multiplication count on complete enumeration of any of the ten frozen eight-target batches: even their branch-only lower bounds exceed the table's measured total. The smallest gap is 15.73 times and the largest is 53.03 times. These are component bounds, not timed hybrid completions, total-operation speedups, or first-hit comparisons.

[Proofs and assumptions](PROOFS.md) · [next experiment](NEXT.md) · [raw](raw.jsonl) · [summary](summary.json) · [audit](audit.json)

## Exact matched-work comparison

All columns below are field-multiplication calls for the same complete eight-target batch at the stated support. Hybrid entries are lower bounds excluding setup and recovery; S3 entries are measured totals including setup and verification. SAT work has no comparable completed multiplication measurement. A timeout has not been converted into a completed result.

18-bit field:

| Variant | prefix d8 | prefix d9 | prefix d10 | f4-stable d8 | f4-stable d10 |
|---|---|---|---|---|---|
| hybrid-reference | >= 3,641,372 | >= 14,622,790 | >= 58,605,596 | >= 3,641,372 | >= 58,605,582 |
| hybrid-filtered | >= 3,641,372 | >= 14,622,790 | >= 58,605,596 | >= 3,641,372 | >= 58,605,582 |
| pair-invariants-s3 | 146,240 | 501,290 | 2,365,916 | 68,666 | 1,951,583 |
| chained-s3 SAT | unmeasured | unmeasured | unmeasured | unmeasured | unmeasured |
| symmetric-S4 SAT | unmeasured | unmeasured | unmeasured | unmeasured | unmeasured |

30-bit field:

| Variant | prefix d8 | prefix d9 | prefix d10 | f4-stable d8 | f4-stable d10 |
|---|---|---|---|---|---|
| hybrid-reference | >= 3,641,372 | >= 14,622,748 | >= 58,605,624 | >= 3,641,414 | >= 58,605,610 |
| hybrid-filtered | >= 3,641,372 | >= 14,622,748 | >= 58,605,624 | >= 3,641,414 | >= 58,605,610 |
| pair-invariants-s3 | 207,523 | 501,965 | 1,391,380 | 231,566 | 1,395,529 |
| chained-s3 SAT | unmeasured | unmeasured | unmeasured | unmeasured | unmeasured |
| symmetric-S4 SAT | unmeasured | unmeasured | unmeasured | unmeasured | unmeasured |

| Bits | Base | d | Hybrid floor / measured S3 multiplication count | Measured S3 / its conservative floor |
|---|---|---|---|---|
| 18 | prefix | 8 | 24.9000 | 2.4345 |
| 18 | prefix | 9 | 29.1703 | 2.1181 |
| 18 | prefix | 10 | 24.7708 | 2.3349 |
| 18 | f4-stable | 8 | 53.0302 | 3.2275 |
| 18 | f4-stable | 10 | 30.0298 | 2.3030 |
| 30 | prefix | 8 | 17.5468 | 2.8588 |
| 30 | prefix | 9 | 29.1310 | 2.0112 |
| 30 | prefix | 10 | 42.1205 | 1.5341 |
| 30 | f4-stable | 8 | 15.7252 | 2.7180 |
| 30 | f4-stable | 10 | 41.9953 | 1.5326 |

Every row is an accounting/architecture diagnostic. The multiplication bounds leave additions, squarings, binary work, allocation and control outside this one axis. They do not provide a calibrated common-operation speedup or rho ratio.

## Yield ceilings

The F4-defined curve orders are exactly 262,926 at 18 bits and 1,073,781,414 at 30 bits. Both are composite; these probabilities are for uniform affine targets in the full curve group. No prime-subgroup probability is inferred.

| Bits | Base | d | K abscissas | Signed triple mass after exclusions | Success ceiling | Uniform attempts lower bound |
|---|---|---|---|---|---|---|
| 18 | prefix | 8 | 127 | 2666966 | 100.000000% | >= 1.000 |
| 18 | prefix | 9 | 256 | 22107836 | 100.000000% | >= 1.000 |
| 18 | prefix | 10 | 534 | 201888892 | 100.000000% | >= 1.000 |
| 18 | f4-stable | 8 | 74 | 518574 | 100.000000% | >= 1.000 |
| 18 | f4-stable | 10 | 488 | 153999096 | 100.000000% | >= 1.000 |
| 30 | prefix | 8 | 140 | 3580640 | 0.333461% | >= 299.885 |
| 30 | prefix | 9 | 263 | 23979288 | 2.233163% | >= 44.780 |
| 30 | prefix | 10 | 505 | 170698078 | 15.896911% | >= 6.291 |
| 30 | f4-stable | 8 | 152 | 4590400 | 0.427499% | >= 233.919 |
| 30 | f4-stable | 10 | 506 | 171714538 | 15.991573% | >= 6.253 |

These are ceilings, not observed success rates. They correct the earlier signed triple count for infinity and target-abscissa collisions. Signed multiplicity differs from projected x-tuples; retaining multiplicity is essential for the proof.

**Negative trace result:** every one of the 54 trace classes has positive admissible triple mass on every larger base. The trace-fiber ceiling equals the corrected global first-moment ceiling in all ten cases. Thus this quotient adds no stronger global success ceiling or whole-target rejection on this panel. It does not prove that finer coefficient constraints are useless.

## Necessary scaling thresholds

To make even the loose counting ceiling reach 50%, one needs 8*C(K,3) >= (N-1)/2. The dimensions below assume as optimistically as possible that every nonzero abscissa in V lifts. Real point support, exclusions and sum collisions can require more. Larger rows are exact derived thresholds, not solver experiments.

| Field bits | Necessary K for 50% ceiling | Necessary binary dimension | Current S3 entries at that K |
|---|---|---|---|
| 18 | 48 | 6 | 2,256 |
| 30 | 740 | 10 | 546,860 |
| 42 | 11,816 | 14 | 139,606,040 |
| 60 | 756,156 | 20 | 571,771,140,180 |
| 90 | 774,301,844 | 30 | 599,543,344,847,498,492 |
| 126 | 3,171,540,345,417 | 42 | 10,058,668,162,604,612,132,558,472 |

## Architecture bound and next step

Under independent uniform targets, complete base scanning per query, and a requirement for a constant multiple of M independent relation rows, the explicit pair-table architecture incurs Omega(M^2 + N/M), minimized at Omega(N^(2/3)). This is a conditional bound for that architecture. It is not a universal index-calculus lower bound, a measured exponent fit, or a full-DLP rho comparison.

The next proposal is **coefficient-block rejection before branch solving**. The general-coefficient normalization in NEXT.md simplifies the target constant to B/r^2 and exposes a restricted linear equation. The existing fixed-branch pullback alone still visits a quadratic number of branches. A block relaxation must prove an entire block inconsistent at less cost than visiting it; whether that happens is unproved. No such block solver was implemented in this round.

On the 30-bit prefix d10 batch, surviving branches with the present seven-call multiplication floor would require at least 97.6259% pruning to match the S3 multiplication total, or 98.1007% pruning to be 20% below it. These necessary thresholds optimistically assign zero cost to pruning, setup and recovery. A different per-branch formula requires a new bound. They are not a forecast that such pruning will succeed; exact rational thresholds are in summary.json.

## Verification and reproduction

Exact signed-triple enumeration and instrumented original branch counts agree on 265 tiny target/space cases, including three fresh spaces. The larger audit covers all 24 cold targets, 80 batch targets and 60 new uniform targets. The normalization passes 2,048 exact residual-identity checks plus 756 exhaustive nonzero-target constant checks across all F4 coefficient choices at GF64. All source/input hashes match and there are zero failures.

`python research/nagao_relations/bounds_01/analyze.py` rechecks saved arithmetic, provenance, input coverage and prior-output agreement and regenerates these tables. To repeat expensive oracle/branch checks, run run.py from a separate checkout of its frozen source commit in publication.json; evidence files refuse overwrite. The identity checker has its own frozen source mapping. No solver performance candidate, timing comparison or full-DLP run was added here.
