# Affine substitution and exact enumeration: primary gains do not pass the universal gate

All four new treatments fail the preregistered requirement that every grouped
comparison have a paired 95% lower bound above 2.0 against the fixed reference
roster. Direct SIMD Gray evaluation passes **16/18** comparisons; the SIMD hybrid
with 16-variable leaves passes **15/18**. Neither result qualifies for confirmation.
The guarded launcher rejects the primary result before creating `run_02` or an
execution protocol. [CONFIRMATION_GUARD.json](CONFIRMATION_GUARD.json) records that
observed refusal. The previously recorded confirmation plan remains unexecuted.

This continues the implemented support-envelope and coefficient-transport work on
standalone generated Boolean systems. It implements general affine substitution,
exact original-variable recovery, bounded degree-preserving Boolean products, and
scalar/SIMD syndrome evaluation. The [contract](README.md) and
[correctness argument](CORRECTNESS.md) explain the algebra and its limits.

## Frozen acceptance results

The reference roster contains thirteen retained methods, an untraced packed
accounting control, the affine finalist pair and all three scalar enumeration
controls. Each candidate is compared with the pointwise fastest eligible reference;
the affine candidate removes itself from that roster. The three new SIMD treatments
are not reference arms in this preregistered experiment. Same-policy backend gates
use a separate 1.05 lower-bound threshold and cannot substitute for the full gate.

| Treatment | Strongest-reference groups passed | Universal >2x gate | Same-policy >1.05x gate |
|---|---:|---|---|
| `affine_sl_basis_fast` | 0/18 | REJECTED | PASS (18/18) |
| `gray_simd` | 16/18 | REJECTED | PASS (18/18) |
| `packed_gray12_simd` | 0/18 | REJECTED | PASS (18/18) |
| `packed_gray16_simd` | 15/18 | REJECTED | PASS (18/18) |

The affine compact implementation is 7.67–12.85x faster than its own list control,
but its grouped median ratio against the fastest reference is only 0.032–0.643.
General substitution's saved branching does not repay its total cost here.
The 12-variable hybrid also passes none of the eighteen strongest-reference gates.

Direct SIMD has grouped median ratios of 2.364–6.599, but the two failed intervals
expose variation hidden by those medians. The remaining failed groups for the two
promising treatments are retained explicitly:

| Treatment | Variables | Family | Split | Paired median ratio | 95% interval |
|---|---:|---|---|---:|---|
| `gray_simd` | 24 | cross_planted | regression | 2.364 | [1.978, 3.319] |
| `gray_simd` | 24 | cross_planted | holdout | 2.873 | [0.236, 5.913] |
| `packed_gray16_simd` | 16 | planted | regression | 1.968 | [1.796, 2.287] |
| `packed_gray16_simd` | 20 | planted | regression | 2.652 | [0.966, 4.795] |
| `packed_gray16_simd` | 20 | cross_planted | regression | 2.610 | [0.904, 4.724] |

In particular, the n24/cross-planted/holdout seed 2097153 has a direct-SIMD
per-fixture median ratio of **0.228740**. The other fresh seed, 20261016, has ratio
**6.839063**. A pooled median does not establish uniform improvement. These two
fixtures remain in the record; no rerun or seed exclusion supersedes the failures.
The intervals describe these fixed fixtures, not a population guarantee.

## Cold complete standalone solve costs

Units are milliseconds, including setup, all solver work, recovery, destruction
and independent result validation. This table pools the two fresh n24 holdout
seeds and 22 balanced arm-position repetitions. Its descriptive ratio uses the
untraced packed baseline on the planted family; it is a ratio of pooled medians,
not the paired strongest-reference acceptance statistic above. All rows have
class **engineering** and correctness PASS. These are standalone solver diagnostics;
calibrated full index-calculus costs are not available.

| Method | Planted (ms) | Cross-planted (ms) | Unplanted (ms) | Untraced packed / method, planted | Correctness |
|---|---:|---:|---:|---:|---|
| Search-only | 79.291375 | 80.110521 | 143.931354 | 0.156 | PASS |
| Degree-3 flat kernel | 261.296500 | 156.635917 | 510.276729 | 0.047 | PASS |
| Degree-3 sparse bucket | 932.561834 | 539.046000 | 1861.392605 | 0.013 | PASS |
| Degree-3 hybrid kernel | 323.089250 | 197.977687 | 826.881042 | 0.038 | PASS |
| Selective degree-2 flat | 90.098479 | 90.386396 | 165.173479 | 0.138 | PASS |
| Selective degree-2 one-word | 78.273562 | 77.630188 | 145.692750 | 0.158 | PASS |
| Ordered-specialization search | 69.558688 | 70.610833 | 126.581771 | 0.178 | PASS |
| Fixed-quadratic state | 28.146959 | 28.916208 | 51.876855 | 0.441 | PASS |
| Retained packed state | 17.568792 | 17.828833 | 32.550312 | 0.706 | PASS |
| Full RREF, list | 237.259479 | 48.409396 | 357.864375 | 0.052 | PASS |
| Full RREF, wide | 53.390625 | 11.743813 | 82.236834 | 0.232 | PASS |
| Echelon/affine tail, list | 125.170687 | 85.458624 | 229.339000 | 0.099 | PASS |
| Echelon/affine tail, wide | 40.682292 | 26.818541 | 70.524333 | 0.305 | PASS |
| Packed, diagnostic hashing removed | 12.401042 | 12.482542 | 22.925854 | 1.000 | PASS |
| Affine substitution and products, list | 266.847854 | 60.621437 | 471.918396 | 0.046 | PASS |
| Affine substitution and products, compact | 30.022771 | 7.068459 | 42.978333 | 0.413 | PASS |
| Direct Gray, scalar | 5.292791 | 13.534583 | 20.737104 | 2.343 | PASS |
| Direct Gray, SIMD | 0.802625 | 2.098875 | 3.210771 | 15.451 | PASS |
| 12-variable leaves, scalar | 12.976146 | 12.147833 | 28.111251 | 0.956 | PASS |
| 12-variable leaves, SIMD | 5.622458 | 5.275979 | 12.244438 | 2.206 | PASS |
| 16-variable leaves, scalar | 9.533771 | 9.073272 | 21.378792 | 1.301 | PASS |
| 16-variable leaves, SIMD | 1.691271 | 1.596583 | 3.761563 | 7.332 | PASS |

[SUMMARY.json](SUMMARY.json) records these costs for all four sizes, every paired
group interval, and individual-fixture ratios. [run_01/RESULT.md](run_01/RESULT.md)
retains the complete fixed-arm comparison. Earlier figures and verdicts remain in
the canonical scoreboard and their original frozen bundles.

## Correctness, accounting and evidence

The primary run completed **58,080 observations**: 120 distinct generated systems,
22 arms and 22 repetitions. The inputs include all 96 earlier systems and 24 fresh
holdouts. Each arm completes 89 SAT and 31 UNSAT fixtures. SAT models are checked on
the original equations; UNSAT is cross-checked with the completed retained search,
not an external proof certificate. Every scalar/SIMD pair preserves its model,
logical work and checksum. The untraced accounting control preserves its logical
search and reports a null trace. Enumeration points and leaves are explicitly
charged; zero tree nodes is not zero work.

The frozen producer receipt records **34 passing Rust tests**, including exhaustive
small affine-map, recovery, Boolean-product and Gray-evaluation checks, full-width
SIMD checks, word-boundary tests, and capped UNKNOWN cases. Thirty Python evidence
checks cover replay, all discovery bundles, changed-source/model/counter rejection,
phase accounting, censoring, confirmation refusal and duplicate-run rejection.
These are producer checks, not independent external review.

Seven source-frozen discovery probes and one separate native SAT feasibility probe
are retained. Native SAT's Python-wrapper timing is not a same-binary comparison.
The complete primary campaign took 2,477.655 seconds. Whole-worker peak RSS was
29,966,336 bytes and includes all arms and reference preparation. Candidate-specific
peak allocation remains unmeasured. See [ARCHIVE_NOTES.md](ARCHIVE_NOTES.md) for an
incidental interpreter cache retained because the original manifest sealed it.

[RUN_LEDGER.json](RUN_LEDGER.json) binds all experiment manifests, source lineage,
summary and the confirmation-guard receipt. Production solver cost, full
index-calculus cost, calibrated-operation ratio and rho ratio remain **null**.
No asymptotic, novelty or cryptanalytic breakthrough is established. The dramatic
gain objective remains open.

## Next bounded hypothesis

A size-based dispatch could select direct enumeration through 20 variables and
16-variable leaves at larger sizes. This is a new policy suggested by primary data,
not an unchanged-source confirmation and not a retroactive pass. Its design is
recorded in [NEXT_EXPERIMENT.md](NEXT_EXPERIMENT.md). It remains unimplemented and
unmeasured. No new run is launched by this report.
