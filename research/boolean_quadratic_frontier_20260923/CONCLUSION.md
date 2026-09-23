# Packed residuals improve the current frontier by 1.56–1.67x

One coefficient word per residual equation, a stack affine basis and exact hash-assisted deduplication improve complete generic Boolean solving by **1.56–1.67x** against the pointwise fastest of all eight retained methods. Every one of the 18 regression/holdout comparisons has a paired lower timing bound above 1.0. The incremental-improvement gate is **PASS**; the stricter current-frontier 2x gate is **REJECTED** (0/18 comparisons pass).

The cumulative comparison with the historical seven methods gives paired medians of **2.63–4.02x**, but its universal 2x gate is also **REJECTED**: 17/18 comparisons pass. The prior two failed confirmation groups now pass that cumulative comparison; a new n16 cross-planted holdout group fails. Neither an older, slower reference nor a favorable subset replaces the declared current frontier.

All **5,832 timed observations** across **72 distinct generated systems** complete with verified results. Search, ordered search, retained quadratic state and packed state have identical outcomes/models, logical counters and canonical trace digests. This is a representation improvement with the same mathematical search; it does not change the search-tree exponent.

## Mechanism and attribution

The [contract](README.md) describes the local coefficient layout and its exact update. Original quadratic terms occupy ordered local bits, followed by common linear-variable coordinates and a constant. Touch masks clear assigned-variable coefficients; XOR moves surviving quadratic contributions into linear or constant coefficients. Degree drops are detected directly. Equations in different local coordinates are compared by their actual monomials, never by raw-word equality or hash equality.

An occupied-pivot mask lets the affine reducer use a fixed stack array, avoiding temporary row vectors and scans of empty pivots. Hash buckets shortlist potential duplicate equations while exact support comparison preserves the original representative order. At most 32 original quadratic terms and q+n<64 fit the word; unsupported inputs use the retained quadratic method. All measured fixtures fit without fallback.

The discovery-only [profile](profile_01/summary.json) attributed median fractions of about 30% to tracing, 21–26% to affine propagation and 29–30% to specialization. Timer overhead perturbs those figures, and they are not promotion evidence. The candidate bundles three changes, so the complete-solve comparison does not isolate each change's individual contribution.

## Complete cold costs at 24 variables

Times below are milliseconds per complete solve plus result validation over the two fresh holdout seeds and nine balanced repetitions. The displayed ratio uses pooled medians of the retained quadratic method. Gates use paired pointwise minima over all applicable controls. Every row has correctness PASS and class engineering.

| Method | Planted (ms) | Cross-planted (ms) | Unplanted (ms) | Quadratic / method, planted |
|---|---:|---:|---:|---:|
| Search-only | 53.764896 | 44.279021 | 127.432626 | 0.357 |
| Degree-3 flat kernel | 253.129771 | 147.513229 | 527.080687 | 0.076 |
| Degree-3 sparse bucket | 1050.400459 | 549.571083 | 1767.499937 | 0.018 |
| Degree-3 hybrid kernel | 447.711125 | 253.588729 | 884.726334 | 0.043 |
| Selective degree-2 flat | 61.827375 | 51.415292 | 149.500770 | 0.310 |
| Selective degree-2 one-word | 56.283229 | 46.981979 | 131.702937 | 0.341 |
| Ordered-specialization search | 47.039854 | 38.542645 | 111.184270 | 0.408 |
| Retained fixed-quadratic state | 19.183375 | 15.664146 | 46.031709 | 1.000 |
| Packed residual coefficients | 11.966916 | 9.869708 | 28.152125 | 1.603 |

## Every retained and fresh comparison

The regression seeds are the full prior confirmation grid, including its failures. Fresh holdout seeds were frozen before measurement. Intervals bootstrap repeated timings on two fixed fixtures per group and do not estimate uncertainty across all Boolean systems. Both references require a lower bound above 2.0 for their dramatic gate.

| Split | Variables | Family | Current frontier ratio [95% interval] | Historical seven ratio [95% interval] | Current / historical 2x gate |
|---|---:|---|---|---|---|
| regression | 16 | planted | 1.641 [1.586, 1.665] | 3.638 [3.158, 4.068] | REJECTED / PASS |
| regression | 16 | cross_planted | 1.564 [1.531, 1.635] | 3.486 [3.264, 3.582] | REJECTED / PASS |
| regression | 16 | unplanted | 1.573 [1.539, 1.596] | 3.963 [3.842, 4.016] | REJECTED / PASS |
| regression | 20 | planted | 1.644 [1.623, 1.659] | 4.009 [3.931, 4.048] | REJECTED / PASS |
| regression | 20 | cross_planted | 1.621 [1.601, 1.635] | 3.909 [3.866, 4.009] | REJECTED / PASS |
| regression | 20 | unplanted | 1.608 [1.584, 1.626] | 3.971 [3.932, 4.022] | REJECTED / PASS |
| regression | 24 | planted | 1.648 [1.636, 1.663] | 3.941 [3.929, 3.951] | REJECTED / PASS |
| regression | 24 | cross_planted | 1.626 [1.607, 1.679] | 3.948 [3.920, 3.974] | REJECTED / PASS |
| regression | 24 | unplanted | 1.670 [1.659, 1.678] | 3.990 [3.972, 4.005] | REJECTED / PASS |
| holdout | 16 | planted | 1.639 [1.603, 1.665] | 3.741 [3.680, 3.812] | REJECTED / PASS |
| holdout | 16 | cross_planted | 1.653 [1.621, 1.694] | 2.632 [1.792, 3.957] | REJECTED / REJECTED |
| holdout | 16 | unplanted | 1.620 [1.591, 1.642] | 4.021 [3.886, 4.088] | REJECTED / PASS |
| holdout | 20 | planted | 1.626 [1.579, 1.653] | 3.942 [3.869, 4.056] | REJECTED / PASS |
| holdout | 20 | cross_planted | 1.631 [1.618, 1.657] | 3.555 [3.511, 4.027] | REJECTED / PASS |
| holdout | 20 | unplanted | 1.609 [1.592, 1.629] | 3.950 [3.885, 4.056] | REJECTED / PASS |
| holdout | 24 | planted | 1.630 [1.607, 1.689] | 3.938 [3.909, 3.954] | REJECTED / PASS |
| holdout | 24 | cross_planted | 1.597 [1.585, 1.621] | 3.911 [3.855, 3.938] | REJECTED / PASS |
| holdout | 24 | unplanted | 1.630 [1.616, 1.651] | 3.918 [3.895, 3.983] | REJECTED / PASS |

The remaining historical-reference failure is informative. On n16/cross-planted/seed20261014, packed search visits 323 nodes with a median cost of 183.083 microseconds. The old one-word inference method visits 59 nodes but costs 321.625 microseconds. The packed representation is faster, yet less than twice as fast as that stronger control on this fixture. The other seed has different costs; pooling them must not hide the failing lower bound.

## Validation and limits

Thirteen Rust tests pass, including every subset of three-variable affine rows in both insertion orders, all first and second partial assignments of 128 small quadratic polynomials, differing local coordinates, cancellations, exact duplicate handling, fallback at the word limit, complete model/trace equality and censored limits. Nineteen evidence tests replay the frozen run and discovery profile, retain prior failures and fresh seeds, reject changed sources/models/traces/counters, enforce the current frontier reference, and keep censored completion costs null. All five CLI configurations preserve their exact legacy arm order.

Each arm completes 52 SAT and 20 UNSAT fixtures. SAT models are checked directly on original equations. UNSAT requires a completed retained search reference, not an external proof checker. Small-system tests independently enumerate Boolean assignments. Trace digests are diagnostics. There are no censored benchmark cells.

Cold timing charges source-table compilation, state construction, propagation, tracing, recursion, solver destruction and result checking. Fixture creation, status-reference preparation and formatting remain outside arm timings and inside worker receipts. The 304-byte mutable state is smaller than the prior 384-byte state, but immutable tables grow. Whole-worker peak RSS reaches 28,950,528 bytes and includes all methods; candidate-specific total memory is unmeasured.

[RUN_LEDGER.json](RUN_LEDGER.json) binds the separate discovery profile and complete-solve run. The canonical scoreboard retains every earlier measurement. Calibrated operation ratios, production solver costs, full index-calculus costs and rho ratios remain null. The broader dramatic-speedup goal remains open.

The next mathematical direction is to preserve a quadratic row basis through assignments so that affine consequences emerge without rebuilding a full Macaulay matrix. Row combinations commute with specialization, but a useful implementation must price basis updates, density growth, tracing and validation. It would need a separate exact reference for its changed inference policy and a complete-solve comparison against packed_state. No benefit from that unimplemented direction is included here.
