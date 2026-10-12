# Row-space transport is correct, but does not beat the complete solver

Both coefficient-row-space policies fail their preregistered strongest-reference improvement and dramatic-gain gates. The retained packed solver remains the best pooled-median method on the fresh n24 families. Full-RREF wide/list backend ratios are **4.27–5.20x**, and echelon/affine-tail wide/list ratios are **3.03–4.59x**, but these matched-backend gains do not transfer into a robust gain against the strongest complete solver.

The full-RREF candidate passes an incremental strongest-reference comparison in 3/18 groups, but no dramatic comparison. Its paired median ratios span 0.199–2.667 with substantial fixture heterogeneity. The source-branching affine-tail candidate passes 0/18 incremental or dramatic comparisons and has paired ratios 0.235–0.952. Every matched-backend >1.05 gate passes. These are negative complete-method results, not an absence of correctness evidence.

## Implemented and tested

The [contract](README.md) establishes that specialization commutes with the GF(2) row span. Two independently implemented reducers carry quadratic polynomial bases through assignments: sorted-list symmetric differences and packed global coefficient rows. The full-RREF policy branches on the reduced basis; the second policy keeps quadratic rows in echelon form, reduces the affine tail, and separately carries original residual equations for the previous frequency heuristic. All of that extra state and work is charged. No polynomial multiplication closure is performed.

Each list/wide pair has identical models, logical counters, width histograms and trace digests. Different inference policies may choose different trees and models. Seventeen Rust tests include 442,368 small row-span/assignment cases, comparison with rebuilding from original equations, independent truth tables, dense word-boundary cases through n36, full solves and censored limits.

The two discovery resource probes are retained. Full RREF increased nodes on two n24/seed17 families. Preserving original-equation branching reduced nodes in the second probe, but still cost more than packed_state. All thirteen arms were then frozen into the full comparison with prior failed inputs and unused holdouts.

## Cold complete costs on fresh n24 holdouts

Units are milliseconds per complete solve plus result validation. The table pools two holdout seeds and thirteen balanced arm-position repetitions. Ratios below compare pooled medians with packed_state; acceptance uses paired pointwise minima over all twelve competing arms. Paired ratio summaries and ratios of pooled medians are different statistics. Every row has correctness PASS and class engineering.

| Method | Planted (ms) | Cross-planted (ms) | Unplanted (ms) | Packed / method, planted |
|---|---:|---:|---:|---:|
| Search-only | 66.217584 | 40.276479 | 77.027376 | 0.222 |
| Degree-3 flat kernel | 306.182646 | 181.077229 | 315.948980 | 0.048 |
| Degree-3 sparse bucket | 1209.731875 | 577.337396 | 1029.723521 | 0.012 |
| Degree-3 hybrid kernel | 564.794188 | 258.394334 | 402.785188 | 0.026 |
| Selective degree-2 flat | 77.067354 | 48.026376 | 89.620125 | 0.190 |
| Selective degree-2 one-word | 68.757750 | 42.481395 | 79.346271 | 0.213 |
| Ordered-specialization search | 58.147124 | 35.623875 | 67.757313 | 0.252 |
| Fixed-quadratic state | 24.245251 | 14.992396 | 27.964833 | 0.605 |
| Retained packed state | 14.673687 | 9.026291 | 16.798105 | 1.000 |
| Full row-space RREF, list | 125.438604 | 51.116333 | 268.847458 | 0.117 |
| Full row-space RREF, wide | 32.977396 | 12.047854 | 64.268354 | 0.445 |
| Echelon/affine-tail, list | 115.883104 | 74.345251 | 135.799709 | 0.127 |
| Echelon/affine-tail, wide | 36.544917 | 23.088125 | 44.602750 | 0.402 |

## Acceptance results

| Candidate | >2x strongest-reference gate | >1x strongest-reference gate | >1.05x same-policy backend gate |
|---|---|---|---|
| basis_wide | REJECTED | REJECTED | PASS |
| tail_wide | REJECTED | REJECTED | PASS |

[run_01/RESULT.md](run_01/RESULT.md) retains every paired interval and the per-size table. The regression split pools four previously used seeds, whereas the fresh holdout split has two. Legacy gates recomputed on this expanded grouping do not supersede earlier frozen failures. In particular, an older method's regrouped PASS is not a claim that each previously failed fixture now crosses its old threshold. No seed, policy or unfavorable cell has been removed from this run.

## A measured unused mathematical degree of freedom

The separate [affine diagnostic](affine_opportunity_01/summary.json) samples only twelve discovery fixtures after the timing campaign. Its instrumented solver exactly matches the uninstrumented tail_wide models, outcomes, every logical counter and trace. It finds **47,962 of 52,864 branching decisions (90.73%)** at consistent mixed nonlinear states with a nonunit affine relation and no remaining unit force. The rank sum at those states is **189,283**, averaging **3.95** affine rows per such state.

These constraints repeat across the search tree; they are not independent collected relations, and the percentage is not a holdout-population estimate. No affine-substitution speedup is measured. The current policies leave these general affine relations unsubstituted while quadratic rows remain.

The next experiment is general affine parametrization with exact recovery maps. The [prospective contract](AFFINE_SUBSTITUTION_NEXT.md) records coefficient updates, Boolean diagonal terms, composition, validation and full-cost requirements. Linear-variable substitution is an established ElimLin step, as described in [the primary literature](https://eprint.iacr.org/2014/893.pdf); no novelty claim is made. The implementation and its cost remain unmeasured. A source-only replay-wrapper correction is retained in affine_opportunity_02, whose raw counts and summary are byte-identical to the original census; the earlier bundle remains unchanged.

## Evidence and limitations

All **16,224 observations** across **96 distinct generated systems** complete with verified results: 71 SAT and 25 UNSAT fixtures per arm. Seventy-two inputs are retained and twenty-four are fresh. SAT assignments are checked on original equations. UNSAT requires a completed retained search reference; it is not an external proof certificate. No observed cell is censored. Resource-censoring tests preserve null completion costs and block promotion.

Twenty-one evidence tests cover frozen replay, both probes, retained failure inputs and fresh holdouts, changed-source/model/trace/counter rejection, strongest-reference selection, policy-group distinctions, censoring and the affine diagnostic. Whole-worker peak RSS is 29,360,128 bytes and includes all thirteen arms; candidate-specific peak memory is unmeasured. Logical operation counts and timing intervals retain their bounded meanings from the protocol.

[RUN_LEDGER.json](RUN_LEDGER.json) binds both discovery probes, the complete comparison and the later diagnostic. The canonical scoreboard preserves the earlier figures and carries this negative result. Production solver cost, calibrated-operation ratio, full index-calculus cost and rho ratio remain null. The dramatic-speedup goal remains open.
