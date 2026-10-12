# Whole-solve measurements favor cheaper search on this grid

The previous specialized-kernel gains do **not** transfer to these complete
generated Boolean solves. Always invoking the degree-3 hybrid kernel is slower
than the flat-kernel solver, despite identical search traces. Search-only often
finishes faster while visiting many more nodes. A specialization improvement
preserves that search exactly and yields a smaller **1.04–1.14x** paired median
gain against the fastest prior method. It does not meet the 2x gate.

All three grids finish every solve with verified results. No censored cell is
discarded. Every universal dramatic full-solve gate remains **REJECTED**. The
broader goal remains open, and no production/index-calculus claim follows.

## What was measured

The standalone driver uses exact specialization, affine propagation, fixed
branch selection and complete Boolean search. It is not the repository's
inherited-F4 implementation. It operates only on generated public systems,
with planted, cross-planted and unplanted families at n=12/16/20/24. Planted
witnesses are fixture construction data and are not supplied to the solver.

Cold totals include the original-system copy, propagation, every search node,
kernel work, solver/cache destruction and result validation. Fixture generation,
the independent search-status reference and diagnostic formatting are outside
arm timing and inside process receipts. SAT assignments are checked directly
against original equations. UNSAT verification requires the search-only
reference to finish UNSAT. Resource limits return UNKNOWN, not UNSAT.

The degree-3 flat/bucket/hybrid methods agree on outcomes, models, every logical
counter and trace digest. The selective degree-2 pair has the same equality
requirement within its policy. Search and ordered-specialization search also
match exactly. Trace digests are reproducibility diagnostics, not cryptographic
certificates; small systems are additionally checked by Boolean enumeration.

## Retained experiments

| Run | Mechanism and fresh comparison | Outcome |
|---|---|---|
| 01 | Search-only versus three always-degree-3 kernel methods | All 48 cells complete; search often wins despite more nodes; no universal 2x gate passes |
| 02 | Degree-2 inference only in at most ten active variables, with general-flat and one-word implementations | One-word method improves its matched policy by 1.13–1.39x but fails every 2x comparison against the strongest prior solver |
| 03 | Order-preserving specialization in search-only | Identical search/model/counters/digest; 1.04–1.14x paired median improvement, not 2x |

The first discovery-only resource probe is retained separately and is not
promotion evidence. Each run uses fresh holdout seeds and preserves every
earlier source bundle, raw sample, result and gate.

The one-word degree-2 kernel represents every possible monomial in at most
`1 + 10 + C(10,2) = 56` bits. It includes all admitted products of affine
generators, not just cancellation of the original rows. Its matched flat
control performs the same inference and search. This isolates representation
cost from a change in pruning policy.

Ordered specialization changes no decisions: zero assignments only filter
terms, and removing one common variable preserves order within that sublist.
A linear symmetric-difference merge replaces re-sorting. Multiple true
assignments use the original canonicalization fallback. The logical
specialization input-term counter is unchanged; it is not an instruction count.

## Latest complete costs at 24 Boolean variables

The table copies `run_03/RESULT.md`. Times are cold solve plus verification
milliseconds, with medians over two held-out seeds and seven balanced arm-
position repetitions. Ratios below compare pooled medians with the flat-kernel
solver. The acceptance gate instead compares the candidate with the pointwise
fastest of **all six** prior whole-solve methods. Every measured row has
correctness PASS and class engineering.

| Method | Planted (ms) | Cross-planted (ms) | Unplanted (ms) | Flat / method, planted |
|---|---:|---:|---:|---:|
| Search-only | 98.194229 | 86.760063 | 146.671084 | 3.033 |
| Degree-3 flat kernel | 297.822917 | 125.040250 | 597.041521 | 1.000 |
| Degree-3 sparse bucket | 1068.500479 | 466.614896 | 2137.548020 | 0.279 |
| Degree-3 hybrid kernel | 429.624542 | 176.332562 | 928.377542 | 0.693 |
| Selective degree-2 flat | 114.024709 | 99.104333 | 171.460291 | 2.612 |
| Selective degree-2 one-word | 99.382437 | 85.962625 | 151.653458 | 2.997 |
| **Ordered-specialization search** | **87.016291** | **76.643687** | **128.637125** | **3.423** |

Against the fastest prior method, the ordered-specialization paired median
ratios at n=24 are 1.129x [1.124, 1.132], 1.117x [1.110, 1.120] and 1.141x
[1.135, 1.148]. These intervals do not support a dramatic improvement. The
larger ratios against flat alone must not replace the stronger reference.

For example, latest planted holdouts use a median 24,182 search nodes versus
2,056 with degree-3 inference, yet search is much faster. The flat-kernel solver
spends about 96.8% of its solve time inside the kernel; the hybrid spends about
97.8%. Reduced node count is not reduced total cost. These findings concern
this square/overdetermined generated distribution, which differs from the
earlier eight-generator kernel panels.

## Verification and limits

The three complete grids contain **144 cells and 4,848 full-solve observations**.
Every arm finishes and is verified. Each grid has 36 SAT fixtures and 12 UNSAT
fixtures per arm; unplanted does not imply UNSAT. Repetitions are not distinct
random problems. Whole-worker peak RSS reaches 24,248,320 bytes and includes
all arms, fixtures and the independent reference.

Eight current Rust tests cover exhaustive small-system status/model agreement,
specialization semantics, censored limits, embedded-variable word-kernel
equality, matched policy traces and 6,912 ordered-specialization cases. Twelve
evidence tests replay all three frozen runs, check lineage and fresh holdouts,
reject corrupted models/traces/profiles, and confirm that synthetic resource
censoring produces null completion costs and blocks promotion. The separate
resource probe retains its source and executable hashes.

`RUN_LEDGER.json` and each run's manifest pin all evidence. No previous kernel
study or production solver is modified. Full IC cost, production inherited-
solver cost and rho ratios remain null.

The next concrete representation question is whether search can retain its
original quadratic coefficients and update only residual affine coefficients
and active-variable masks, avoiding per-node monomial-list reconstruction.
That mechanism is not implemented or measured here. It would need exact
specialization/deduplication, model and trace checks against the current fastest
search control before any performance claim.
