# Exact pairing transfer priced against the frozen extension-field ECDLP

**Status:** bounded two-cell native solve and independent scalar replay complete; performance evidence is shared-host L0. The [protocol](PROTOCOL.md), [V2 amendment](AMENDMENT_V2.md), [V1 freeze](FREEZE.md), and [V2 freeze](FREEZE_V2.md) preceded their respective runs. This report uses the [V2 run](reports/run-20261015T0555Z/RUN_REPORT.md); the [V1 run and its coarser phase ledger](reports/run-20261015T0546Z/RUN_REPORT.md) remain intact. The V2 source and binary hashes, host, exact commands, exits, and raw outputs are in the frozen records, [`RUN_V2.sh`](RUN_V2.sh), [`runs_v2/`](runs_v2/), [`SUMMARY_V2.json`](SUMMARY_V2.json), and [`SHA256SUMS_V2`](SHA256SUMS_V2).

## Exact transfer and answer check

The public source was `E:y²=x³+1` over `F_(p²)` with the exact frozen generator and target coordinates in [`INPUTS.jsonl`](INPUTS.jsonl). The two prime-order subgroup cases were `(p,ν,r)=(61,2,19)` and `(271,3,61)`; their source ICV1/EC1 identities remain in the [fixture](../extension_field_20261014/FIXTURES.jsonl). The solver input omitted `answer_scalar`. It constructed `K=F_(p²)[z]/(z^9-2)` or `F_(p²)[z]/(z^5-2)` and checked each polynomial with the Frobenius/gcd irreducibility criterion before field inversion. The evaluation-field degrees over `F_p` are 18 and 10. At `p=61`, the root-of-unity output lies in a degree-9 subfield, so the degree-18 arithmetic is a priced implementation choice rather than an optimized destination-field solver.

Deterministic point search chose the first nondegenerate `r`-torsion pairing input at candidate x number 2 or 4. The curve-order calculations gave `r`-adic valuation 2 in both fields; the selected H, modulus, irreducibility receipt, α=`t_r(G,H)`, and β=`t_r(Q,H)` are pinned in [`SELECTION_61.json`](SELECTION_61.json) and [`SELECTION_271.json`](SELECTION_271.json). Each timed row rebuilt those values, checked `α^r=β^r=1` and `α≠1`, then used native baby-step/giant-step to solve `β=α^k` without reading the answer. The six returned scalars were `k=12` and `k=10` in their respective cells. A separate native process with independent `F_p`, `F_(p²)`, and curve formulas checked `[k]G=Q` for every row; the fixture scalar matched only in the post-solve audit. Checks on the public multiples 2, 3, and 5 passed 18/18. Rho and rho A/A returned and verified the same six target scalars in all 12 controls.

The Rabin and group computations are exact finite checks in the recorded implementation. This report introduces no new general pairing theorem. It does not claim a proof of a new algebraic identity or a Lean formalization of the implementation.

## One unit, one-target cold boundary

The fixed reference is same-target negation rho. The formal two-class generic walk-step floor is `sqrt(pi/4)=0.886226925` per `sqrt(r)` in counted group operations. The table's unit is `S_CPU/GAE = current-thread CPU ns / (paired rho's measured ns per source addition × sqrt(r))`. This is a local CPU-equivalent comparison: the formal counted-operation floor is **not** a hardware-time lower bound. Ratio-to-rho is the median of same-round complete CPU ratios; ratio-to-floor uses each row's median `S_CPU/GAE`. Each `[min,max]` is the range of three rounds, not a confidence interval. All rows have the same public `G,Q,r`; every row's scalar answer was verified.

| `(p,r)` | Variant and cost interval | `S_CPU/GAE` median [min,max] | Median / floor | Median paired CPU / rho | Correct |
|:--|:--|--:|--:|--:|:--|
| `(61,19)` | Negation rho, cold | 64.01 [39.01, 98.06] | 72.23 | 1.00 | 3/3 |
| `(61,19)` | Negation rho A/A, cold | 50.12 [34.61, 71.05] | 56.56 | 1.11 | 3/3 |
| `(61,19)` | Pairing, cold solve + scalar replay | 163,807.34 [161,020.04, 170,683.25] | 184,836.79 | 2,666.48 | 3/3 |
| `(61,19)` | Pairing + all post-solve checks | 170,737.19 [169,493.47, 178,065.56] | 192,656.29 | 2,781.81 | 3/3 |
| `(271,61)` | Negation rho, cold | 32.76 [28.97, 33.27] | 36.97 | 1.00 | 3/3 |
| `(271,61)` | Negation rho A/A, cold | 40.90 [28.42, 74.12] | 46.16 | 1.23 | 3/3 |
| `(271,61)` | Pairing, cold solve + scalar replay | 21,978.75 [21,607.04, 28,981.92] | 24,800.37 | 660.71 | 3/3 |
| `(271,61)` | Pairing + all post-solve checks | 23,193.81 [22,872.96, 30,304.74] | 26,171.41 | 698.15 | 3/3 |

![Paired full-validation CPU ratios](GRAPH_V2.svg)

The fully validated workflow includes an extra post-solve bilinearity check used for research validation. The cold DLP method row charges construction, every failed and successful H candidate, both pairings, the finite-field DLP, and the independent source-scalar replay. Its setup reuse count is one target. Both rows are shown so the optional research check cannot move cost invisibly. The finite-field solve is an exact native BSGS on order 19 or 61; it is not an asymptotically optimized finite-field DLP. The measured one-target cost is above rho in both cases, including without the extra check. The six-class automorphism-aware rho already measured on the same targets in [V4](../extension_field_20261014/REPORT.md) is contextual, since it was run in a different session and is not a paired control here; no cross-session speed ratio is taken.

## Phase and primitive ledger

The following are **median current-thread microseconds**, with each phase sorted independently across three rounds. The complete total therefore need not equal the sum of displayed medians. Full row-by-row counts, timer residuals, and process memory are in `SUMMARY_V2.json`.

| `(p,r)` | Setup | H search | Miller | Final exponent | Pairing order checks | Field DLP | Scalar replay | Extra checks | Full total |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| `(61,19)` | 30.167 | 101,723.959 | 969.000 | 759.625 | 35.000 | 48.917 | 2.167 | 4,497.833 | 108,490.417 |
| `(271,61)` | 38.125 | 23,873.166 | 209.833 | 143.333 | 13.001 | 19.334 | 2.875 | 1,362.042 | 25,681.333 |

Across every round of a cell, the exact primitive totals were identical: `16,491,269 F_p` multiplications and `46 F_p` inversions at `p=61`, or `5,602,853 F_p` multiplications and `64 F_p` inversions at `p=271`, counting the extra checks. The full sequence used `29,538` or `41,840` extension multiplications, including squares. The H-search phase used `26,469/27,169` or `39,115/39,673` extension multiplications within selection; Miller used 386 or 322, and final exponentiation used 300 or 216. The corresponding selection-phase `F_p` multiplication counts were H search `15,050,300` or `5,188,085`; Miller `205,647` or `52,480`; final exponentiation `169,770` or `37,500`. The paired rho median used 949 or 1,557 `F_p` multiplications. Each `F_p` inversion is also reported as its own primitive; its internal arithmetic is priced by CPU, not silently treated as zero. The candidate-search phase includes square-root tests, failed x candidates, cofactor projection, and exact torsion/order checks.

The H-search charge also includes the integer group-order recurrence and cofactor calculation. It used 148 or 293 extension inversions; Miller used two in each cell. These subphases remain in the cold total even though the group-order calculation is labeled inside H search rather than the field-only setup row.

Pairing peak process RSS had medians 6,586,368 and 6,602,752 bytes; rho peak medians were 6,389,760 and 6,406,144 bytes. All runs stayed below the 512 MiB cap. The 10-million-extension-multiplication, 4,096-x, 512-acceptable-H, and 15-minute current-thread CPU caps did not trigger. The V2 native accountant checks every phase's multiplication and inversion subtotals against the aggregate selection count. V1 and V2 retained exactly the same deterministic answer, field-operation totals, and rho counted GAE in all six cells; their CPU timings differ as descriptive repeated measurements, not an algorithm change.

## Interpretation and unresolved tests

**Observation:** H search dominates these two cold runs, while the multiplicative DLP itself took median 48.917 or 19.334 µs. The result is an executable transfer to a finite-field DLP, with a complete verified scalar solve on both frozen public targets. **Interpretation:** these two cells do not show an end-to-end cost improvement over rho. The phase ledger motivates a separate reusable-H/α experiment, but making H free would change the cold single-target contract; any amortized result must state an actual target count and charge setup once. Even the measured post-search work needs a new run with precomputed H to establish a warm per-target cost. No exponent or deployed-curve inference follows from orders 19 and 61.

The still-missing nonfactorable χ polynomial decomposition oracle is a different torsion-fiber test: implement it in native Rust, exhaustively compare returned decompositions with independent pair enumeration at the smallest frozen instance, then complete relation collection, recovery, and matched controls on all four frozen targets. The [χ gate](../extension_field_20261014/NEXT_GATE.md) remains open. Isogeny attribution also remains open: the present degree-three destination has a public scaling isomorphism over every `F_(p²)` in this family, so a new explicit base-field non-isomorphic source–destination pair and paired solved relations are needed. Those are separate experiments, not results of this pairing run.

The current CPU measurements exclude input parsing, output serialization, and the source-add calibration from both candidate and rho. Shared-host L0 rho CPU varied substantially under A/A (`0.353–1.285` paired ratio at `p=61`, `0.981–2.262` at `p=271`). The table reports exact observed ratios but no isolated wall-performance claim. A scaling claim would need at least four increasing subgroup sizes, comparable minimal-field destination solvers, isolated paired timing, and an independent host. This finite two-cell result remains a measured correctness and cost boundary for the implemented route. No canonical index-calculus scoreboard or prior V4 graph was changed: they record their own frozen relation runs, and this pairing study is a separate method with no qualifying ledger promotion.

## Requirement status

| Requested item | Status and evidence |
|:--|:--|
| Exact H and nondegenerate pairing on `p=61,r=19`, then `p=271,r=61` | Verified in frozen selections and six V2 timed replays. |
| Multiplicative DLP without passing the answer, independent `[k]G=Q` | Verified six times, plus post-solve fixture comparison and 18 pairing checks. |
| Fully charged cold setup, search, Miller, final exponent, field solve, verification, memory, rho | Verified at shared-host L0 in `SUMMARY_V2.json` and both run reports; source-add calibration and file I/O are disclosed exclusions. |
| Nonfactorable χ polynomial oracle with exhaustive control and complete four-target relations | Pending native implementation and run under the separate frozen χ gate. |
| Non-isomorphic extension-field source–destination attribution pair | Pending construction, map verification, and paired solved relation rows. |
| Growing-size, isolated, independent-host speed claim | Not attempted; the exact acceptance gate is in the protocol. |
