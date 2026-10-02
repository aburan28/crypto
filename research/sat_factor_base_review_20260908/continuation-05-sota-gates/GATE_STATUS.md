# Koblitz index-calculus SOTA gate status

Current through Stage 189, 2026-10-02.

Stages 162–174 are the native-Boolean-F4 branch that culminated in the
size-gated contiguous-M4RI result on the already-opened public
`n=59, ell=9, m=3` true-negative. Stage 175 reconciles that branch with current
`main` and deliberately restores the repository's newer five-column
`BlockTables` engine rather than merging the older F4 fork. On the identical
target, the current engine is correct and exhaustive but regresses the Stage
174 medians by 1.3157x wall, 1.8071x CPU and 1.4622x RSS. Same-binary direct
MITM is much faster still. The result is retained as a negative current-engine
baseline, not rewritten as progress.

Stages 176–180 reject fixed-X1 schedule tuning, hash-indexed and dense-indexed
symbolic reducers, and the built-in four- and six-column table alternatives.
Every arm preserves the frozen equation fingerprint and exact exhaustive UNSAT.
The dense reducer index cuts its lookup count by 97.5 percent but misses the
frozen CPU threshold; four and six columns perform more actual XORs than the
five-column default.

Stages 181–182 add a separately controlled full-matrix block-8 M4RI path to the
current engine and certify its rank, pivot columns and row space. It reduces
actual elimination XORs from 147,794,583,858 to 99,192,937,526 and improves
multi-worker CPU, but misses the frozen paired wall gate and is slower on one
worker. It therefore remains available only through `F4_F2_FULL_M4RI=1`.

Stages 183–185 replace the current quadratic new-pair selection with exact
dense LCM groups and submask cover lookup on Boolean domains through 20
variables. The full-M4RI panel has median paired dense/quadratic ratios of
0.7368 wall and 0.8023 CPU. The decisive replay on current five-column
`BlockTables` has median paired ratios of 0.8078 wall, 0.8313 CPU and 0.9481
RSS. Stage 185 then builds selection commit `8014149a2` from a detached
checkout, passes ten F4 and three backend tests in both selected and quadratic
control modes, and verifies that an unset selector routes 1,011,275 updates
through the dense path. Dense selection is now the repository default;
`F4_F2_DENSE_PAIR_SELECT=0` retains the exact quadratic control. The clean
replay also caught the missing ignored lockfile. An initial root-tracking repair
was rejected because historical workflows deliberately copy and remove lock
files at that path before clean-tree assertions. The exact supplied lock is
instead force-added inside the Stage-185 evidence packet with SHA-256
`4f17b356fa7bac392b6d801d1c74fb9e36b6517f9465c8ebc19bb9a2792a84c5`;
the measurements and selected implementation are unchanged.

Stage 186 rechecks elimination after the pair-selector selection. Full M4RI
still reduces XORs and RSS but takes 1.8478x wall and 1.0182x CPU in the frozen
screen, so no confirmation is run. The selected Phase-B configuration is the
current repository F4 with five-column `BlockTables` and dense exact pair
selection; full M4RI remains a research control.

Stage 187 composes the current audit without double-counting inherited setup.
Stages 175–186 add 96 uniquely metered components, 9,028.458755 wall-seconds
and 22,105.703870 core-seconds. The campaign's measured lower bound through
Stage 186 is 508 components, 21,808.881290 wall-seconds, 52,387.262218
core-seconds and 6,310,576,128 bytes maximum RSS. Complete campaign cost stays
`null` because interactive compile/test work outside process meters remains.
The public reproduction issue contains 36 comments, all from `aburan28`, and
zero unaffiliated comments. Gates 1–3 remain partial, gates 4–5 retain finite
coverage, gate 6 remains false and gate 7 remains open. All-seven and SOTA
remain false. Stage 187 audit SHA-256 is
`1778162909b1722fe006b2721f6c42fb89793113158727ba4b9f97657fd2092d`.

Stage 188 re-tests the exact Stage 178 dense symbolic-reducer index after dense
critical-pair selection became the current F4 default. A Rust-native runner
replaces the legacy Python measurement path for the new experiment and owns
fresh-process CPU/wall/RSS accounting, watchdogs, hashing, terminal validation,
composition and replay. The index again cuts divisor operations by 97.5294
percent, from 4,190,633,182 linear tests to 103,532,494 exact submask probes.
The frozen three-pair confirmation nevertheless has median paired ratios
1.076402 wall, 0.972531 total core and 0.979483 RSS. It misses the unchanged
0.97 wall-and-core gate, so the candidate is rejected and reverted; the
selected runtime remains dense exact pair selection, linear reducer scan and
five-column `BlockTables`. Stage 188 adds a measured lower bound of 21
components, 898.725862 wall-seconds and 3,270.706466 core-seconds, taking the
cumulative lower bound to 529 components, 22,707.607152 wall-seconds and
55,657.968684 core-seconds at the unchanged 6,310,576,128-byte maximum RSS.
Complete cost remains `null` because the explicitly recorded native-meter
bootstrap work was not outer-metered. Final native verification passes 31/31
with result SHA-256
`8c3e1d1fc09f341434f71e28ff25b29119b4a0b9854d05f4568992d8a5bbdeb3`.
No SOTA gate changes.

Stage 189 keeps the twelve-way outer fixed-X1 batch but serializes product,
symbolic, packing, table-build and row-reduction sections inside each F4 call.
The exact current/outer-only screen has ratios 1.003837 wall, 0.953393 total
core and 0.946659 RSS. All algebraic and structural counters agree, but the
candidate misses the pre-registered below-1.00 wall gate; confirmation is
therefore prohibited, and the scheduling knob is reverted. Stage 189 adds a
measured lower bound of 11 components, 283.337375 wall-seconds and 1,258.115755
core-seconds. The cumulative lower bound becomes 540 components,
22,990.944527 wall-seconds and 56,916.084439 core-seconds at the unchanged
6,310,576,128-byte maximum RSS. Complete cost remains `null`. Final native
verification passes 19/19 with result SHA-256
`9b5e957686d9af59a8fad3bfaf0f0a58715876ab15019a9cd7476fbc31d05967`.
No SOTA gate changes.

The historical optimization chain is
`stage-99-optimization-chain-20260913/verification.json`. Stage 108 adds the
machine-replayable four-shard and direct-routing chain, five host-identified
routing comparisons, the selected five-pair `n=53` panel, and the refreshed
public unknown-scalar run. Stage 109 composed the current seven-gate audit on
2026-09-13 without changing any prior frozen artifact. Stage 124 re-sealed that
audit over the boundary ledger as updated on 2026-09-21 by the index-calculus
boundary ledger (operation-counted rows for the prime and binary regimes,
oracle pricing and the whole-process count on the Koblitz rows;
`research/notes/index-calculus/RESEARCH_IC_BOUNDARY_LEDGER.md`), chaining to
the Stage 109 seal. Stage 125 re-seals it again over that ledger's Round 2
of the same day (folded pair tables, walk targets, the exact counting
ceiling and a balanced Koblitz base; the note's §10, with the first-round
records kept in the ledger's history), chaining to the Stage 124 and Stage
109 seals. Stage 126 re-seals it once more over that ledger's Round 3, also of
the same day (a walk restart that costs one group operation instead of
thirty-two scalar multiplications, bases sized at the derived family optimum
on the prime and binary ladders, and the family shape law reported per row;
the note's §11), chaining to the Stage 125, Stage 124 and Stage 109 seals.
Stage 127 re-seals it a fourth time over that ledger's Round 4, also of the
same day: the unit's conversion from native counters to group additions now
uses ratios pinned in `docs/ic/calibration.json` rather than factors measured
on the host at the start of each run (the note's §12). That round reprices
every row by up to 9.5 per cent without moving a single native counter on any
of the 537 rows compared, and chains to the Stage 126, 125, 124 and 109 seals.
Stage 128 re-seals it a fifth time over that ledger's Round 5, on 2026-09-22:
the walk's sixteen restart offsets are drawn on first use instead of at setup,
so a row is charged for the offsets it took rather than for a pool it may never
reach (the note's §13). Because the offsets come off their own random stream and
are taken in a fixed cycle, the before and after walk bit-identical
trajectories, and all 537 rows are identical in every counter but the pool's
with zero trajectories moved; the saving is confined to the small rungs. It
chains to the Stage 127, 126, 125, 124 and 109 seals. No Koblitz gate fact moved
at any of these re-seals and no prior frozen artifact changed.

Stage 129 re-seals it a sixth time, also on 2026-09-22, and moves a gate fact
for the first time since Stage 109. The gate-6 online wall crossing holds only
against the fixture rho control, not against a batched signed-Frobenius rho
(see "The n=53 crossover against a batched rho" below), so gate 6 is now false
outright. The Stage 108 measurements themselves are unchanged. Stage 129 chains
to the Stage 128, 127, 126, 125, 124 and 109 seals, and no prior frozen
artifact changed.

Stage 130 supersedes the current audit without mutating those snapshots. It
hash-pins the frozen 160-input Phase-B matrix and the current scalar-blind
`n=53` workflow evidence, and keeps them explicitly separate: Phase B is the
common `n, ell, m` Semaev solver matrix, while the optimized `n=53` run uses a
subgroup-orbit factor base that has no equivalent WDSat/Magma `ell` encoding.
It also incorporates the current one-core refresh, native projected-predicate
operation counts, every accepted and rejected optimization process, and the
unchanged Magma and unaffiliated-review gaps.

Stages 130 onward were composed on a parallel branch as Stages 127 onward and
were renumbered by three on merge, because Stages 127 to 129 above were already
sealed. Their measurements are unchanged; only their stage ids, predecessor
links and seals moved, and Stage 130 chains to the Stage 129 seal.

Stage 131 adds a current unified `n=41, ell=6, m=3` single-target series over
the fixed public two-torsion-saturated Frobenius-union base.  It charges the
same process from algebraic materialisation through projected predicate,
relations, verification, sparse linear algebra, descent and automorphism rho,
under five default-thread and five one-worker repetitions.  It replaces the
older online-only reading with a current full-cost loss in both thread classes.

Stage 132 materialises the exact non-Frobenius-closed standard
`n=59, ell=9, m=3` Phase-B factor base in the end-to-end workflow.  A bounded
1.2-million-probe attempt reaches zero natural relations, as the public counting
ceiling predicts, while same-target signed-Frobenius rho completes in five
metered repetitions.  The result is censored at its cap: it is larger-regime
execution evidence, not a completed IC run or an UNSAT certificate.

Stage 133 counts the public cofactor classes of every standard `n=59` width
from `ell=9` through `16`, charges those counters and the width constructions,
and binds the result to the exact compact-table memory formula.  It establishes
a finite no-go frontier for materialized standard bases under the current 4 GiB
budget; other algebraic families, implicit bases and larger storage remain open.

Stage 134 completes one public `n=59` unknown-scalar workflow with the
cofactor-projected image of the standard `ell=14` base.  The recipe uses only
public cofactor multiplication and point identity deduplication; it retains no
scalar preimages or discrete-log labels and does not enumerate the target
subgroup.  Both the default-thread and one-worker runs recover and verify the
same scalar from the same 26,796 relations.  Fully charged IC remains much
slower than automorphism rho in both thread classes.

Stage 135 widens that same algebraic recipe to `ell=15` on the same target,
solver, collection window, sparse filter and rho seed.  It halves probes and
summand scans, reducing default-thread IC wall by 4.81 percent and one-worker
wall by 17.23 percent.  Peak RSS rises by about 4.8 times, and full IC still
loses to rho by three orders of magnitude on one worker.

Stage 136 holds the selected `n=59, ell=15` instance fixed while tuning its
collector.  Equal 102.4-million-summand arms reject windows 512 and 2048 in
favour of 1024.  Two source-pinned buffer-reuse pilots preserve every relation
but regress relation-unit wall, so that implementation is rejected and archived.

Stage 137 adds an explicit witnessed-compact pair table on that same instance.
Each four-byte compact rest carries a four-byte packed pair witness; every
candidate is re-added in the exact fast group, then exact candidates are sorted
before selection.  Both complete modes reproduce Stage 135's relation hash and
reduce full wall and core cost while raising peak memory by about 1.5 times.

Stage 138 caches each first-pass packed pair-sum key in row-major order and
derives its witness from the row and offset.  Scatter then avoids the second
quadratic group-addition pass.  Full time and CPU fall again in both thread
classes, while the temporary 4.34 GB key cache raises construction peak RSS to
about 10 GB.

Stage 139 doubles the witnessed presence filter from four to eight nominal bits
per pair in two exact pilots.  Both preserve the relation hash and slightly
reduce CPU, but relation-unit and build wall regress, so the wider filter is
source-pinned, charged and rejected; Stage 138 remains selected.

Stage 140 replays exact prefixes of the selected uniform relation stream.  Units
22--25 leave 5, 4, 1 and 1 projected columns absent.  After 25 units, the
workflow fixes a public factor-base point from missing column 15726 and asks the
same exact pair table for the remaining summands.  One relation from 10,000
direct lookups completes and certifies all logs.  Matched default and one-worker
full runs improve, with every discovery and rejected policy charged.

Stage 141 continues on the same public seed-59001 target.  Once raw column
coverage is complete, it ranks projected columns by their current public
relation incidence and forces exact relations through the 64 least represented
columns per round.  Fifteen uniform units plus 2.95 million direct pair lookups
derive all logs.  Against the selected Stage 140 result, default-thread full IC
wall falls 38.31 percent and whole-process CPU 29.77 percent; one-worker wall
falls 31.49 percent and CPU 31.51 percent.  The 12- and 14-unit successes and
both capped 10-unit failures are retained and charged.

Stage 142 corrects the mechanism attribution without changing the performance
result.  A rank-disabled replay emits the same 31,798-relation stream.  Every
one of the selected run's 295 forced-column attempts occurs while at least one
raw column is still uncovered, so the 15-unit winner uses Stage 140's public
uncovered-column mechanism.  The new 64-column rank fallback is exercised by
the 10-, 12-, and 14-unit frontier controls only.  The selected speedup is
therefore attributed to the 15-unit parameter choice, with the lower-unit rank
frontier retained as charged support.

Stage 143 tests stopping each raw-coverage attempt after the first deterministic
2,048-probe batch containing an exact relation.  It reduces targeted lookups
8.01 percent and preserves a complete verified unknown-scalar solve, but a
matched default-thread ABBA control regresses median full wall 4.24 percent and
whole-process CPU 0.19 percent.  The candidate is archived and rejected; the
Stage 142 selected implementation and result remain current.

Stage 144 delays sparse solving to every fourth rank-tail round on the 14-unit
frontier.  Repeated linear algebra falls from 14.81 to 4.93 seconds, but the
cadence overshoots to 20.03 million targeted lookups and full IC takes 122.78
seconds.  This is more than twice the selected absolute receipt, so the cadence
patch is archived and rejected without changing the selected result.

Stage 145 parallelises the 64 fixed-column attempts of each rank round.  A
nested-Rayon form and a single-layer form with serial per-column scans preserve
the same 29,948-relation stream and verified scalar, but full IC rises to
118.34 and 138.52 seconds.  Random access to the 10 GB pair table is
memory-bandwidth bound; both scheduling patches are rejected.

Stage 146 exposes homogeneous block-Wiedemann null-vector support as a possible
rank-tail target.  No failed prefix-14 solve returns such a witness, so every
round falls back to the existing global-incidence policy and reproduces the
Stage 145 relation hash.  The diagnostic patch is archived and rejected without
changing the selected implementation.

Stage 147 doubles witnessed-table presence-filter prefetch lookahead from 32 to
64 on a source-pinned four-run-per-arm panel.  All eight processes emit the same
2,157 relations, but candidate median relation-unit wall regresses 32.06 percent
and whole-process wall 12.04 percent.  Lookahead 32 remains selected.

Stage 148 halves the same lookahead to 16.  The four-run pilot is nearly neutral,
and a complete ABBA panel appears 2.36 percent faster overall.  The only changed
path, however, is 1.12 percent slower in relation-unit wall; the apparent win
comes from unchanged table-build and logs variation.  The candidate is rejected
on the causal path and lookahead 32 remains selected.

Stage 149 skips the duplicate presence-filter check after a window probe has
already been admitted.  A one-unit panel is strongly wall-positive, but the
minimal six-run full panel has 1.22 percent higher median relation-unit wall and
a paired-ratio median above one, despite small CPU reductions.  The helper and
minimal variants are both archived and rejected under the predeclared wall rule.

Stage 150 divides the existing presence-filter allocation into two independent
halves and requires both bits.  Exact relations and memory are unchanged, but
the second atomic set and filter load raise median unit wall 2.50 percent, unit
CPU 6.02 percent and build CPU 6.72 percent.  The split filter is rejected.

Stage 151 caches each block's presence-filter hash for both prefetch and
admission.  The extra 8 KiB scratch stream costs more than recomputing the
multiply-mix: median unit wall regresses 0.58 percent and CPU 0.31 percent.
The hash-cache patch is rejected.

Stage 152 separates admitted window probes from bucket consumption, prefetches
their bucket-index entries in FIFO order, and pins the unchanged exact relation
verifier out of line.  All paired default-thread relation-unit comparisons win;
median full IC wall falls 12.46 percent and whole CPU 4.36 percent.  A matched
one-worker pair improves full IC wall 10.97 percent and whole CPU 7.89 percent.
The exact relation stream and recovered scalar are unchanged, so this source is
selected.

Stage 153 adds a second prefetch for the compact-rest payload after the bucket
index arrives.  Its one-unit panel improves wall 30.25 percent, but complete
workflow relation-unit wall regresses 19.53 percent and exact verification
again slows roughly threefold.  The payload-prefetch patch is rejected and
Stage 152 remains selected.

Stage 154 halves the selected FIFO's internal block from 1,024 to 512 entries.
Median relation-unit wall falls 6.82 percent, but charged unit CPU rises 2.33
percent because smaller batched inversions spend more arithmetic.  The candidate
is rejected under the no-extra-work rule and Stage 152 remains selected.

Stage 155 doubles that internal block to 2,048 entries.  The larger working set
regresses median relation-unit wall 5.39 percent and CPU 0.95 percent, so it is
also rejected and the selected 1,024-entry block is retained.

Stage 156 retunes presence-filter lookahead on the selected FIFO from 32 to 64.
Median unit CPU falls 0.69 percent, but unit wall rises 2.87 percent and whole
wall 2.85 percent.  The wider lookahead is rejected and 32 remains selected.

Stage 157 retunes the same lookahead from 32 to 16.  Unit wall falls only 0.57
percent while CPU rises 3.25 percent, so the shorter lookahead is rejected and
32 remains selected.

Stage 158 skips the duplicate filter admission inside selected FIFO
consumption.  The minimal existing-body fast path preserves exact relations but
regresses median unit wall 1.99 percent and CPU 5.52 percent.  It is rejected.

Stage 159 halves bucket-index entries from `2^26` to `2^25`, doubling target
run length from about 16 to 32 while saving 128 MiB theoretically.  Unit wall is
neutral and CPU regresses 6.21 percent, so the current bucket width is retained.

Stage 160 doubles bucket-index entries to `2^27`, halving target run length to
about eight at a theoretical 256 MiB index cost.  Unit CPU falls 1.95 percent,
but wall rises 1.79 percent and median RSS about 9.3 percent.  It is rejected.

Stage 161 retests the eight-bit witnessed presence filter on the selected FIFO.
Median unit wall rises 6.06 percent, unit CPU 1.55 percent, build wall 12.71
percent and median RSS about 3.6 percent.  The four-bit filter remains selected.

The campaign has target-independent algebraic factor bases; matched native-XOR,
WDSat, CryptoMiniSat, direct-MITM, GGMP, and signed-Frobenius-rho controls; a
balanced 160-instance PDP panel through `n=59`; public unknown-scalar end-to-end
runs at degrees 23, 31, 41, 53, and 59; and same-target known-answer and
scalar-blind construction, rank, solve, and rho comparisons at `n=53`. It has not passed all seven gates and
does not establish a new state of the art.

| Gate | Status | Current evidence | Remaining requirement |
|:--|:--|:--|:--|
| 1. Charge every stage and resource | **Partial overall** | The corrected Phase-B matrix charges 21,038.596137 core-seconds. The current `n=53` campaign retains 318 measured processes / 1,395.640868 sequential wall-seconds / 3,192.181904 core-seconds. The Stage 141--161 frontier and optimization controls charge 168 processes: 8,805.864018 sequential wall-seconds, 46,940.856427 core-seconds and 10,620,682,240 B maximum RSS. These incremental charges do not replace Stage 140's 23 processes. One inherited preliminary harness process has 15.614915 s observed wall but no retained CPU/RSS receipt and remains explicitly excluded. | Execute licensed Magma with complete resources and independently recover or repeat the missing preliminary receipt; preinstalled OS/toolchain acquisition remains excluded. |
| 2. WDSat, CryptoMiniSat, Magma F4, MITM, GGMP | **Partial** | Native XOR SAT, WDSat, CryptoMiniSat, and direct MITM ran on the exact 160-input packet. Standard `n=31`/`n=41` and GGMP `n=31` are represented. The Stage-32 successor removes the two original WDSat buffer errors without rewriting Stage 26. | Execute all 160 Stage-22 Magma inputs on a licensed host under the frozen one-thread/no-retry contract, seal the return before truth scoring, and report F4 resources. |
| 3. Single-core, core-seconds, memory, conflicts, wall | **Partial because Magma is absent** | The selected bucket-prefetch panel has median default-thread IC 78.619846 s / 504.697434 core-seconds / 9,299,828,736 B median RSS and a matched one-worker candidate at 437.950836 s / 425.370001 core-seconds / 7,217,905,664 B RSS. Against their same-binary baselines, default IC wall improves 12.46% and CPU 4.36%; one-worker IC wall improves 10.97% and CPU 7.89%. | Supply the same fields for licensed Magma F4. SAT conflicts remain inapplicable to exact pair-table arms and are reported as null. |
| 4. Scale through `n=31`, `n=41`, and a larger PDP regime | **Satisfied for finite execution coverage including one completed larger IC run** | Phase B covers `n=31`, GGMP `n=31`, `n=41`, and `n=59`. Stages 132--133 retain the standard `n=59` cap and width frontier. Stages 134--141 complete and optimize the cofactor-projected `n=59, ell=15, m=3` public unknown-scalar workflow through persisted coverage and sparse-rank tails. | The evidence is finite and toy-sized; it is not an asymptotic scaling law or a literature-scale speed record. |
| 5. Unknown scalar with no constructed factor-base logs | **Satisfied for finite degrees 23, 31, 41, 53, and 59** | Stage 108 archives the `n=53` public hash-seed-53001 run. Stage 141 derives all 16,344 `n=59` logs from 31,727 uniform relations plus 71 public-matrix-selected relations, recovers `d=17861472351607`, and verifies `[d]G=Q` in default and one-worker modes. Neither target scalar nor factor-base logs are supplied. | Repeat on independent public seeds and obtain unaffiliated replay; these strengthen rather than replace the finite gate-5 execution. |
| 6. Full cost against automorphism-optimized Pollard rho | **False. The online wall gate passed only against the fixture rho control; full-cost/core gate false** | Stage 108 selects the four-shard direct route after five direct/mixed wins with 0.904540 median wall and 0.912199 median core ratios. Its identified EPYC 9V74 panel has five direct/rho wins and 0.839190 median wall ratio. Median direct core remains 2.261565 times rho, retained support is 738,197,504 B, and fresh build plus direct is 17.636692 times rho. The rho control in that panel (`koblitz_rho_fixture`) inverts once per step, canonicalises by squaring chains, stores every point and runs on one thread, about 9.5 µs a step. Against a batched signed-Frobenius rho on the same target (cryptanalysis `suite/examples/koblitz_batched_rho.rs`, M4 Pro, portable arithmetic for all arms), the selected direct takes a median 17.3 s against 0.49 s for 1-thread and 0.20 s for 4-thread rho. That is 46.3M support queries against about 0.56M rho steps. | Beat an automorphism-optimized rho (shared inversion, orbit key, distinguished points), not the fixture control, first on the EPYC gate host. The query count `r/(n·|F|)` grows as `r^{2/3}` against rho's `r^{1/2}`, so constant-factor improvements (including an orbit-folded support table) cannot close this gate toward larger degrees. |
| 7. Independent external reproduction and novelty review | **Missing** | Issue [#97](https://github.com/aburan28/crypto/issues/97) and [mtrimoska/EC-Index-Calculus-Benchmarks#1](https://github.com/mtrimoska/EC-Index-Calculus-Benchmarks/issues/1) now include the five-run selected panel, current-head source pin, unknown-scalar result, exact verifier boundary, and the `CONCUR` / `QUALIFIED` / `BREAKS` format. Stage 11 archives GitHub Actions run [34428320022](https://github.com/aburan28/crypto/actions/runs/34428320022) as a project-authored Linux degree-23 reproduction (`STAGE11_RESULTS.md`); it is not independent review. Reviewers bind each gate to per-instance records with [the evidence guide](EXTERNAL_REVIEW_EVIDENCE.md) and [review template](EXTERNAL_NOVELTY_REVIEW_TEMPLATE.json); empty fields remain requests, not completed review. | An unaffiliated reviewer must return a sealed reproduction and source-pinned novelty/correctness assessment. Project-authored CI and replays do not satisfy independence. |

The local measurements above are limited to their stated fields. The meter's
`single_core_seconds` field equals `total_core_seconds`, both user plus system
CPU; it is not a separate measurement of single-core elapsed time. Requested
thread counts and maximum per-process RSS do not by themselves establish
observed single-core execution or aggregate parallel peak memory. Inclusive
outer receipts must not be added to their nested process charges.

Field-size coverage also does not establish literature-instance equivalence.
The n41 ell5 m3 cell differs from WDSat's prominent n41 ell20 m2 experiment.
The n31 standard/GGMP comparison uses different curve coefficients and actual
base sizes (a=1, 31 points versus a=0, 63 points), so its timing ratio cannot
isolate the construction's effect. See [the parameter and sampling comparison](LITERATURE_MATCHING_20260910.md).

External reviewers can bind each gate assessment to individual backend and
instance observations using [the evidence guide](EXTERNAL_REVIEW_EVIDENCE.md)
and [review template](EXTERNAL_NOVELTY_REVIEW_TEMPLATE.json). The packet includes
WDSat, GGMP and 2025 SATIC prior-art questions; empty fields and invitations
remain requests for evidence, not completed external review.

## Current Phase-B matrix

The corrected 480-row SAT view has 141 true positives, 120 true negatives, 219
inconclusive outcomes, no false classifications, and no solver errors. Standard
`n=31` and `n=41` resolve all 40 targets under all three SAT backends. The
`n=31` GGMP cell yields 2 native, 1 WDSat, and 18 CryptoMiniSat true positives;
the rest are inconclusive. All three SAT backends are inconclusive on the 40
`n=59` targets under their caps. Direct MITM classifies all 160 inputs: 80 true
positives and 80 true negatives.

| Backend | Core-seconds | Summed process wall | Peak RSS | Conflicts or operations |
|:--|--:|--:|--:|--:|
| Native XOR SAT | 958.931674 | 1,278.114043 s | 140,333,056 B | 9,233,732 conflicts |
| WDSat, corrected view | 7,404.749336 | 9,524.938771 s | 24,252,416 B | 7,599,249 conflicts |
| CryptoMiniSat | 6,237.273520 | 8,137.581155 s | 122,703,872 B | 19,480,743 conflicts |
| Direct MITM | 759.710741 inclusive outer | 766.538655 s one-CPU elapsed | 85,479,424 B tree | 4,804,252 additions; 4,768,800 pair entries |

The original Stage-26 cells remain immutable. Stage 32 adds 426.718011 charged
core-seconds to turn the two `n=59` WDSat buffer assertions into clean
timeout-inconclusive terminals. The classifications do not change.

## Current n=41 and n=53 end-to-end results

Stage 39 fixes the `n=41` factor base algebraically as the two-torsion
saturation of the Frobenius union generated by masks `[1,2,4,8,16,32]`.
Construction uses zero target samples or discrete-log labels: 2,380 abscissae,
4,759 rational points, 60 signed orbits, and 29 projected columns.

It derives and certifies all 29 logs from 35 relations over 114,688 probes.
Five domain-separated public hash targets solve in 0.259641 seconds versus
0.909281 seconds for signed-Frobenius rho. Precomputation is 4.310245 seconds,
so the online ratio crosses while the amortized ratio remains 5.025825. The
four-core ABBA control reduces median precomputation to 2.010618 seconds and
complete five-target IC wall to 2.479844 seconds, while increasing median
process CPU by 28.7 percent.

Stage 68 retains the current point-defined `n=53` base with 9,964 points and 94
orbit columns, selected without scalar labels or target-subgroup enumeration.
Packed witnesses reduce the support table to 738,197,504 bytes. Rank-aware
collection reaches augmented rank 95 with exactly 95 relations, and ordered,
pipelined support expansion plus fixed/fused PCLMUL arithmetic preserves every
relation hash and group check. On the same public point
`Q=(2565091273463387,5885236316843894)`, fused direct takes 6.039837 wall
seconds / 11.269810 core-seconds / 1,051,987,968 B RSS versus 5.106080 seconds
for signed-Frobenius rho: a 1.182872 loss. Fresh build plus direct is 19.474746
times rho.

Stage 72 then retains the optimized unrelated public hash target. It constructs
no scalar, derives all factor-base logs from 95 relations, and has direct IC and
rho independently recover `7892094459170` with `[d]G=Q`. Unknown-scalar direct
takes 6.571510 seconds and remains 1.467611 times rho; fresh build plus direct
is 25.359925 times rho. The corrected Linux single-core receipt pins parent and
child to CPU 0: direct takes 10.945707 wall / 8.576129 core-seconds versus rho
5.490656 wall / 4.334277 core-seconds, a 1.993515 wall and 1.978676 CPU loss.

Stage 89 retains the post-optimization selected stack. It keeps generic initial
pair sums and generic Frobenius support expansion while selecting split Bloom
indices, insert-hash reuse, shared dual-sign denominators, 16-byte compact query
scratch, combined n=53 dispatch, and the exact direct-PCLMUL pair-query batch.
The specialized query batch improves collection by 17.66--19.31 percent across
three same-binary A/B runs while preserving every relation hash and counter.
Five selected-stack process ratios against the same signed-Frobenius rho target
are 0.946483, 1.019315, 1.094974, 0.932363, and 1.018776: two wins, three losses,
and a 1.018776 median loss. The current-head process reports 5.227982 wall
seconds / 9.193539 core-seconds / 1,053,990,912 B RSS versus 5.131633 seconds
rho. Fresh build plus direct remains 20.262764 times rho. The optimized public
unknown-scalar process takes 5.629163 seconds and remains 1.279112 times rho.

Stage 92 refreshes the forced-single-core comparison on that selected stack.
Linux affinity restricts both parent and child to CPU 0, and the direct arm uses
one query thread with the exact specialized pair batch. Direct takes 9.405716
wall seconds / 6.868316 core-seconds / 863,047,680 B RSS versus 6.945058 wall
seconds for rho: a 1.354303 wall and 1.340180 core loss. Fresh build plus direct
remains 14.691007 times rho.

Stage 99 archives the subsequent optimization chain. Direct low/high
coordinate windows replace the mixed Bloom hash; candidates that already pass
the filter skip a redundant second filter read; and an Itoh--Tsujii addition
chain reduces each degree-53 inversion from 105 field products to 59. The
matched Stage 94 arm reduces wall by 5.09 percent and core cost by 6.96 percent
with identical relation hashes. Stage 95 retains width 4096. Stage 96 rejects a
one-word blocked filter despite reducing exact misses by 69.27 percent because
its query and process wall regress. Stage 97 then records five online direct
wins, but its 0.969444 median misses the predeclared 0.95 gate. Stage 98 reduces
repeated evidence output by 98.13 percent and wall by 1.20 percent while
retaining the base hash, all relation hashes, and every mathematical and
reference validation. Host-to-host rho variation remains material.

Stage 108 archives the four-shard support-table and direct-routing chain.
Sharding reduces matched process wall by 17.32 percent and support setup by
32.55 percent. Replacing the per-query mixed shard hash with xor routing then
wins all five host-identified comparisons, with 0.904540 median wall and
0.912199 median core ratios against mixed routing. The selected EPYC 9V74
direct/rho panel ratios are 0.847095, 0.858937, 0.839190, 0.829081, and
0.838930. The public hash-derived unknown-scalar run takes 4.307977 seconds
versus 4.362556 rho and reconstructs the target without a supplied scalar or
constructed base logs. The exact support still occupies 738,197,504 bytes;
selected median direct core is 2.261565 times rho and fresh build plus direct
is 17.636692 times rho.

Stage 130 adds the distinct current `ic workflow` route on public hash seed
53001.  Its factor base is the target-independent subgroup-orbit recipe with
36,464 points and 344 projected columns; no factor-base logs or target scalar
are constructed.  The selected one-unit stream has 588 relations with canonical
SHA-256 `d8e18605b4a27f307101ecbb5ba9973b50fc01ba272d6f39d216c4fd6a615d37`
and recovers `7892094459170`.  Predicate construction is charged once and emits
344 cofactor multiplications, 1,896,128 derived-coordinate squarings, 18,232
negations and 72,928 canonical-orbit coordinate squarings.  Five matched
one-worker runs have median IC 5.195076 seconds versus 1.580426 rho, a 3.29
times loss; five default-thread runs have median IC 0.891978 versus 1.580089
rho, a 1.777988 rho/IC wall advantage.  Whole-process CPU remains dearer, so
this does not close the full-cost gate.

Stage 131 refreshes one target from the fixed-algebraic `n=41` holdout under
the current workflow.  The public recipe saturates the Frobenius union generated
by masks `[1,2,4,8,16,32]` with two-torsion; it materialises 4,759 points and
29 projected columns from zero target samples or scalar labels.  Fourteen
8,192-probe units yield the same 35 relations in all runs, canonical SHA-256
`8e51f7adab62f4bc351bcb35225fd6522d7cd56eeeac4f162bc00ef6c9b8058a`.
Public hash seed 41301 has no constructed scalar; relation-derived logs recover
`281099696942` and verify `[d]G=Q`.  Five default-thread processes have median
full IC 0.503511 seconds versus 0.118236 rho (4.258367 times slower), 3.336178
whole-process core-seconds and 13,221,888 B RSS.  Five one-worker processes have
median IC 2.829142 seconds versus 0.118650 rho (23.817013 times slower), 2.823121
whole-process core-seconds and 10,944,512 B RSS.  The twelve-process retained
series charges 21.784403 sequential wall-seconds, 37.073066 core-seconds and
14,352,384 B maximum RSS.  This current full-cost result fails the rho gate even
though the older five-target experiment had a fast online descent after shared
precomputation.

Stage 132 adds the exact standard Phase-B `K_1/F_(2^59)` base to the workflow:
the polynomial span `x = sum(c_i z^i), 0 <= i < 9` has 512 abscissae, 483
rational points, 242 negation columns and 231 public projected columns.  It is
not Frobenius-closed, so the implementation refuses the folded table and uses
an exact 116,886-entry compact table.  Public hash seed 59001 constructs no
scalar.  Eight 150,000-probe units scan 307,200,000 summands and return zero
relations; logs records zero solve attempts and a censored rank failure.  The
full capped process costs 18.155840 wall seconds / 78.395618 core-seconds /
14,696,448 B RSS.  There are at most `C(485,3)=18,896,570` unordered triples
against subgroup order 25,179,555,920,633, so even ideal uniform coverage is at
most `7.5047e-7`; 1.2M probes have expected yield at most 0.9006 and 231
relations need at least 307.8M probes before window, collision, dependency and
rank losses.  Five metered same-target rho runs all recover `17861472351607`;
median rho is 3.728387 wall / 3.719220 core-seconds and maximum RSS is
88,522,752 B.  The spent IC cap is already 4.87 times median rho wall and 21.08
times median rho core without reaching one relation.

Stage 133 makes the cofactor obstruction exact across standard widths.  In the
small cyclic quotient of order 22,894, only `4.354e-5` to `4.593e-5` of
unordered triples cancel for every measured `ell=9..16`, close to the exact
`1/22894` penalty.  The best class-aware probe lower bounds among bases fitting
the 4 GiB compact-table budget are 25.73B at `ell=13`, 6.45B at `ell=14`, and
1.582B at `ell=15`; the last table already occupies 2,709,116,566 bytes.  At
`ell=16`, the lower bound falls to 399.24M probes but 2,157,161,086 pair entries
need 9,975,660,343 bytes, over budget.  A charged `ell=14` pilot builds the
132,999,895-entry table and scans 102.4M summands in 85.949930 wall /
146.766002 core-seconds / 896,630,784 B RSS with zero relations, consistent with
its at-most 0.125 expected yield.  The full frontier charges 16 processes,
119.367630 sequential wall-seconds, 185.617302 core-seconds and 896,630,784 B
maximum RSS.  This is a finite no-go result for materialized standard bases
under this budget, not for implicit bases or index calculus in general.

Stage 134 removes the standard base's cofactor-class obstruction without using
secret labels: mapping every `ell=14` parent point through public `[h]` produces
16,308 prime-subgroup points and 8,094 projected columns.  Each complete run
collects 26,796 verified relations over 5.3 million probes and 5.4272 billion
summand scans.  Sparse filtering reduces this to an 802-column core; block
Wiedemann uses 611 products, all 8,094 logs are group-certified, and both modes
recover `17861472351607` for public hash seed 59001.  The default-thread workflow
uses 132.782432 IC wall seconds / 1,387.731459 whole-process core-seconds /
895,205,376 B peak RSS versus 0.553665 seconds rho.  One worker uses
1,219.842388 IC wall seconds / 1,190.875026 whole-process core-seconds /
850,935,808 B peak RSS versus 0.546318 seconds rho.  Online descent takes
0.024512 and 0.023714 seconds respectively, but those online ratios exclude the
large precomputation and do not satisfy the full-cost gate.

Stage 135 selects `ell=15` for this one target.  The public cofactor image has
32,934 points, 16,467 negation orbits and 16,344 projected columns; its exact
compact table stores 542,340,645 pairs.  Both complete modes collect the same
54,749 relations, canonical SHA-256
`e2004ac6e81979caf984e4fa745dc1a1ee99d13892a3f60615a5662179e3ad99`,
over 2.6 million probes and 2.6624 billion scans, then recover and verify the
same scalar as Stage 134.  Against `ell=14`, default IC wall falls from
132.782432 to 126.398334 seconds and whole CPU from 1,387.731459 to
1,194.982194 seconds.  One-worker IC wall falls from 1,219.842388 to
1,009.627483 seconds and whole CPU from 1,190.875026 to 1,000.180872 seconds.
Peak RSS grows by 4.888 times default and 4.798 times one-worker.  The selected
full-cost ratios remain 229.710 and 1,811.815 times rho.

Stage 136 runs two equal-scan neighboring windows and two repetitions of a
source-pinned buffer-reuse candidate.  Window 512 reaches 95.52 percent and
window 2048 reaches 91.47 percent of the selected window-1024 relation
throughput.  Reusing window scratch buffers preserves the 2,157-relation pilot
hash but takes 1.1235 and 1.1725 times the selected relation-unit wall.  The four
fresh rejected processes charge 115.333875 sequential wall-seconds,
1,078.171690 core-seconds and 4,378,951,680 B peak RSS.  The selected Stage 135
source remains unchanged.

Stage 137 stores one packed pair witness beside each of the 542,340,645
compact rests under an 8 GiB ceiling.  Exact group re-addition rejects truncated
rest collisions, and sorting exact candidates makes the parallel build
deterministic.  Default IC falls from 126.398334 to 85.747666 seconds and whole
CPU from 1,194.982194 to 853.002930 seconds, with RSS increasing from
4,375,724,032 to 6,565,560,320 B.  One-worker IC falls from 1,009.627483 to
724.574469 seconds and CPU from 1,000.180872 to 720.258242 seconds, with RSS
increasing from 4,082,794,496 to 6,253,510,656 B.  Both modes reproduce canonical
relation SHA-256
`e2004ac6e81979caf984e4fa745dc1a1ee99d13892a3f60615a5662179e3ad99`
and recover the same verified unknown scalar.  Full cost remains 152.615 and
1,329.575 times rho.

Stage 138 caches 542,340,645 packed keys (4,338,725,160 temporary bytes) from
the counting pass.  Default build wall falls from 21.389824 to 11.712773 seconds
and full IC from 85.747666 to 76.342293 seconds; one-worker build falls from
198.245643 to 101.933057 seconds and full IC from 724.574469 to 648.780325
seconds.  Both full modes preserve the 54,749-relation hash and verified scalar.
Peak RSS rises to 10,072,309,760 B default and 10,044,899,328 B one-worker.
Full cost remains 140.153 and 1,158.501 times rho.

Stage 139 allocates a 1,073,741,824-byte witnessed presence filter in place of
the selected 536,870,912-byte filter.  Both 2,157-relation pilots preserve
canonical SHA-256
`5d486165d6e822e83795c876b7eb8bf62d2f1a8e281276d18aa0fd0667468054`,
but relation-unit wall ratios are 1.0824 and 1.0659.  The two rejected processes
charge 35.782623 wall-seconds, 328.790212 core-seconds and 10,076,028,928 B peak
RSS.  Stage 138 remains selected.

Stage 140 finds that uniform units 22, 23, 24 and 25 leave 5, 4, 1 and 1
columns uncovered.  Missing column 15726 has public factor-base point index
21212.  Fixing that point and performing 10,000 direct pair lookups yields the
exact relation `(a, points) = (6472497976388, [3990,21212,22115])`, canonical
SHA-256
`c83ba65ae1727ae5aac4b8506a7ffe2d96e81fb06805b3bc32ea971fb2ffaadc`.
It completes a 52,635-row certified log system and recovers the same unknown
scalar.  Matched default IC falls from 182.728457 to 129.954516 seconds and
whole CPU from 791.419437 to 772.360819 seconds.  Matched one-worker IC falls
from 636.707717 to 630.392787 seconds and CPU from 630.367697 to 626.276231
seconds.  The retained Stage-140 campaign charges 23 processes, 2,134.246808
wall-seconds, 6,492.868869 core-seconds and 10,088,611,840 B peak RSS; one
preliminary 15.614915-second process lacks CPU/RSS and remains an explicit gap.

Stage 141 adds a public sparse-rank fallback and selects a 15-unit stream with
31,727 uniform relations from 1.5 million probes and 71 targeted relations from
2.95 million direct lookups.  Its canonical combined SHA-256 is
`32cb1609d40513a12750254ba291fb2a8a262a6eaeb05e598b28212a53fbf3aa`.
Default full IC is 58.027800 seconds / 526.896900 core-seconds / 10,072,866,816
B RSS; one-worker full IC is 431.866038 seconds / 428.938946 core-seconds /
10,043,244,544 B RSS.  Both recover and verify the same unknown scalar.  The
incremental frontier and selected campaign charges 19 processes, 996.612125
sequential wall-seconds, 3,748.178032 core-seconds and 10,091,528,192 B peak
RSS.  Full cost still loses to same-target signed-Frobenius rho by 103.215 and
773.321 times.

Stage 142 adds a rank-disabled 15-unit control.  Its 295 forced-column attempts,
71 targeted relations, sparse report and combined relation hash exactly match
the selected stream.  The attempt counts by round are all below 64, proving the
selected execution always followed the raw-uncovered-column branch and never
entered the rank fallback.  The rank phase is retained as a lower-prefix
frontier mechanism, while the selected speedup is attributed to uncovered-column
targeting after 15 uniform units.  The two correction processes add 36.873650
wall-seconds, 194.727064 core-seconds and 10,085,711,872 B peak RSS.

Stage 143 evaluates a deterministic first-hit coverage tail.  Thirty-two
64-probe runs execute in parallel per batch; after the first successful batch,
only its earliest exact relation is retained.  Candidate lookups fall from
2,950,000 to 2,713,712 and targeted relations from 71 to 58.  Both matched
candidate runs emit combined relation SHA-256
`dcbf0aad8eabce6f22241801f70ade6f4d948e2d1bc67d39a5dd1d2fd140e76c`,
derive every log, and recover the same public-target scalar.  Across the
default-thread ABBA pair, median whole wall rises from 78.506868 to 81.833509
seconds and median CPU from 540.796846 to 541.830480 core-seconds.  The nine
charged processes consume 633.598448 sequential wall-seconds, 3,806.440475
core-seconds and 10,080,059,392 B peak RSS.  The source patch is archived but
not applied to the selected code.

Stage 144 runs the 14-unit, 64-column rank tail with a sparse solve after every
four rank rounds.  It retains 29,562 uniform and 420 targeted relations, uses
20,030,000 direct pair lookups, and finishes 27 linear-algebra attempts in
4.934907 seconds.  All logs and the same unknown scalar verify, but full IC is
122.781740 seconds / 676.071984 core-seconds / 9,350,201,344 B RSS, versus the
selected 58.027800-second receipt.  The one fully charged process is retained,
and the cadence patch is archived without entering the selected source.

Stage 145 keeps the 14-unit mathematical policy fixed while moving independent
fixed columns onto the Rayon pool.  Both the nested-Rayon pilot and the final
single-layer serial-column candidate emit the same 29,948 relations and recover
the same unknown scalar.  Their full IC times are 118.335409 and 138.522551
seconds, with 651.606000 and 658.228767 core-seconds.  Together they charge
258.133232 sequential wall-seconds, 1,309.834767 core-seconds and
10,082,500,608 B peak RSS.  The 10 GB random-access table saturates memory
bandwidth, so the final patch is archived and rejected.

Stage 146 records original-column support when block Wiedemann returns a
homogeneous kernel vector with zero homogenising coordinate.  The prefix-14
candidate receives no such vector in any failed attempt, falls back to the same
global-incidence columns, and exactly repeats combined relation SHA-256
`1b97d125b7f08a1404ee14bcc794e7e008aa2729482fa98d1ffe61a380356b4e`.
The complete workflow takes 116.732654 seconds of IC wall, 644.323580
core-seconds and 10,077,847,552 B RSS.  The one process is charged and the
inactive diagnostic policy is rejected.

Stage 147 measures prefetch lookaheads 32 and 64 over four source-pinned runs
per arm.  Every process builds the full witnessed table, executes one
100,000-probe / 102.4-million-scan unit, and emits the same 2,157-relation hash.
Candidate median unit wall rises from 3.523893 to 4.653538 seconds and unit CPU
from 22.564842 to 23.887295 core-seconds; median whole wall rises from
30.139010 to 33.767785 seconds.  The eight processes charge 265.845801
sequential wall-seconds, 1,291.854067 core-seconds and 10,080,190,464 B peak
RSS.  The one-line patch is archived and rejected.

Stage 148 measures lookahead 16 against 32 in eight one-unit processes and four
complete unknown-scalar workflows.  All pilot relations and all 31,798-row full
streams match exactly.  In the full panel, candidate median whole wall falls
from 76.843062 to 75.030543 seconds and CPU from 508.435756 to 500.166199
core-seconds, but the only changed path—relation-unit wall—rises from 40.041194
to 40.490959 seconds.  The favorable aggregate comes from unchanged build and
logs stages.  The twelve processes charge 499.765728 sequential wall-seconds,
3,287.061575 core-seconds and 10,080,911,360 B peak RSS; the patch is rejected.

Stage 149 removes a duplicate filter hash and random filter load from admitted
window probes.  The eight-run one-unit panel improves median unit wall 14.45
percent.  A helper variant then regresses verification, and the minimal boolean
variant is extended to six complete runs per arm.  Its median whole CPU improves
0.58 percent and unit CPU 0.47 percent, but median unit wall regresses 1.22
percent and the paired unit-wall ratio median exceeds one.  All full runs emit
the same 31,798 relations and verified scalar.  The 24 processes charge
1,347.433898 sequential wall-seconds, 9,477.238660 core-seconds and
10,084,401,152 B peak RSS.  Both variants are rejected.

Stage 150 keeps the 512 MiB filter allocation but places one independent bit in
each half.  The expected false-admission probability falls, while every stored
pair pays two atomic bit sets and probes may load two filter words.  Four
source-pinned processes emit the same 2,157 relations.  Candidate median unit
wall rises from 1.932228 to 1.980474 seconds, unit CPU from 23.218455 to
24.615690 core-seconds, and build CPU from 139.615756 to 148.997076
core-seconds.  They charge 72.543322 wall-seconds, 696.357834 core-seconds and
10,075,766,784 B peak RSS.  The no-extra-memory split is rejected.

Stage 151 computes one filter hash per key block and reuses it for both
lookahead prefetch and admission.  Four source-pinned one-unit processes emit
the same 2,157 relations.  Candidate median unit wall rises from 1.733608 to
1.743590 seconds and unit CPU from 22.506273 to 22.575571 core-seconds.  They
charge 63.849557 wall-seconds, 662.355731 core-seconds and 10,072,932,352 B
peak RSS.  The 8 KiB block scratch costs more than the avoided hash and the
candidate is rejected.

Stage 152 first checks the presence filter for a 1,024-key block, stores the
admitted offsets in order, prefetches their random `bucket_start` entries, then
consumes that FIFO without changing witness order.  The exact relation verifier
is pinned out of line to isolate it from query-loop code layout.  Eight one-unit
pilots, sixteen complete default-thread processes across the preliminary and
selected variants, and a matched one-worker pair retain one exact relation hash
and scalar.  In the selected four-run-per-arm default panel, median unit wall
falls 19.51 percent, unit CPU 6.39 percent, full IC wall 12.46 percent and whole
CPU 4.36 percent; every paired unit comparison wins.  One-worker IC wall falls
from 491.894708 to 437.950836 seconds and CPU from 461.794250 to 425.370001
core-seconds.  The 26-process campaign charges 2,659.078039 wall-seconds,
10,455.178415 core-seconds and 10,084,040,704 B peak RSS.  Full cost still loses
to rho by roughly 140 and 774 times.

Stage 153 adds a middle pass over admitted offsets: after `bucket_start` has
been prefetched, it reads the run start and prefetches the first compact rest.
Eight one-unit processes improve median unit wall 30.25 percent and CPU 4.18
percent.  Four complete workflows reverse that result: candidate median unit
wall rises 19.53 percent, full IC wall 12.87 percent, and verification 199.10
percent.  All hashes and scalars still match.  The twelve processes charge
717.712397 wall-seconds, 3,307.549828 core-seconds and 10,084,892,672 B peak
RSS.  The patch is archived and rejected.

Stage 154 changes only the selected window loop's internal block constant from
1,024 to 512.  Eight source-pinned one-unit processes emit the same 2,157
relations.  Candidate median unit wall falls from 2.904745 to 2.706591 seconds,
while unit CPU rises from 22.235606 to 22.753293 core-seconds.  The panel charges
201.720277 wall-seconds, 1,318.438451 core-seconds and 10,085,335,040 B peak RSS.
The wall-for-work trade is rejected.

Stage 155 changes the same constant from 1,024 to 2,048.  Four source-pinned
processes emit the same 2,157 relations.  Candidate median unit wall rises from
2.512955 to 2.648489 seconds and CPU from 22.704948 to 22.920353 core-seconds.
They charge 87.236383 wall-seconds, 668.786782 core-seconds and 10,084,302,848 B
peak RSS.  The larger cache footprint is rejected.

Stage 156 changes only the selected FIFO presence-filter lookahead from 32 to
64.  Eight processes emit the same 2,157 relations.  Candidate median unit wall
rises from 2.512537 to 2.584736 seconds, while unit CPU falls from 22.580665 to
22.425630 core-seconds.  They charge 182.501493 wall-seconds, 1,331.369772
core-seconds and 10,082,418,688 B peak RSS.  The wall regression rejects it.

Stage 157 changes the selected FIFO presence-filter lookahead from 32 to 16.
Eight processes retain the same 2,157 relations.  Candidate median unit wall
falls from 2.562383 to 2.547881 seconds, while CPU rises from 22.743259 to
23.482047 core-seconds.  They charge 177.409112 wall-seconds, 1,353.602444
core-seconds and 10,080,960,512 B peak RSS.  The extra work rejects it.

Stage 158 adds an `already_admitted` flag to the existing lookup body and uses
it only from the selected FIFO consumer.  Four processes retain the same 2,157
relations.  Candidate median unit wall rises from 2.600805 to 2.652481 seconds
and CPU from 22.485871 to 23.728118 core-seconds.  They charge 93.739907
wall-seconds, 686.047697 core-seconds and 10,085,679,104 B peak RSS.  The fast
path is rejected.

Stage 159 lowers compact bucket bits from 26 to 25.  Four processes emit the
same 2,157 relations.  Candidate median unit wall changes from 2.628147 to
2.625274 seconds, while CPU rises from 22.421583 to 23.814860 core-seconds.
They charge 88.843412 wall-seconds, 672.866766 core-seconds and 10,080,976,896 B
peak RSS.  The longer compact runs reject the memory trade.

Stage 160 raises compact bucket bits from 26 to 27.  Four processes emit the
same 2,157 relations.  Candidate median unit wall rises from 2.589146 to
2.635581 seconds, CPU falls from 23.064888 to 22.614644 core-seconds, and median
peak RSS rises from 9,715,990,528 to 10,618,896,384 B.  They charge 90.389257
wall-seconds and 677.655218 core-seconds.  The wall-and-memory trade is rejected.

Stage 161 doubles the witnessed presence filter from 536,870,912 to
1,073,741,824 bytes on the selected FIFO stack.  Four processes retain the same
2,157 relations.  Candidate median unit wall rises from 2.497084 to 2.648392
seconds, unit CPU from 22.804185 to 23.158461 core-seconds, build wall from
16.090470 to 18.134880 seconds, and median RSS from 9,729,105,920 to
10,080,403,456 B.  They charge 91.836085 wall-seconds and 674.917285
core-seconds.  The wider filter is rejected again.

## The n=53 crossover against a batched rho

The Stage 108 online-wall crossing was measured against `koblitz_rho_fixture`.
That control is correct but slow. It inverts once per step, canonicalises both
coordinates by 52 squarings each, stores every point, and runs on one thread:
about 9.5 µs a step. On 2026-09-22, `suite/examples/koblitz_batched_rho.rs` in
`aburan28/cryptanalysis` (PR #54) ran the same signed-Frobenius walk with one
inversion shared across 64 walks, the normal-basis orbit key, and only
distinguished points stored. It was run on the same target (`d=476811900269`)
and the same host as the selected direct route, and every run recovered the
scalar and checked `[d]G=Q`:

| Arm | Runs | Median wall |
|:--|--:|--:|
| Selected direct (Stage 106/108 flags), 4 threads | 3 | 17.3 s (about 16 s CPU) |
| `koblitz_rho_fixture`, 1 thread | 5 | 3.03 s |
| batched rho, 1 thread | 7 | 0.49 s |
| batched rho, 4 threads | 7 | 0.20 s |

The host is an Apple M4 Pro, which has none of the x86 PCLMUL paths both
fixtures specialise, so these are not EPYC gate-host ratios and the gate host
needs its own rerun. The operation counts do not depend on the host: the direct
route spends 46.3M support queries against about 0.56M rho steps. The batched
walk still favours index calculus in one respect, because it moves `y` onto the
orbit representative by up to `n-1` squarings.

The gap widens with degree. A run spends about `r/(n·|F|)` queries (40M
predicted at `n=53`, 46.3M measured), and the base rule `|F|^3 ≈ 6r·η` makes
that grow as `r^{2/3}` against rho's `r^{1/2}`. At ECC2K-130 that is about 2^80.6
queries against 2^60.9 rho steps. Folding the support table by the signed
Frobenius orbit (738 MB to about 15 MB) would pass the memory requirement and
change neither verdict. Per-run records are in
`suite/docs/ic/runs/koblitz-n53-rho-control-20260922.json` in that repository.

The narrow supported conclusion is unchanged: this is strong internal
engineering and finite public toy-research evidence. The known SAT-based
point-decomposition and Frobenius-invariant-factor-base literature remain the
prior-art baseline. The finite `n=41` online descent is faster after
precomputation. The selected `n=53` panel passes the predeclared online wall
threshold on the recorded host only against the fixture rho control. Against a
batched signed-Frobenius rho it loses by one to two orders of magnitude, and its
query count grows faster than rho's with degree. Its memory, core, and
fresh-build costs remain above rho. Licensed Magma and unaffiliated
reproduction/novelty review remain open. This is not a new Koblitz
index-calculus SOTA result.

## Parallel-branch native-F4 history (Stages 159–170)

These stage numbers were assigned on the native-F4 branch before `main` used the same numbers for the n=59 collector chain. The collector Stage-159 audit is preserved under `stage-159-bucket25-current-gate-audit-20260922`; the text below retains the native-F4 numbering required by its sealed Stage-162–170 predecessor chain.

The native chain includes Stages 161 and 162 trusted-mask hashing and grouped
critical-pair selection, plus Stage 163 shape-selected block-4 M4RI elimination,
followed by the dense pipeline, target-independent ordering, product
cancellation, deterministic parallel batches, fresh-target holdout and
parallel fixed-X1 construction documented below.

## Stage 159 native-F4 single-target supplement

Stage 159 runs the repository's full Boolean `f4-f2` implementation on one
authenticated `n=59, ell=9, m=3` Phase-B target. The solver receives no truth
class, witness, discrete-log label, or target-subgroup enumeration. A direct
three-summand S4 expansion exceeded the watchdog at 5.25--6.31 GB RSS in two
censored attempts. The retained formulation fixes `X1`, runs full F4 over the
remaining 18 Boolean variables, and checks every returned root by exact curve
lifting and group addition.

The first complete arm took 205.702959 wall seconds, 204.982191 core-seconds,
and 950,878,208 bytes RSS. Skipping only `X1` masks with no rational curve lift
reduced it to 100.179380 wall seconds, 99.920352 core-seconds, 921,681,920
bytes RSS, 65 F4 calls, and 87,513,949,370 elimination word XORs while
returning the same exact witness. This is a 2.05x solver-stage engineering
gain. A trace-homomorphism equation was retained as a rejection: it raised
wall by 18.2 percent and word XORs by 58.8 percent despite lowering RSS.

On the same local host and exact target, direct MITM completed in 2.999661
seconds / 2.991745 core-seconds / 45,842,432 bytes, so selected F4 remains
33.40x slower by wall and CPU and 20.11x larger by RSS. Native XOR SAT reached
its 100,000-conflict cap in 8.906472 seconds. The authenticated Linux receipts
for this target leave WDSat and CryptoMiniSat censored at their 120-second
watchdogs. Native F4 is a distinct implementation and formulation; it does
not satisfy the still-missing licensed Magma F4 gate. Its conflict field is
`null` because F4 does not expose SAT conflicts.

The Stage 159 measured research lower bound, including five distinct clean
builds, six F4 attempts, local controls, failed builds, and dependency-fetch
receipts, is 1,596.064123 sequential wall-seconds, 1,565.839519 core-seconds,
and 6,310,576,128 bytes maximum process RSS. Complete campaign cost remains
`null` because some early incremental checks lacked an outer meter. The full
160-input native-F4 panel, licensed Magma, end-to-end IC/rho cost, and
unaffiliated review remain open. The canonical result is
[`stage-159-native-f4-single-target-20260922/result.json`](stage-159-native-f4-single-target-20260922/result.json).

## Stage 160 fixed-X1 constant-linear construction

Stage 160 preserves the exact Stage 159 factor base, blind target, X1 schedule,
65 F4 calls, 3,864,601 Boolean terms, equation fingerprint, 87,513,949,370
elimination word XORs, algebraic roots, and group-verified witness. It changes
only how each fixed-X1 S4 system is constructed. The shared symbolic `X2*X3`
product is formed once, while fixed `X1` and target powers act through exact
linear multiplication maps in the polynomial basis. A unit gate compares the
complete specialized and generic polynomial vectors across 60 field/value
cases.

The clean blind run reduces fixed-X1 construction from 26.323585 seconds to
0.973558 seconds, a 27.04x construction improvement. Whole-target F4 wall
falls from 100.179380 to 73.757990 seconds and core cost from 99.920352 to
73.585054 seconds, a 1.36x wall improvement. Peak RSS falls slightly from
921,681,920 to 910,934,016 bytes. The clean one-job build plus custody run is
241.262254 wall seconds / 236.554363 core-seconds / 1,683,472,384 bytes peak
RSS.

Direct MITM remains 24.59x faster by same-host wall, 24.60x cheaper by CPU,
and 19.87x smaller by RSS. Four other metered ideas are retained as rejected:
global matrix ordering, trailing-zero row trimming, a hashed reducer index,
and a dense reducer index. Through this stage the measured campaign lower
bound is 2,319.762066 sequential wall-seconds, 2,281.647854 core-seconds, and
6,310,576,128 bytes maximum process RSS. Complete campaign cost remains
`null` because the incremental compiles and focused tests were not outer
metered.

This remains one finite toy-PDP target. Licensed Magma, the full native-F4
panel, end-to-end IC/rho cost, and unaffiliated reproduction remain open. The
canonical additive result is
[`stage-160-fixed-x1-constant-specialisation-20260922/result.json`](stage-160-fixed-x1-constant-specialisation-20260922/result.json).

## Stage 161 trusted-mask hashing

A charged 20-second sample of the Stage 160 binary collected 3,334
top-of-stack samples: 1,143 in dense echelon elimination, 533 in critical-pair
installation, 238 in column packing, 225 directly in SipHash writes, and 118
in queued-pair retention. Stage 161 replaces SipHash only for internal maps
and sets keyed by solver-constructed `u64` monomial masks. The deterministic
SplitMix64 hash changes bucket placement only; hashbrown still resolves every
collision by exact key equality. Tuple keys and other data retain their prior
hashers.

The clean blind run preserves the equation fingerprint, 3,864,601 terms, all
65 F4 calls, matrix dimensions, pair counters, 87,513,949,370 elimination word
XORs, roots, and group-verified witness. Symbolic/matrix build falls from
30.523303 to 17.204156 seconds, a 1.77x phase improvement. Whole-target F4
wall falls from 73.757990 to 60.645599 seconds and core cost from 73.585054 to
60.498754 core-seconds, a 1.22x improvement. Peak RSS rises 5.85 percent to
964,263,936 bytes.

Same-host direct MITM remains 20.22x faster by wall and CPU and 21.03x smaller
by RSS. The clean one-job build plus custody run costs 228.347900 wall seconds,
224.535510 core-seconds, and 1,692,057,600 bytes peak RSS. Through this stage
the measured campaign lower bound is 2,684.534288 sequential wall-seconds,
2,641.558033 core-seconds, and 6,310,576,128 bytes maximum process RSS.
Complete cost remains `null` because one sandbox-denied profiler attempt was
interrupted before its receipt and incremental compiles and focused tests were
not all outer metered.

The result remains one finite toy-PDP target. Licensed Magma, the full
native-F4 panel, end-to-end IC/rho cost, and unaffiliated reproduction remain
open. The canonical additive result is
[`stage-161-fast-mask-hash-20260923/result.json`](stage-161-fast-mask-hash-20260923/result.json).

## Stage 162 grouped critical-pair selection

A charged post-hash profile put 695 of 3,327 top-of-stack samples in
`State::insert`. Binary disassembly placed its dominant samples in the
quadratic scan that tests whether any other new-pair LCM divides the current
LCM. Stage 162 groups equal LCMs and finds proper divisor LCMs by exact submask
lookup. It then emits the same lowest-index duplicate representative in the
same pair order as the previous UPDATE implementation. A 3,600-case randomized
test compares selected pairs, order, and chain/product counters directly with
that quadratic reference.

The clean blind run preserves every equation, term, F4 call, matrix dimension,
pair counter, 87,513,949,370 elimination word XORs, root, and exact witness.
Whole-target wall falls from 60.645599 to 57.862894 seconds and core cost from
60.498754 to 57.736423 core-seconds, a 1.048x improvement. Peak RSS changes by
0.40 percent to 968,146,944 bytes.

Same-host direct MITM remains 19.29x faster by wall, 19.30x cheaper by CPU,
and 21.12x smaller by RSS. The clean build plus custody run costs 225.353725
wall seconds, 221.602144 core-seconds, and 1,692,352,512 bytes peak RSS. The
measured campaign lower bound through Stage 162 is 3,029.348391 sequential
wall-seconds, 2,981.632960 core-seconds, and 6,310,576,128 bytes maximum process
RSS. Complete campaign cost remains `null`.

This remains a finite one-target solver-stage result. Licensed Magma, the full
native-F4 panel, end-to-end IC/rho cost, and unaffiliated reproduction remain
open. The canonical additive result is
[`stage-162-grouped-pair-selection-20260923/result.json`](stage-162-grouped-pair-selection-20260923/result.json).

## Stage 163 shape-selected block-4 M4RI

Stage 163 routes matrices with at least 128 rows and 256 columns and at most
four times as many columns as rows through block-4 Method of Four Russians
elimination. Other shapes retain streaming elimination, and
`PQ_F4_DISABLE_M4RI=1` provides a same-binary control. Four large synthetic
matrix shapes across word boundaries return identical pivot-column sets and
canonical row spaces under the two kernels.

The same-binary control takes 58.370052 wall seconds / 57.840795 core-seconds;
M4RI takes 54.538915 / 54.413164, a 1.070x wall and 1.063x CPU improvement.
The clean selected run takes 54.773779 wall seconds, 54.662604 core-seconds,
and 915,668,992 bytes peak RSS. It routes 196 matrices through M4RI and reduces
counted elimination-and-table work from 87,513,949,370 to 45,879,309,338 word
XORs. The equation fingerprint, target classification, roots, and exact
group-verified witness are unchanged; a few intermediate basis-path row and
reducer counters differ slightly and are retained rather than called
identical.

Direct MITM remains 18.26x faster by wall, 18.27x cheaper by CPU, and 19.97x
smaller by RSS. The selected clean build plus run costs 223.248497 wall
seconds, 218.415925 core-seconds, and 1,689,485,312 bytes peak RSS. The measured
campaign lower bound through Stage 163 is 3,589.954839 sequential wall-seconds,
3,532.297699 core-seconds, and 6,310,576,128 bytes maximum process RSS.
Complete campaign cost remains `null`.

The result remains a finite one-target solver-stage improvement. Licensed
Magma, a full native-F4 panel, end-to-end IC/rho cost, and unaffiliated
reproduction remain open. The canonical additive result is
[`stage-163-m4ri-echelons-20260923/result.json`](stage-163-m4ri-echelons-20260923/result.json).

## Stage 164 dense native-F4 pipeline

Stage 164 keeps the same algebraic factor base, blind `n=59, ell=9, m=3`
target, 65 fixed-X1 F4 calls, 3,864,601 Boolean terms, equation fingerprint,
two algebraic roots, and exact three-point curve witness. It adds dense
epoch-stamped LCM grouping on the 18-variable systems, exact restricted
submask covers, reusable pair-selection storage, indexed symbolic reducers,
dense monomial-to-column maps, order-preserving batch UPDATE, cache-local
M4RI tables, consecutive-pivot extraction, and trailing-zero table trimming.
Eleven focused tests include 3,600 randomized UPDATE comparisons, 230,400
indexed-reducer comparisons, randomized batch-vs-sequential state equality,
100,000 DegRevLex comparator comparisons, and exact large-matrix row-space
equality against streaming elimination.

The clean blind F4 process takes 36.450481 wall seconds, 36.319425
core-seconds, and 946,388,992 bytes peak RSS. It charges 27,264,366,281
elimination-and-table word XORs. The same clean binary takes 48.447073 seconds
with M4RI disabled and 38.698749 seconds with batch insertion disabled; all
arms retain the same equations, pair-criterion totals, roots, and witness.
Relative to Stage 163, selected wall improves by 1.503x and core cost by
1.505x. A separately metered lockfile dependency fetch costs 0.316796 wall
seconds / 0.102986 core-seconds, while dependency acquisition plus the clean
one-job build and custody run cost 214.044317 wall seconds / 208.873448
core-seconds / 1,689,780,224 bytes peak RSS.

Direct MITM remains 12.15x faster by wall, 12.14x cheaper by CPU, and 20.64x
smaller by RSS. The measured native-F4 campaign lower bound through Stage 164
is 6,464.947902 sequential wall-seconds, 6,383.990738 core-seconds, and
6,310,576,128 bytes maximum process RSS across 132 components. Complete cost
remains `null` because some development checks were not outer-metered.

This remains one finite toy-PDP solver-stage improvement. Licensed Magma, the
full native-F4 panel, end-to-end IC/rho cost, and unaffiliated reproduction and
novelty review remain open. The canonical additive result is
[`stage-164-native-f4-dense-20260923/result.json`](stage-164-native-f4-dense-20260923/result.json).

## Stage 165 target-independent sparse fixed-X1 order

Stage 165 changes only the order in which the symmetric direct-S4 solver fixes
the first factor-base coordinate. Coefficient masks are sorted by Hamming
weight and then numeric value. This permutation is defined solely by the
published polynomial-basis coefficients: it does not inspect the target,
truth class, witness, subgroup, discrete-log labels, or a known scalar.
`PQ_F4_X1_ORDER=ascending` retains the previous same-binary order, and a unit
test checks exact permutation and ordering for dimensions zero through twelve.

On the clean blind target, the selected process visits 85 masks, constructs
and completes 39 rational fixed-X1 F4 systems, and returns the same exact three
curve points as the ascending schedule. It takes 21.512089 wall seconds,
21.484441 core-seconds, 948,502,528 bytes peak RSS, and 16,821,055,616 charged
word XORs. The same clean binary in ascending order visits 138 masks, completes
65 systems, and takes 34.446599 wall / 34.412474 core-seconds, for a 1.601x
wall improvement. Three interleaved development pairs have 21.368530 versus
34.230726 wall-second medians, also 1.602x.

The schedule itself is target-independent, but it was selected post hoc after
inspecting this target. The result is therefore target-specific engineering,
not evidence that the order improves expected time over fresh or adversarial
targets. Direct MITM remains 7.17x faster by wall, 7.18x cheaper by CPU, and
20.69x smaller by RSS. The clean build plus custody run costs 200.076092 wall
seconds / 195.904177 core-seconds / 1,723,154,432 bytes peak RSS.

The measured native-F4 campaign lower bound through Stage 165 is 6,902.647655 sequential wall-seconds, 6,817.038349 core-seconds, and 6,310,576,128 bytes
maximum process RSS across 150 components. Complete cost remains `null`.

Licensed Magma, the full native-F4 panel, end-to-end IC/rho cost, and
unaffiliated reproduction and novelty review remain open. The canonical
additive result is
[`stage-165-x1-weight-order-20260923/result.json`](stage-165-x1-weight-order-20260923/result.json).

## Stage 166 dense F4 product cancellation

Stage 166 changes the construction of F4 rows after monomial multiplication.
For Boolean systems with at most twenty variables, one generation-tagged
`u32` array now records the parity of every mapped output mask. The first
occurrence stores a compact key whose natural order is descending degree and
ascending mask. The solver sorts only odd-parity survivors instead of sorting
all mapped terms and cancelling adjacent duplicates afterward.
`PQ_F4_DISABLE_DENSE_MUL=1` restores the exact previous path. A randomized
test compares the two products through twenty variables, including generation
wraparound, and the existing F4/Buchberger and M4RI/streaming certificates
continue to pass.

The clean blind process takes 19.887753 wall seconds, 19.858236 core-seconds,
and 900,169,728 bytes peak RSS. Across 818,171 monomial products it maps
651,502,431 input terms to 549,909,017 surviving terms, cancelling 101,593,414
terms before sorting. The exact clean disabled control takes 21.582929 wall /
21.552617 core-seconds / 939,032,576 bytes, for a 1.085x wall improvement.
Three interleaved development pairs have 20.107567 versus 21.526022
wall-second medians, a 1.071x improvement. Equation fingerprint, F4 calls,
matrix and algebraic counters outside the new product counters, roots, and the
exact curve-verified witness agree.

Direct MITM remains 6.63x faster by wall, 6.64x cheaper by CPU, and 19.64x
smaller by RSS. The clean build plus custody run costs 193.162520 wall seconds
/ 187.303839 core-seconds / 1,703,051,264 bytes peak RSS. The failed default
offline-cache build is retained and charged; it stopped at vendoring because
`redis 0.32.7` was absent and did not fall back to the network.

The measured native-F4 campaign lower bound through Stage 166 is 7,400.453979 sequential wall-seconds, 7,307.628306 core-seconds, and 6,310,576,128 bytes
maximum process RSS across 180 components. Complete cost remains `null`.

Licensed Magma, the full native-F4 panel, end-to-end IC/rho cost, and
unaffiliated reproduction and novelty review remain open. The canonical
additive result is
[`stage-166-dense-product-cancellation-20260923/result.json`](stage-166-dense-product-cancellation-20260923/result.json).

## Stage 167 deterministic parallel fixed-X1 batches

Stage 167 exploits only the independence already present between rational
fixed-X1 systems. The backend gathers the next deterministic batch in the
algebraic schedule, solves every launched system through a bounded Rayon
pool, waits for the complete batch, charges every result, and then inspects
results in schedule order. It does not receive target labels, a known scalar,
the witness, subgroup enumeration, or factor-base logarithms.

The clean thirteen-thread F4 process takes 5.589951 wall seconds, 36.302718
total core-seconds, and 3,033,563,136 bytes peak RSS. The same clean binary at
batch one takes 19.828459 wall seconds, 19.805719 single-core seconds, and
904,085,504 bytes. Parallel wall improves 3.547x while total CPU rises 1.833x
and RSS rises 3.355x. Three interleaved development pairs give 5.534149 versus
19.708587 wall-second medians, a 3.561x speedup. Both clean arms visit 85
masks, complete the same 39 F4 systems, charge 16,821,055,616 word XORs, have
the same equation fingerprint, roots and exact witness, and execute no
speculative system after the successful one.

Batch width thirteen was selected post hoc: this target's first valid
relation is in rational system 39, so three batches end exactly there. The
mechanism is label-blind, but this width result is target-specific rather than
fresh-target or expected-time evidence. Direct MITM remains 1.86x faster by
wall, 12.13x cheaper by CPU, and 66.17x smaller by RSS. Clean build plus the
parallel custody run costs 173.190335 wall seconds / 199.814249 core-seconds /
3,033,563,136 bytes peak RSS. The full-cost gate remains false.

The measured native-F4 campaign lower bound through Stage 167 is 7,948.233270 sequential wall-seconds, 8,144.184517 core-seconds, and 6,310,576,128 bytes
maximum process RSS across 215 components. Complete cost remains `null`.

Licensed Magma, the full native-F4 panel, fresh-target validation, end-to-end
IC/rho cost, and unaffiliated reproduction and novelty review remain open. The
canonical additive result is
[`stage-167-parallel-x1-batches-20260923/result.json`](stage-167-parallel-x1-batches-20260923/result.json).

## Stage 168 one deterministic batch of 39 systems

Stage 168 puts the same selected 39 rational fixed-X1 systems into one
deterministic batch, removing the two barriers in Stage 167. Twelve Rayon
workers execute the full batch dynamically; every system is charged before
the completed results are inspected in algebraic schedule order. The
factor base, equations, roots, work counters, and exact witness are unchanged.

The clean F4 process takes 4.987010 wall seconds, 28.179809 total core-seconds,
and 2,763,997,184 bytes peak RSS; `single_core_seconds` is `null`. Relative to
Stage 167's clean batch-13 arm, wall improves 1.121x, CPU falls to 0.776x, and
RSS falls to 0.911x. Three twelve-thread development runs have a 4.744985
wall-second median. The selected source and clean build are unchanged, and
the build cost is reused rather than charged twice.

Batch size 39 is exactly the already known successful prefix and was selected
post hoc. It is not a generic stopping rule: an unknown fresh target could
require additional batches or no relation within the cap. Direct MITM remains
1.66x faster by wall, 9.42x cheaper by CPU, and 60.29x smaller by RSS. Clean
build plus the selected custody run is 172.146798 wall seconds / 191.677868
core-seconds / 2,763,997,184 bytes peak RSS. The full-cost gate remains false.

The measured native-F4 campaign lower bound through Stage 168 is 7,988.745489 sequential wall-seconds, 8,381.199829 core-seconds, and 6,310,576,128 bytes
maximum process RSS across 223 components. Complete cost remains `null`.

Licensed Magma, the full native-F4 panel, fresh-target validation, end-to-end
IC/rho cost, and unaffiliated reproduction and novelty review remain open. The
canonical additive result is
[`stage-168-single-batch-39-20260923/result.json`](stage-168-single-batch-39-20260923/result.json).

## Stage 169 preregistered fresh-target holdout

Stage 169 selects the first presentation-order `n=59, ell=9, m=3` blind
instance after excluding the Stage-168 target. The selected ID, batch size 39,
twelve Rayon threads, budgets, no-retry policy, and publication of every
terminal were committed before execution or truth scoring. The selection did
not consult the target class or a witness.

Native F4 exhaustively visits all 512 coefficient masks, skips 270
non-rational values, constructs and completes 242 rational systems in seven
batches, processes 14,278 equations and 14,515,915 terms, and finds no roots.
It returns UNSAT in 20.539819 wall seconds / 169.094100 core-seconds /
2,872,541,184 bytes peak RSS while charging 99,199,976,264 word XORs. Only
after the run seal was frozen did scoring identify the target as
nondecomposable, producing one true negative.

On the same exported instance, direct MITM is also a true negative in 2.953129
wall seconds. Native XOR is inconclusive in 8.986403 seconds at 100,000
conflicts. WDSat and CryptoMiniSat reach 120.002468 and 120.007540-second
watchdogs. Licensed Magma remains unexecuted, and GGMP is not a same-instance
`n=59` construction. F4 is 6.96x slower than direct MITM by wall, 57.35x by
CPU, and 63.80x by RSS.

The fresh negative takes 4.119x the wall, 6.001x the CPU, and 6.205x the F4
systems of the post-hoc selected positive target. This rejects treating the
4.99-second Stage-168 result as generic fresh-target time. One holdout remains
finite evidence rather than an expected-time distribution or full panel.

The measured native-F4 campaign lower bound through Stage 169 is 8,261.489665 sequential wall-seconds, 8,801.978114 core-seconds, and 6,310,576,128 bytes
maximum process RSS across 229 components. Complete cost remains `null`.

Licensed Magma, the full native-F4 panel, end-to-end IC/rho cost, and
unaffiliated reproduction and novelty review remain open. The canonical
additive result is
[`stage-169-fresh-target-holdout-20260923/result.json`](stage-169-fresh-target-holdout-20260923/result.json).

## Stage 170 parallel fixed-X1 construction

Stage 170 keeps the Stage-169 algebraic factor base, target, Hamming-weight
mask order, batch size 39, twelve-thread pool, F4 implementation, and exact
verification path. It changes how each deterministic batch is prepared:
rational fixed-X1 S4 systems are constructed through an indexed Rayon
parallel iterator and collected in their original mask order before the
unchanged parallel F4 solve. `PQ_F4_DISABLE_PARALLEL_CONSTRUCTION=1` restores
serial construction in the same binary.

The clean selected process returns the same exhaustive true-negative UNSAT
terminal in 17.131544 wall seconds / 176.171665 total core-seconds /
2,904,489,984 bytes peak RSS. The clean same-binary serial-construction control
takes 20.796178 wall seconds / 172.001422 core-seconds / 2,839,740,416 bytes,
for a 1.214x wall speedup, 1.024x CPU ratio, and 1.023x RSS ratio. Construction
itself falls from 4.629543 to 0.528026 wall seconds, an 8.768x phase speedup.
Three alternating development pairs have 20.694938 versus 17.163728-second
wall medians, a 1.206x speedup.

Every selected and serial arm retains all 512 masks, 270 non-rational skips,
242 constructed and completed systems, 14,278 equations, 14,515,915 terms,
equation fingerprint
`02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb`,
99,199,976,264 word XORs, zero roots, and the same true-negative
classification. A clean 14-thread/batch-64 arm reaches 14.661977 seconds, but
it was selected after seeing this negative target, raises RSS to 3,535,716,352
bytes, and can do extra work on decomposable targets. It is a post-hoc ceiling,
not the selected generic policy.

Selected F4 remains 5.80x slower than direct MITM by wall, 59.75x by CPU, and
64.51x by RSS. Clean build plus the selected run costs 186.059675 wall seconds
/ 340.897134 core-seconds / 2,904,489,984 bytes peak RSS. The measured
native-F4 campaign lower bound through Stage 170 is 8,808.374612 sequential
wall-seconds, 12,423.234481 core-seconds, and 6,310,576,128 bytes maximum
process RSS across 264 components. Complete cost remains `null`.

The optimization was developed after the Stage-169 truth was opened, so it is
same-target engineering rather than another fresh holdout. Licensed Magma, a
larger preregistered native-F4 target distribution, full index-calculus/rho
cost, and unaffiliated reproduction and novelty review remain open. The
canonical additive result is
[`stage-170-parallel-fixed-x1-construction-20260923/result.json`](stage-170-parallel-fixed-x1-construction-20260923/result.json).
