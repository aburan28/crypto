# ECC2K-130 index-calculus research program: next measured decisions

Status: research plan and hypothesis ledger, 2026-09-24. This note sets
experiments, not a claim to solve the Certicom challenge. Any number from a
local unmerged worktree is provisional until its source and evidence have a
reviewable PR. Update this ledger in the PR that resolves each gate.

Accepted update, 2026-09-25: PR #735 establishes n=53 full-rank point-only
recovery; PR #738 finds no cold fixed-stream crossover against matched
signed-Frobenius batched rho at one or eight public Q points; PR #739 gives a
conditional counting/memory no-go for directly porting the current
materialized-root, exhaustive-miss design to n=131. PR #737 adds six clean
same-Q n=37/41 full-rank controls, and PR #743 certifies that the first
saved degree-263 kernel descends and quantifies the lost cheap geometric
Frobenius action on its leaf. PR #750 verifies the normalized degree-263
dual on all eight saved maps, including a preserved checker failure and
fully charged map construction from supplied torsion generators. PR #753
adds independently replayed disjoint-orbit toy factor-base cells and a bounded
actual-challenge-point degree-263 smoke with no n=131 m=2 relation or rank.
PR #747 measures sparse n37/41 compact bases, and PR #757 independently
certifies complete group-law four-sum support on held-out toy target
streams: 164/512, 10/512 and 3/128 oracle members exactly match extractor
hits. These finite samples favor a support/representation gate over root-index
repair. PR #762 then passes an exact n=13 rotated normal-basis support gate:
rotated m=5 bases cover 1,591/2,003 projected targets versus 61 for
same-size repeated bases; rotated m=6 covers 2,003/2,003 for both frozen
normal generators. Its n=131 SHA-mask lift-density sample admits m=6,d=21
and m=5,d=25 through a necessary count screen, without measuring target
support or an S6 solver. PR #763 audits the required multiword/streaming
solver corpus and keeps solver timing a separate gate. PR #767 then freezes
four planted and four complete-oracle-certified negative projected PDP targets
at n=13,m=5 and n=19,m=6, with no solver timing. PR #766 independently
certifies the Frobenius-compressed relation-row identity on all 17,917 saved
positive n=13 rotated full-point witnesses, explicit torsion/repeated-column
controls and eight SHA-fixed n=131 public synthetic planted tuples. Its 2/3
canonical columns per toy arm are a structural fact, not rank. PR #768
then excludes all rows for 16 SHA-fixed positive Q holdouts per arm and finds
rank 2/2, 2/2, 2/2 and 3/3, with only 2/2/2/3 independent rows; an explicitly
charged saved-row oracle lets 64/64 point-only scalars replay by group law.
None is an ECC2K-130 challenge-log or whole-attack speed claim. The canonical wall table is in
`docs/index-calculus-scoreboard.html`.

## Decision rule and common controls

The target is a verified discrete logarithm for the public ECC2K-130
instance using a complete, cofactor-correct index-calculus pipeline. Until
then, use exact small-field Koblitz instances as controlled proxies. Record
cold total operations from base construction through scalar verification,
plus CPU and peak RSS as secondary measures. Compare against the same
target workload using automorphism-aware batched rho; independent rho per
target is an invalid batch reference. Record single-target and batch sizes
1, 8, 32, and full-rank separately. Retain every failed target, timeout,
OOM, dependent row, and verification failure.

The primary response variable is cost per new independent, verified relation
and then total cost to a verified scalar. A target relation is counted only
after its group equation is checked. A matrix row is useful only if it raises
rank. A completed run checks [k]P = Q. All claims use the boundary, table,
ratio and frozen-regression rules in AGENTS.md. No stage-only result is
called a DLP speedup.

Freeze source commit, field polynomial, curve and subgroup, cofactor,
factor-base constructor and actual useful point set, target seeds, target
order, solver budgets, machine, and resource limits before comparing
variants. Use the same natural targets for every policy and separately
marked planted controls. Hold out whole target orbits and at least one
field size; do not select a winning variant on its evaluation targets.
Record a machine-readable manifest and exact rerun command with each PR.

## Ranked experiments and stop conditions

| Priority | Decision and evidence | Next controlled experiment | Gate |
| --- | --- | --- | --- |
| 1 | PR #735 independently replayed 512/512 n=53 training relations on 220 orbit columns, first full rank at relation 461, 51 post-rank labels, and 8/8 point-only logs; its old producer's SAT accounting defect is isolated. PR #738 independently checked 27 matched signed-Frobenius rho outputs on the same one/eight-point Q prefixes. Charged IC wall was 159.916/170.521 s versus rho ranges 1.081–3.435/4.468–19.256 s. All measured blocks lost, with host-load and post-panel L1 caveats; no general interval or n=131 timing is inferred. PR #739 shows the direct n=131 root-index port has at most 2.12e-13 coverage under an optimistic 1 TiB packed-root cap; PR #747 tightens the exact unordered four-sum ceiling to 8.82e-15 at that cap. PR #757 finds no oracle-positive extractor gap in 1,152 base-target evaluations over 1,024 fresh target points (R12 reuses an R8 prefix): the exact oracle and extractor agree at 164/512, 10/512 and 3/128. PR #762's complete n=13 cofactor-correct gate finds projected rotated/repeated support of 1,591/61 at m=5 and 2,003/85 or 2,003/369 at m=6 for two frozen beta choices; 81,920 n=131 lift samples pass independent replay and keep the m=6,d=21 count screen open. | PR #737 completed six clean same-Q n=37/41 full-rank controls: median rho/IC whole-process wall was 0.670 and 0.0634, respectively; n=41 IC peaked at 1,711,898,624 bytes RSS. PR #747 completed the nonsaturated base sweep; PR #757 completed its exact fresh-stream four-sum membership gate. PR #762 passes the toy rotated-support screen, while its n=131 density supports only a conditional tuple-count admission. PR #766 passes the exact compressed-row semantic gate on 17,917 archived n=13 witnesses and eight synthetic n=131 tuples, without measuring rank or new relations. PR #767 admits fixed positive/negative projected PDP corpora at n=13,m=5 and n=19,m=6; no solver has run. PR #763 finds a 64-bit Boolean representation limit before an S6/S7 portfolio can run. PR #768 computed exact held-out modular-q toy rank: only 2/2/2/3 independent rows span the four 2/2/2/3-column arms, despite 17,917 valid identities. Its 64/64 point-only logs use a charged saved-row archive oracle, not a new PDP solver. Next preregister n=131 compressed-base cardinality/rank cost alongside the canonical multiword/streaming PDP corpus, complete positive and negative semantics, and an implicit memory/work gate; do not spend on compact root-index repair solely from these matched finite streams. PR #754 completed three L32 and three L128 true shared-log n=53 point-only batches after one rank-220 solve: all 480 IC/rho scalar checks passed; IC child-only lower/rho was 2.688–3.034 and 1.541–1.730 per block. The L128 training-once three-block portfolio bounds straddle rho (0.757 lower, 1.305 verified upper). PR #760 then ran the exact combined L384 union with one IC index and one rho table: 384/384 same-Q scalar replays passed; operational IC/rho wall was 353.649/154.799 s (2.285x), and even the incomplete IC two-child lower was 1.787x rho. Gate further batch-size work on a pre-frozen full-cost model and common-unit calibration; prioritize n=131 low-memory supported-target admission. | Stop the current materialized-root path if its coverage/memory and fully charged cost remain behind rho; reopen only for a reviewed low-memory, low-miss-work search that survives independent rank and point-log replay. Solver-stage savings alone cannot pass. |
| 2 | Earlier K1 degree-5/8 comparison was a solver-stage signal only. PR #716's one-stream degree-7 F_(2^21) 16-point toy native base lost to original at first verified rank. Preregistered PR #753 tested two independent equal-useful-size, same-source-orbit bases on disjoint 512-target signed-τ holdouts. All 16 policies reached rank 17 and recovered scalar 339; independent replay checked 8,192 outcomes, 2,137 target draws, four bases, rank and cold ledgers. Native first-rank attempts were 97/121 versus original 159/148 in seed 2511, but 126/289 versus 115/138 in seed 2512. Its cold field-multiplication component was 2.087–3.700x original in all four cells and beat transported in only two. Degree 7 is ramified, not an n=131 yield model. | Stop the 16-point m=2 native selector. Preregister an implicit m≥3 original/leaf producer only after a supported-target-mass and memory admission bound predicts observable independent rank. Use held-out natural targets, cofactor-correct lifts, actual or tagged-conjugate leaf action, all misses, dual/map setup and same-Q rho. | Reopen native selection only if equal-useful-size held-out rank and fully charged cost beat original after action, construction, solver, projection and recovery. The toy field-multiplication component fails; modular-r work has no calibrated conversion and no whole-ECDLP speed ratio is asserted. |
| 3 | PR #703 established the degree-263 structure and PR #717 the full-point oriented map. PR #743 independently certified the first saved kernel line descends, with two horizontal and 262 descending lines; source discriminant is -7 and leaf discriminant -484183. The least degree of a non-scalar leaf endomorphism is 121046. In the measured Python code, source Frobenius takes two squarings, one map takes 1053 multiplications/400 squarings/one inversion, and generic leaf scalar action takes 25476 multiplications/25604 squarings/193 inversions. These are stage costs, not an attack speed ratio. PR #750 independently verifies 64 forward/reverse full coordinates on all eight saved maps, 4,192 twist-kernel exceptions, and the sign-corrected dual identity; per-representative four-map setup from saved generators costs 82,638 multiplications/80,750 squarings/605 inversions, and a positive dual composition costs 2,106 multiplications/800 squarings/two inversions. PR #753's bounded exact smoke checked 32 public Certicom P,Q combinations on each of two nonconjugate descending degree-263 lines with source, native, transported and pullback 16-point bases. It charged dual/map construction, cofactor projection, full map covariance and one native leaf-orbit neighbor per selected point; independent replay checked 64 targets and 256 policy-target cases. All four m=2 arms had 0/32 hits and rank 0 on both lines, as the 32×136/r expectation is at most 6.40e-36. Prior torsion discovery is a separate 3.267-second receipt without compatible field-operation metering. This is representation feasibility, not PDP yield. | Gate larger exact-leaf work on an implicit m≥3 producer whose supported-target-mass and memory bound make useful rank plausible. Charge torsion discovery, both maps, actual or tagged-conjugate leaf action, cofactor, misses, solver, rank, recovery and same-Q rho. | Reject geometric selection if verified relation yield and total cost do not improve after map construction, orbit action, cofactor and matrix-rank costs. Do not count Frobenius-equivalent kernels as independent trials. |
| 4 | PR #706's retrospective m=2 affine S3 contradiction retained every toy hit, but PR #751's preregistered chained m=3 implicit-orbit pilot is negative after full charging. Its 2,048 cofactor-masked full-group targets gave 7/2,048 training-base and 4/2,048 held-base group hits; the double-held-out cell had 2/649 independent hits and verified the toy scalar. The screen safely skipped 18,753 held-out source residuals, yet confirmation field multiplications rose 601,634 to 1,382,384 (2.30x), and transported-base work rose to 10,491,134 (17.44x). Both bases reached verified global rank; every source/transported relation trajectory matched. | Stop placing this affine screen ahead of a complete finite-base pair lookup. Test it again only where the complete alternative is a materially costlier algebraic solver, with disjoint orbits, failed work, feature setup, useful rank and scalar replay charged. Retain the preserved saturation and zero-hit controls as base/target-distribution warnings. | Reopen only on a held-out reduction in complete charged solve cost with zero lost S3 roots/group hits and at least equal verified rank. Skip count or lookup count alone fails. No n=131 or rho inference. |
| 5 | Inherited F4 splitting has a strong stage-engineering signal, but native F4 terminal cases remain slower than matched direct MITM. F5 toy tests in PR #701 find no consistent certified-neighbor advantage. | Test adjacent-assignment and target-batch row-space reuse with exact inherited-ideal controls; benchmark first verified relation and exhaustive UNSAT separately. | Continue only if paired CPU/operation cost beats direct MITM on held-out systems, not merely a prior F4 version. |

## Immediate execution order

1. PRs #735, #738 and #739 resolved the compact n=53 completion, matched
   fixed-stream rho, and direct n=131 architecture transfer gates. The
   current materialized-root four-sum implementation is a no-go as the next
   n=131 port. Preserve the old `claim_relation_yield.json` as historical
   pending-validation input; the accepted rank/point receipts are separate.
2. PR #737 merged six clean same-Q n=37/n=41 full-rank pairs with durable
   raw evidence and independent replay; the point-defined method lost on the
   panel. PR #747 completed the nonsaturated base sweep and PR #757 matched a
   complete group-law four-sum oracle to its compact extractor on fresh
   disjoint targets in all three tested bases. PR #762 then passes the exact
   n=13 rotated m=5/6 support gate and samples n=131 lift density; it does
   not solve an n=131 PDP. PR #763 audits the direct S6/S7 solver interfaces
   and requires a multiword/streaming exporter with group-law validation.
   PR #767 supplies exact positive/negative projected PDP corpora and #766
   certifies signed Frobenius-compressed rows on archived toy witnesses and
   synthetic n=131 tuples. PR #768 resolves the toy rank gate: 2/2/2/3
   independent rows after frozen Q-level exclusions and 64/64 point-only
   logs through a charged archive oracle. It does not measure n=131 solver
   yield. Next bound n=131 compressed-base cardinality and rank cost;
   prioritize the multiword/streaming solver and bounded implicit memory/work
   gate; retain the
   n=131 materialized four-sum no-go. PR #754 completed L32/L128 true shared
   n=53 controls, and PR #760 completed the combined same-Q L384 panel:
   even its incomplete two-child IC lower exceeded rho by 1.787x.
3. PRs #717, #743 and #750 establish full-point degree-263 transport,
   kernel directions, leaf order and an exact positive dual. PR #753's two
   nonconjugate descending lines also pass bounded public-target transport,
   dual pullback and native orbit-neighbor checks, but the 16-point m=2
   pair tables have no exact-target hit or rank. Require an admitted implicit
   m≥3 PDP producer before spending on a larger leaf-native attack run.
4. PR #716's one-stream native-base diagnostic is superseded by #753's
   two-construction, disjoint-orbit holdout. All toy scalars verified but the
   native cold field-multiplication component lost to original in every
   cell. PR #751's held-out chained affine screen also lost after full
   charging. Stop these selectors in their tested placement; prioritize the
   supported-target-mass and memory admission gates.

## Alternative ideas kept testable

A descended order has a different conductor, but isogeny preserves the
rational group order and transports target points. The discriminant
difference alone predicts no easier PDP. Test whether the map produces a
geometric base with larger independent-row yield per construction and
solver cost. At the same time, count Frobenius-equivalent kernels as one
orbit, not independent opportunities. A low first-fall degree is not a
solving-degree bound. At degree 131 the favorable invariant subspaces of
smaller composite degrees are unavailable; validate an actual degree-131
base constructor before projecting any small-field trend.

For solver work, compare concrete workloads rather than names. Gray-code
FES is a control for suitable quadratic systems; it does not cover the
cubic chained systems. SAT, F4/F5, crossbred, and msolve enter a portfolio
only with an identical encoding, correctness gate, cost unit and resource
limit. A solver that proves UNSAT cheaply may be useful even if it rarely
finds relations, but its screening cost and missed producers belong in
the total.

## Required PR record

Each experiment PR includes a protocol before the outcome, code and pinned
dependencies, frozen input/target manifest, baseline and candidate raw
records, independent or differential correctness checks, a single-unit
cost table, negative results and limitations, the decision to continue or
stop, and the updated canonical scoreboard when a comparable end-to-end
measurement changes it. Large data use the AGENTS.md content-hash manifest
rule. The plan itself changes when a gate is resolved; it is not a queue of
unowned speculative tasks.

## Evidence anchors

- Tracked base-selection note:
  research/notes/index-calculus/RESEARCH_FACTOR_BASE_SOLVE_COST.md.
- Tracked inherited-F4 note:
  research/notes/ecc2k130/RESEARCH_INHERITED_F4.md.
- Tracked FES/SAT note:
  research/notes/ecc2k130/RESEARCH_WDSAT_IC_UNIFICATION.md.
- Certified toy F5 neighbor comparison: PR #701.
- Degree-263 structural preflight and representative replay: PR #703.
- Full-point oriented Vélu transport with exceptional inputs and order checks:
  merged PR #717.
- Frozen degree-7 equal-useful-cardinality factor-base pilot: PR #716;
  512 paired targets, all group lifts/scalars verified, native cold cost worse
  on this one toy stream. This is a diagnostic, not exact-target evidence.
- Toy PDP affine-screen grouped holdout: PR #706 contains the reproducible
  source, raw cases, charged native counters and decision.
- Compact orbit producer provenance: PR #714 merged to
  cursor/ic-boundary-experiments-d111; PR #720 promoted its bounded producer,
  frozen n=53 base and replay script to main. The PR validation replayed
  12/12 archived public synthetic witnesses, while the old
  `claim_relation_yield.json` is historically `PENDING_INDEPENDENT_VALIDATION`
  and pins a dirty feature-worktree source. It is not a main-based full-rank
  or matched-rho receipt.
- Main-based compact full-rank and point-only recovery: merged PR #735
  (commit `d92439080c2d5c3a858f71abfc8e749e273125bb`) froze 512/512
  n=53 relations, first rank 220 at 461, 51 post-rank labels, and 8/8
  point-only logs with independently checked [d]G = Q. Training wall was
  158.127 s; the cold eight-point process was 12.394 s. This is public
  synthetic completion, not a rho or ECC2K-130 win.
- Matched n=53 signed-Frobenius rho: merged PR #738 (commit
  `76869fd7e4979f68b601c22df8c61c21884bdd0d`) verified all 27 rho
  output scalars and the IC one-point addendum. On the fixed one/eight-Q
  streams, fully charged compact IC lost every measured cold wall block;
  host load and lack of contemporaneous pairing limit inference.
- Current-design n=131 transfer: merged PR #739 (commit
  `b0428f1c1e5a34e24b9c7a5d2bbfd02556b42db6`) bounds four-sum coverage
  and materialized-root memory with exact integer/rational calculations.
  At 1 TiB optimistic packed roots, p <= 2.12e-13; permitting 1% coverage
  already projects 239 PB of roots. This is a conditional architecture no-go,
  not an n=131 timing or a claim against all index calculus.
- n=37/n=41 same-Q full-rank controls: merged PR #737 (commit
  `ab4a25ea46fcbef8b5f0f3d77cfe66196912df06`) completed six clean
  pairs and independently replayed 820 relations, 38,334 orbit labels and
  all scalar checks. Median rho/IC wall was 0.670 at n=37 and 0.0634 at
  n=41. The point-defined full-rank fixture repeats its solve for every Q;
  shared-log batching is now measured at n=53 in PR #754; per-block
  child-only lower cost exceeds matched rho in all six fixed streams, while
  the training-once L128 three-block portfolio spans rho. Raw archive SHA-256 is
  `5e01c894c05919a573f2716bbdeace53952cdd009a85969af56aaad591c2c9ed`.
- Degree-263 endomorphism ring and orbit action: merged PR #743 (commit
  `a123997a3da30999002d097858ccb4235a1d72d8`) independently checked
  both saved seeds and all 264 lines per seed. The first map descends; the
  leaf has conductor 263 and discriminant -484183. Its least non-scalar
  endomorphism degree is 121046; the measured generic subgroup scalar action
  is a stage-cost control, not a lower bound on future implementations.
- Chained m=3 affine-screen orbit holdout: PR #751 freezes projected, raw, torsion-masked and 2,048-target power protocols before their respective measurements. On the masked full-group stream, source/transported factor bases agree on all group witnesses and rank; held double holdout gives 2/649 independent group hits and exact toy scalar recovery. The certified screen is slower after setup and failed-target work; compact receipts and independent replay are in research/ecc2k130_pdp_chain_holdout_20260925/. No ECC2K-130/rho claim.
- Disjoint-orbit factor-base replication and exact smoke: PR #753 preserves
  preregistration commit `d1e4e56dcd0838f7819a6d38894601cc8a8359a2`, all
  16 degree-7 toy rank/scalar cells, 8,192 independently replayed outcomes,
  two nonconjugate exact degree-263 descending maps, 64 independently replayed
  public Certicom P,Q target coordinates, and 0/32 m=2 hits per policy/line.
  The latter has no rank and no ECC2K-130 DLP or rho implication.

- n=53 true shared-log L32/L128 batch control: PR #754 independently replayed
  512 training relations and 480 IC point logs plus 480 same-Q batched-rho
  scalars across six cold paired blocks. Training first reached rank 220 at
  relation 461. Exact IC source SHA-256 is
  `c2bc8b05087df69bef9593363e9d7c112e843ef16da122da50eb29ab22115f09`;
  accepted 119-file archive SHA-256 is
  `6100f9b47da4a5fb6fa3ee152f4ec7376368fe09749db1aa196e4df204d4dfeb`.
  The first failed checker attempt remains in Git separately. All six
  per-block IC lower/rho ratios exceed one; the L128 three-block portfolio
  interval spans rho, so no general many-target, common-unit S, n=131 or
  ECC2K-130 crossover is claimed. That combined one-table question is now
  measured in PR #760 below.
- Complete four-sum support versus compact S3 extraction: PR #757 froze the
  #747 point sets and disjoint target files before outcomes, archived six
  matched process receipts plus three invalid/contended attempts, and
  independently rebuilt all 24,753 / 215,496 / 484,620 unordered pair sums.
  The complete oracle found 164/512, 10/512 and 3/128 members, exactly the
  unchanged extractor hit counts, with zero oracle-positive misses. All 237
  distinct positive four-multiset witnesses and 16 exhaustive negative
  samples per arm replayed. The SHA256SUMS ledger is
  `eee6dec54fccb590091abffff954d42e349becae94e68288a73729c5548bd840`.
  This is a finite toy Gate 1 result, not global S3 completeness or n=131 yield.
- n=53 combined L384 one-table control: preregistered PR #760 used the exact
  384 distinct Q values from PR #754's three L128 blocks, one fresh 512-target
  rank-220 training solve, one compact point index and one same-Q rho table.
  Independent audit checked 512 training relations, 23,320 orbit labels,
  220 base logs and all 384 IC/rho scalars. Fully charged operational IC was
  353.649 s (358.745 s CPU, 396.6 MiB peak) versus rho 154.799 s
  (154.773 s CPU, 76.8 MiB peak); even IC's incomplete two-child lower
  cost 276.592 s exceeded rho. The 30-file raw archive SHA-256 is
  `bef32e5d3993fe0a6cb5cc6be7c92b8b094cbaa24ff01031400d3f7e3df583b6`.
  A post-outcome archive checker used a 64-hex file-digest predicate on a
  40-hex Git SHA; the failed CI run, exact old checker and correction are
  preserved, and the corrected archive-only CI independently replayed all logs. Total panel elapsed 659.421 s, below the 2,700 s ceiling.
  This is a fixed-stream no-crossover with pre-certified base construction
  uncharged; no common-unit S, general batch limit, n=131 or ECC2K-130
  speedup is inferred.
- The normal-x degree-53 local artifacts still need main-based review before
  promotion from provisional evidence.
