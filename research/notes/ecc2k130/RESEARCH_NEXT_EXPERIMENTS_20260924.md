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
fully charged map construction from supplied torsion generators. None is
an ECC2K-130 challenge-log or whole-attack speed claim. The canonical wall table is in
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
| 1 | PR #735 independently replayed 512/512 n=53 training relations on 220 orbit columns, first full rank at relation 461, 51 post-rank labels, and 8/8 point-only logs; its old producer's SAT accounting defect is isolated. PR #738 independently checked 27 matched signed-Frobenius rho outputs on the same one/eight-point Q prefixes. Charged IC wall was 159.916/170.521 s versus rho ranges 1.081–3.435/4.468–19.256 s. All measured blocks lost, with host-load and post-panel L1 caveats; no general interval or n=131 timing is inferred. PR #739 shows the direct n=131 root-index port has at most 2.12e-13 coverage under an optimistic 1 TiB packed-root cap. | PR #737 completed six clean same-Q n=37/41 full-rank controls: median rho/IC whole-process wall was 0.670 and 0.0634, respectively; n=41 IC peaked at 1,711,898,624 bytes RSS. Freeze at least three nonsaturated base sizes per rung and measure unique-root occupancy, failed-query S3 work, useful rank, cold memory and same-Q batched rho. Test a shared-log batch with one training solve rather than repeating full rank per Q. | Stop the current materialized-root path if its coverage/memory and fully charged cost remain behind rho; reopen only for a reviewed low-memory, low-miss-work search that survives independent rank and point-log replay. Solver-stage savings alone cannot pass. |
| 2 | Earlier K1 comparison measured 2,405,397 versus 53,895,954 elimination word XORs for degree-5 versus degree-8 bases across 12 verified runs each (22.41-fold solver-stage difference, not full cost). PR #716 froze and ran a degree-7 crater-edge proxy over F_(2^21) at 16 useful points per base and 512 paired targets. Original/transported had 131 hits and first full rank at attempt 135; descendant-native/pullback had 117 hits and rank at 202. Cold native field multiplications to verified rank were 360,907 versus 89,948 original (4.012×), including map construction and cofactor costs. All scalars/lifts verified. One toy stream gives no exact-target or rho claim. | Repeat with disjoint target streams and several equal-cardinality/equal-orbit bases before any selector, then test exact degree-263 native bases only after the same group/rank gates and dual-map checks. Charge field, construction, failed targets, solving, lifting, rank and final recovery. | Continue only if held-out total cost to equal verified rank falls; stop a selector if native benefit vanishes when geometry and base construction are charged. |
| 3 | PR #703 established the degree-263 structure and PR #717 the full-point oriented map. PR #743 independently certified the first saved kernel line descends, with two horizontal and 262 descending lines; source discriminant is -7 and leaf discriminant -484183. The least degree of a non-scalar leaf endomorphism is 121046. In the measured Python code, source Frobenius takes two squarings, one map takes 1053 multiplications/400 squarings/one inversion, and generic leaf scalar action takes 25476 multiplications/25604 squarings/193 inversions. These are stage costs, not an attack speed ratio. PR #750 independently verifies 64 forward/reverse full coordinates on all eight saved maps, 4,192 twist-kernel exceptions, and the sign-corrected dual identity; per-representative four-map setup from saved generators costs 82,638 multiplications/80,750 squarings/605 inversions, and a positive dual composition costs 2,106 multiplications/800 squarings/two inversions. | Compare original, descendant-native and pullback bases on blind paired targets at matched useful size and orbit composition. Charge the actual leaf orbit action, or use an explicitly charged tagged-conjugate alternative. Use PR #716's negative degree-7 result as a control. | Reject geometric selection if verified relation yield and total cost do not improve after map construction, orbit action, cofactor and matrix-rank costs. Do not count Frobenius-equivalent kernels as independent trials. |
| 4 | The new retrospective toy m=2 feature panel reran 716 cases from main: 640 natural and 76 planted. Affine row-span contradictions safely rejected 88/640 natural cases and retained all 111 hits; the double holdout by base and Frobenius target orbit rejected 4/55 and retained 13/13 hits. Charged confirmation field-multiplication ratio was 0.9415 and row-XOR ratio 0.9806. The rank priority t in W was inconsistent across fields and is rejected. This is a stage diagnostic, not an ECDLP gain. | Test the certified affine early exit on naturally sampled implicit-orbit and chained m=3 PDPs with cached per-base span, then measure first verified relation, failed work and independent row gain against the same producer without the screen. Only if complementary solvers show held-out headroom, test a charged portfolio. | Retain the safe screen only if matched total cost falls without losing useful relations. Do not train a rank selector on these toy cells; a selector needs a fresh orbit-held-out hindsight advantage over the best constant policy after feature cost and producer misses. |
| 5 | Inherited F4 splitting has a strong stage-engineering signal, but native F4 terminal cases remain slower than matched direct MITM. F5 toy tests in PR #701 find no consistent certified-neighbor advantage. | Test adjacent-assignment and target-batch row-space reuse with exact inherited-ideal controls; benchmark first verified relation and exhaustive UNSAT separately. | Continue only if paired CPU/operation cost beats direct MITM on held-out systems, not merely a prior F4 version. |

## Immediate execution order

1. PRs #735, #738 and #739 resolved the compact n=53 completion, matched
   fixed-stream rho, and direct n=131 architecture transfer gates. The
   current materialized-root four-sum implementation is a no-go as the next
   n=131 port. Preserve the old `claim_relation_yield.json` as historical
   pending-validation input; the accepted rank/point receipts are separate.
2. PR #737 merged six clean same-Q n=37/n=41 full-rank pairs with durable
   raw evidence and independent replay; the point-defined method lost on the
   panel. Run the preregistered three-base-size sweep per rung at nonsaturated
   hit rates, charging failed queries, root occupancy, RSS, rank and matched
   rho. A true shared-log batch needs one base-log solve plus new Q relations,
   not one full-rank solve per Q.
3. PRs #717, #743 and #750 establish full-point degree-263 transport,
   saved kernel directions, leaf order and a sign-corrected exact-target dual.
   The first line descends in both seeds; cheap coordinate squaring does not
   remain a self-map of that leaf. Run paired leaf-native versus transported
   bases with orbit-action, torsion discovery, both maps and setup charged
   before counting a native pullback as an exact-target attack primitive.
4. PR #716's degree-7 equal-useful-cardinality pilot gave no native-base
   advantage on its one toy stream. Repeat several base constructions and
   disjoint target streams matched by useful cardinality and orbit
   composition; test the legally placed PR #706 affine screen on held-out
   chained PDPs. Charge screen setup, all failures, map/cofactor work,
   independent row gain and final scalar recovery.

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
  shared-log batching remains unmeasured. Raw archive SHA-256 is
  `5e01c894c05919a573f2716bbdeace53952cdd009a85969af56aaad591c2c9ed`.
- Degree-263 endomorphism ring and orbit action: merged PR #743 (commit
  `a123997a3da30999002d097858ccb4235a1d72d8`) independently checked
  both saved seeds and all 264 lines per seed. The first map descends; the
  leaf has conductor 263 and discriminant -484183. Its least non-scalar
  endomorphism degree is 121046; the measured generic subgroup scalar action
  is a stage-cost control, not a lower bound on future implementations.
- The normal-x degree-53 local artifacts still need main-based review before
  promotion from provisional evidence.
