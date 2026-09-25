# ECC2K-130 index-calculus research program: next measured decisions

Status: research plan and hypothesis ledger, 2026-09-24. This note sets
experiments, not a claim to solve the Certicom challenge. Any number from a
local unmerged worktree is provisional until its source and evidence have a
reviewable PR. Update this ledger in the PR that resolves each gate.

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
| 1 | PR #720 promoted the compact-orbit S5 producer, arithmetic dependency, archived 23,320-point n=53 base and replay script to `main` from feature-branch PR #714. Its PR validation replayed 12/12 archived public synthetic witnesses. The old `claim_relation_yield.json` still says `PENDING_INDEPENDENT_VALIDATION` and pins a dirty feature-worktree source; treat it as historical input, not a main-based rank or speed receipt. | On main, freeze an independent full-rank/point-only recovery receipt and compare explicit-table, rebuilt-index and reused-index policies at cold batch sizes 1/8/32/full rank. First correct batch-only SAT/cost labels and validate cancellation cases. Independently replay every relation, row, recovered log and failure; then run matched signed-Frobenius batched rho on the same target workload. | Advance only if complete lifts/rank/scalars pass and a memory-constrained useful base or measured cold batch break-even emerges. A stage relation hit or amortized lookup alone is insufficient. If reuse cannot amortize construction and the larger base gives no useful-rank gain, stop this branch. |
| 2 | Earlier K1 comparison measured 2,405,397 versus 53,895,954 elimination word XORs for degree-5 versus degree-8 bases across 12 verified runs each (22.41-fold solver-stage difference, not full cost). PR #716 froze and ran a degree-7 crater-edge proxy over F_(2^21) at 16 useful points per base and 512 paired targets. Original/transported had 131 hits and first full rank at attempt 135; descendant-native/pullback had 117 hits and rank at 202. Cold native field multiplications to verified rank were 360,907 versus 89,948 original (4.012×), including map construction and cofactor costs. All scalars/lifts verified. One toy stream gives no exact-target or rho claim. | Repeat with disjoint target streams and several equal-cardinality/equal-orbit bases before any selector, then test exact degree-263 native bases only after the same group/rank gates and dual-map checks. Charge field, construction, failed targets, solving, lifting, rank and final recovery. | Continue only if held-out total cost to equal verified rank falls; stop a selector if native benefit vanishes when geometry and base construction are charged. |
| 3 | PR #703 established the exact-target degree-263 orientation structure (two horizontal loops, 262 descendant j-invariants); PR #717 merged a normalized full-point Vélu map with infinity, kernel and exceptional inputs, independent full-point replay and order checks. No descendant PDP yield or DLP advantage has been measured. | Verify the dual/composition on representative exact-target lines, then compare original, descendant-native and pullback bases on blind paired targets at matched useful size and orbit composition. Use PR #716's negative degree-7 result as a control, not an extrapolation. | Reject geometric selection if verified relation yield and total cost do not improve after map construction, cofactor and matrix-rank costs. Do not count Frobenius-equivalent kernels as independent trials. |
| 4 | The new retrospective toy m=2 feature panel reran 716 cases from main: 640 natural and 76 planted. Affine row-span contradictions safely rejected 88/640 natural cases and retained all 111 hits; the double holdout by base and Frobenius target orbit rejected 4/55 and retained 13/13 hits. Charged confirmation field-multiplication ratio was 0.9415 and row-XOR ratio 0.9806. The rank priority t in W was inconsistent across fields and is rejected. This is a stage diagnostic, not an ECDLP gain. | Test the certified affine early exit on naturally sampled implicit-orbit and chained m=3 PDPs with cached per-base span, then measure first verified relation, failed work and independent row gain against the same producer without the screen. Only if complementary solvers show held-out headroom, test a charged portfolio. | Retain the safe screen only if matched total cost falls without losing useful relations. Do not train a rank selector on these toy cells; a selector needs a fresh orbit-held-out hindsight advantage over the best constant policy after feature cost and producer misses. |
| 5 | Inherited F4 splitting has a strong stage-engineering signal, but native F4 terminal cases remain slower than matched direct MITM. F5 toy tests in PR #701 find no consistent certified-neighbor advantage. | Test adjacent-assignment and target-batch row-space reuse with exact inherited-ideal controls; benchmark first verified relation and exhaustive UNSAT separately. | Continue only if paired CPU/operation cost beats direct MITM on held-out systems, not merely a prior F4 version. |

## Immediate execution order

1. PR #720 has promoted the compact producer to main; PR #714 is its
   feature-branch provenance. Preserve the archived claim file's historical
   pending-validation status. The next main-based receipt must freeze
   full-rank and point-only recoveries, correct batch-only SAT/cost accounting,
   replay group lifts and rows independently, and record fresh RSS. PR #735
   is open for this gate as of 2026-09-25; its result is pending until merged.
2. PR #717 merged the oriented full-point degree-263 map, including valid
   exceptional inputs and exact-order checks. Prove and replay a normalized
   degree-263 dual/composition on representative lines before treating
   descendant pullback as an exact-target attack primitive.
3. PR #716's frozen degree-7 equal-useful-cardinality pilot reached verified
   rank 17 but gave no native-base advantage on its one toy stream. Repeat
   several base constructions and disjoint target streams matched in useful
   cardinality and orbit composition. Keep all failed targets, cofactor work,
   map setup and matrix-rank cost; promote a selector only after a held-out
   lower total to the same verified rank.
4. The compact n=37/41/53 cold 1/8/32/full-rank and matched batched-rho
   panel proceeds on its own correctness and full-cost gates; it does not
   depend on a positive degree-7 native-base result. Compare identical target
   workloads and charge setup, failed extraction, lift, rank and recovery.
   Keep S, rho ratio and speedup null until both arms verify every scalar.
   A negative native-base result may instead close that base-selection branch.

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
- Main-based compact rank and point-only recovery: PR #735 is open as of
  2026-09-25, with claimed 512/512 n=53 relations, rank 220 at attempt 461
  and eight point-only logs. Do not promote those numbers to accepted evidence
  until its final head, independent replay and applicable CI pass and the PR
  merges. Matched rho and n=37/n=41 remain open either way.
- The normal-x degree-53 local artifacts still need main-based review before
  promotion from provisional evidence.
