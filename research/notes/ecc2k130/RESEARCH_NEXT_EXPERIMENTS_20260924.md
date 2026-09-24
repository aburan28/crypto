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
| 1 | Compact orbit extraction may buy feasible base size, but the degree-53 local prototype used about 601 MB and about 625 ms median extraction versus a historical 7.49 ms warm explicit-table lookup. Its positive witnesses have an independent 12/12 replay; the source and replay are not yet merged on main. | The producer file is absent from main and currently lives on a dirty branch descended from PR #355, which targeted another feature branch. Record and resolve that upstream dependency first; do not import a 6,000-line example as a small fix. Then land its independent replay, freeze deterministic traversal and exact small-field support/cancellation tests. On one fixed base and target stream compare explicit table, rebuilt compact index, and reused compact index at cold batch sizes 1/8/32/full rank. | Advance only if all lifts and rank/scalar controls pass and a memory-constrained useful base or a measured cold batch break-even emerges. If reuse cannot amortize construction and the larger base gives no useful-rank gain, stop this branch. |
| 2 | Selecting a base by relation trials ranked candidates badly: a tracked K1 comparison used 2,405,397 versus 53,895,954 elimination word XORs for degree-5 versus degree-8 bases across 12 verified runs each, a 22.41-fold solver-stage difference. This is not a full-cost ratio. | Pre-register a cost-per-new-rank base selector. Compare equal useful point cardinality, orbit count, and subgroup handling for ordinary orbit unions and any degree-263 native or pullback bases. Charge construction, unsuccessful targets, solver, lifting, rank and final recovery. | Continue only if held-out total cost to equal rank falls with valid answers; a trial-only or solver-only gain is diagnostic. |
| 3 | PR #703 merged the exact-target degree-263 structural preflight and sampled oriented replay: 264 kernel lines, two horizontal loops, and 262 descending j-invariants. It did not measure PDP yield or implement exceptional-input transport. | Complete a generic normalized full-point Vélu map, including infinity, kernel and exceptional inputs; verify homomorphism, codomain, dual/composition and transported target order. Then blind-compare original, descendant-native and pullback bases at matched useful size. | Reject geometric selection if verified relation yield and total cost do not improve on held-out targets after map construction and cofactor costs. Do not extrapolate from the 262 descendant labels as independent trials. |
| 4 | The new retrospective toy m=2 feature panel reran 716 cases from main: 640 natural and 76 planted. Affine row-span contradictions safely rejected 88/640 natural cases and retained all 111 hits; the double holdout by base and Frobenius target orbit rejected 4/55 and retained 13/13 hits. Charged confirmation field-multiplication ratio was 0.9415 and row-XOR ratio 0.9806. The rank priority t in W was inconsistent across fields and is rejected. This is a stage diagnostic, not an ECDLP gain. | Test the certified affine early exit on naturally sampled implicit-orbit and chained m=3 PDPs with cached per-base span, then measure first verified relation, failed work and independent row gain against the same producer without the screen. Only if complementary solvers show held-out headroom, test a charged portfolio. | Retain the safe screen only if matched total cost falls without losing useful relations. Do not train a rank selector on these toy cells; a selector needs a fresh orbit-held-out hindsight advantage over the best constant policy after feature cost and producer misses. |
| 5 | Inherited F4 splitting has a strong stage-engineering signal, but native F4 terminal cases remain slower than matched direct MITM. F5 toy tests in PR #701 find no consistent certified-neighbor advantage. | Test adjacent-assignment and target-batch row-space reuse with exact inherited-ideal controls; benchmark first verified relation and exhaustive UNSAT separately. | Continue only if paired CPU/operation cost beats direct MITM on held-out systems, not merely a prior F4 version. |

## Immediate execution order

1. Record the compact producer's upstream dependency: current main lacks
   examples/koblitz_s5_sat_instance.rs, and the local implementation is dirty
   on a feature branch descended from PR #355. Land that precursor through a
   separate reviewed PR before changing its traversal or benchmarking it as
   a reproducible main-based candidate. Keep its existing 12/12 replay
   provisional until source and receipts share a reviewable lineage.
2. PR #703 merged the degree-263 scope correction and replay with exact
   source/seed hashes. Build the exceptional-case map before any claimed
   descendant-base advantage.
3. Run a paired, first-verified-relation pilot for equal-cardinality bases.
   Keep n=7/9/13/17 morphology separate from subgroup difficulty: the
   degree-17 sample has prime order 239 and cofactor 548, whereas degree-13
   has prime order 2003 and cofactor 4. Use explicit cofactor clearing and
   compare equal verified workloads.
4. Only after the pilot passes correctness and shows headroom, run the
   pre-registered cold 1/8/32/full-rank cost panel and matched rho. The
   scoreboard and previous baseline move in the same PR as the measurement.

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
- Toy PDP affine-screen grouped holdout: PR #706 contains the reproducible
  source, raw cases, charged native counters and decision.
- The compact orbit and normal-x degree-53 local artifacts still need
  main-based review before they are promoted from provisional evidence.
