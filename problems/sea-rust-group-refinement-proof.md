# SEA Rust-to-group refinement proof is not yet established

- **Status:** open
- **Kind:** proof gap; implementation-verification gap
- **Affected work:** [SEA verification map](../research/sea_point_count_20261016/VERIFICATION.md), PR #1923
- **Next action owner:** unassigned follow-on formal-verification work
- **Last updated:** 2026-10-16

## Problem

The checked abstract BSGS and refinement proofs do not yet establish that the
actual Rust bigint, field and curve operations satisfy their group assumptions
and that every SEA input congruence retains the true trace.

This gap applies to the standalone `isogeny_algos` implementation at revision
44f8d17f1c5b8988a57a16e9dddc8e5ce099d325, for prime fields of characteristic
above three and the supported multiword representations. The formal work
requested for PR #1923 is partially delivered by model checking and abstract
proofs. A claim of complete implementation verification remains open.

## Evidence and uncertainty

- **Observed:** positive TLC configurations, failing legacy/sample-only
  controls, and successful Lean builds with inspected axioms are retained in
  [formal logs](../research/sea_point_count_20261016/formal/logs/) and
  [final Lean logs](../research/sea_point_count_20261016/formal/logs/final/lean-build.txt).
- **Derived or proved:** `SeaProofs.bsgs_complete` holds for arbitrary counts
  and positive widths over a lawful additive group. Refinement preservation,
  singleton soundness, and interval uniqueness have checked proofs with their
  explicit premises. The code/model correspondence is documented but not
  proved as a machine-checked refinement.
- **Suspected explanation:** none required; the missing proof connection is
  explicit, rather than an unexplained failing arithmetic result.
- **Missing capability:** formal definitions and semantics connecting the
  Rust representations and exceptional-case formulas to the group, and
  proofs for SEA congruences, point-order divisibility and CM construction.
- **What this does not establish:** universal correctness of compiled Rust,
  certified analytic Hilbert polynomials, deterministic primality, or timing
  resistance. Exhaustive tests at small fields do not close the general gap.

## What would unblock it

- [ ] Define the concrete prime-field representation, non-singular curve
  points and scalar action; prove `padd` and `pmul_big` refine that action,
  including infinity, doubling, inverse points, normalization and word bounds.
- [ ] Refine the actual baby table and giant-step loop to the checked `Hits`
  predicate and connect returned indices to SEA's signed trace progression.
- [ ] Prove trace parity, Elkies/Atkin constraints and point-order divisor
  accumulation under explicit primality/curve hypotheses; connect them to
  preservation and certificate-uniqueness theorems.

These are sufficient for a stronger SEA search/count correctness claim under
the specified arithmetic and compiler assumptions. CM's analytic polynomial
construction and side-channel properties would still require separate proofs.

## Next discriminating step

- **Action:** formalize one existing small prime-field backend and every
  branch of `padd`, then prove its scalar action obeys `GroupOps`.
- **Inputs and controls:** the exact source revision above, the existing
  infinity/low-order exhaustive tests, and the retained legacy counterexample.
- **Execution limits and stop condition:** no new large search or paid
  infrastructure; finish and kernel-check the representation bridge or report
  the exact unproved branch/lemma, without admitting it as an axiom.
- **Outcomes:** a checked bridge allows a Rust-refinement proof; a genuine
  counterexample becomes a bug and regression; a missing lemma remains open.
- **Recheck trigger:** a checked arithmetic bridge or a change to the relevant
  Rust operations, congruence algorithms or certification assumptions.

## Resolution and history

2026-10-16: abstract proofs and finite-state checks added. The complete
implementation-refinement obligation remains unresolved.
