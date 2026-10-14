# Plan: index calculus on the `F_q` plane quartic against `rho.negation` (not executed)

2026-10-12. Conductor T-2280. **This directory is a frozen protocol only.
No solver is implemented, no ecbench arm exists, no run has happened, and
no result exists.**

## Status

- [PROTOCOL.md](PROTOCOL.md) freezes the question, the eleven workload
  identities, the two index-calculus variants (`ic.quartic_full`,
  `ic.quartic_2lp`), the `rho.negation` reference, the cost unit `S`, the
  pinned charges (including hash lookups and rho canonicalization), the
  per-variant counting floors, the predictions, the localizing checks and
  the per-workload decision rule, before any measurement.
- **Not done:** the divisor-class arithmetic on a smooth plane quartic
  that the solver needs, the ecbench curve and method extension, and every
  measurement. These remain the open obligation recorded in
  [`problems/plane-quartic-jacobian-class-arithmetic-missing.md`](../../problems/plane-quartic-jacobian-class-arithmetic-missing.md)
  (opened with PR #1761).

## Scope

The plan is scoped to the **eleven instances PR #1761 built and
validated** (`q ∈ {5, 7, 11}`; PR #1761 head
`adf5d9b3c99443af74954b6fa0e77c735cba2e68`), named by their field
indices, with field basis, `E0` model, generator `G = [4]P*` and
per-target scalars fixed deterministically in the protocol and in a
manifest the follow-on PR commits with hashes. It does not claim
instances above `q = 11`; a larger panel would need its own build and
validation.

The transfer needs `φ` only at auxiliary points whose preimages are
`K`-rational (`(u·u')(t0)` a nonzero square). PR #1761 left the
`F_{q⁶}`-conjugate points untested and no `F_{q⁶}` arithmetic exists, so
the plan restricts the auxiliary points to that eligible set, charges the
search for them, reports the eligible fraction per instance, and makes no
unrestricted-transfer claim.

## How it conforms to the repository's rules

- **ecbench arm, not a standalone program.** The comparison is specified
  as an ecbench curve construction plus two registered IC methods against
  `rho.negation`: one-unit, one-target, sealed, replay-audited
  (AGENTS.md section 12; skills `ecbench-extend`, `ecbench-measure`).
- **Strongest single-target reference.** `rho.negation`, with its
  canonicalization charged; no "without negation" arm.
- **One target per row.** The verdict is per paired one-target workload
  (instance, target, variant); failures stay in the table; any batch mean
  is descriptive only (AGENTS.md, "Primary comparison uses one target").
- **No bound-versus-bound ratios.** Lookups and canonicalization have
  pinned conversion costs, preregistered; a row still `bounded` stays a
  counted-work diagnostic and cannot earn "below rho".
- **Honest boundaries.** A derived counting floor per variant; the
  `(q+1)·c_pt` setup term is relabelled a setup-phase prediction.
- **Correctness.** A wrong scalar invalidates the composite pipeline only;
  separately preregistered homomorphism, order and relation checks
  localize a transfer failure.
- **Rust only.** No Python anywhere in execution, checking or result
  generation; native SVG from recorded values and a Typst PDF.

## Why this is a plan and not a run

Implementing divisor-class arithmetic on a non-hyperelliptic genus-3
curve, the ecbench extension, relation collection and linear algebra is
substantial. The protocol is committed first, per the repository's rule
that a preregistered experiment belongs in a PR before an outcome exists
(AGENTS.md, "Research work belongs in pull requests"). Execution, if it
happens, lands in a follow-on PR that links this one and the problem
record; merging this plan does not close the workstream.

## Acceptance conditions for the follow-on

1. A native routine deciding linear equivalence of two `F_q`-rational
   degree-0 divisors on a smooth plane quartic, passing check C1 on the
   eleven quartics of PR #1761.
2. Checks C2 and C3: `φ(P+Q) ~ φ(P) + φ(Q)`, `ℓ·φ(P) ~ 0`,
   `φ(P) ≁ 0` for a generator of `E0(K)[ℓ]`, at least one case per
   instance with eligible `ℓ`-subgroup points; C4 on emitted relations.
3. The ecbench extension and spec of PROTOCOL.md run as a session with
   `ecbench verify --replay`, one row per paired one-target workload,
   every returned answer verified externally by the runner, with the
   committed session, report, graph and PDF the `research-visuals` skill
   requires.

Conditions 1 and 2 establish an executed transfer. Condition 3 is the
cost comparison, read per workload; no crossover is predicted in advance
for either variant and none is claimable from three values of `q`.
