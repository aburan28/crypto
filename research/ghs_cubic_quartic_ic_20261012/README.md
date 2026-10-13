# Plan: index calculus on the `F_q` plane quartic against rho (not executed)

2026-10-12. Conductor T-2280. **This directory is a frozen protocol only.
No solver is implemented, no run has happened, and no result exists.**

## Status

- [PROTOCOL.md](PROTOCOL.md) freezes the question, the two index-calculus
  variants (`ic.full`, `ic.2lp`), the rho reference, the cost unit `S`,
  the charges, the boundary, the predictions, and the decision rule,
  before any measurement.
- **Not done:** the linear-equivalence / divisor-class arithmetic on a
  smooth plane quartic that the solver needs, the solver itself, and every
  measurement. These remain the open obligation recorded in
  [`problems/plane-quartic-jacobian-class-arithmetic-missing.md`](../../problems/plane-quartic-jacobian-class-arithmetic-missing.md)
  (opened with PR #1761).

## What this depends on

The instances are the smooth `F_q` plane quartics built and checked in
`research/ghs_cubic_genus3_build_20261012/` (PR #1761): for each,
`#Jac(C)(F_q) = #E0(F_{q³})` and every transfer divisor tested is
`F_q`-rational. That study established the transfer as an `F_q`-rational
map on points; it did **not** execute it as a group homomorphism, because
no divisor-class arithmetic exists. This plan is what would close that
gap and turn the derived cost comparison into a measured one.

## Why this is a plan and not a run

Implementing divisor-class arithmetic on a non-hyperelliptic genus-3
curve, plus the relation collection and linear algebra of an index
calculus, is a substantial piece of work. The protocol is committed first,
per the repository's rule that a preregistered experiment belongs in a PR
before an outcome exists (AGENTS.md, "Research work belongs in pull
requests"). Execution, if it happens, lands in a follow-on PR that links
this one and the problem record; merging this plan does not close the
workstream.

## Acceptance conditions for the follow-on (from the problem record)

1. A native routine deciding linear equivalence of two `F_q`-rational
   degree-0 divisors on a smooth plane quartic, validated on the nine
   quartics of PR #1761 by `#Jac · D ~ 0` for random `D` and by
   associativity.
2. With it: `φ(P+Q) ~ φ(P) + φ(Q)` and `ℓ·φ(P) ~ 0`, `φ(P) ≁ 0` for a
   generator `P` of `E0(K)[ℓ]`, on at least one case per field.
3. The index calculus of PROTOCOL.md run and counted against rho in the
   unit `S`, every returned answer verified against the planted scalar.

Conditions 1 and 2 establish an executed transfer. Condition 3 is the
cost comparison; only a mean `S` below rho's at the three largest fields,
with every answer correct, earns "below rho on these instances". No
crossover is predicted in advance for either variant.
