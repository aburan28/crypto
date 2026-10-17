# Cached quintic Frobenius: derivation and discriminating experiment

Status: written derivation, counterexamples, primary-source review and a blocked
execution card. No new solver implementation, execution, timings, machine-checked
proof, external-tool win or 8x claim. The research program remains active.

The proposed characteristic-two operator reuses X^6 mod H and X^8 mod H when
evaluating the domain annihilator in K[X]/H. It retains coefficient squaring:
the map is Frobenius-semilinear over K. A GF4 counterexample rules out treating
it as K-linear. The construction works for repeated H as well as squarefree H.

- [Exact contract and derivation](CONTRACT.md)
- [Derived counterexamples](counterexamples.json)
- [Primary-source and comparator review](PRIMARY_SOURCE_REVIEW.md)
- [Explanatory diagram](flow.svg), [complete PDF source](REPORT.typ)
- [Next bounded experiment](../../backlog/XP-20261017-quintic-support-frobenius/PROTOCOL.md)
- [Open all-roots/8x obligation](../../problems/generalized-rootfinding-all-roots-8x.md)

The dense operation model predicts a modest reduction, not an 8x guarantee.
For 18 recurrence steps plus three nilpotence squares it changes the specified
count from 105S+510M to 105S+315M, including cache setup. Its model ratio is at
most 34/21. This conditional arithmetic model is not measured wall time, a
lower bound on every competing algorithm, or an impossibility result.

Execution limitation observed 2026-10-17: the managed executor reports offline.
The earlier native-source upload was interrupted; its expected remote tree was
not present when checked. No native checkpoint PR or new run is claimed here.
The saved local artifact was previously recorded as
native_quintic_resultant_20261016_v1.zip, SHA256
2f6748a2b20bccfd99efbbd2850d276b09d6c358be019888ce6dd09117b539c7;
it is currently inaccessible and has not been revalidated in this cycle.

Conductor, cargo tests, backlog validation, Lean, identity-certificate issuance
and PDF rendering could not run without the executor. The SVG and complete
Typst source are supplied; there is no rendered PDF. This is the fallback in
[research-visuals](../../.agents/skills/research-visuals/SKILL.md), not a claim
that those checks passed. No software was published to external compute,
no task was posted to cairn and no compute was rented.

Canonical IC/rho timelines, leaderboards and exponent graphs remain unchanged:
this cycle adds no measured cost, verified pipeline ratio, curve identity or
exponent. The new diagram marks the operator as proposed and unbenchmarked.
Merging this derivation/plan would not complete the 8x objective.
