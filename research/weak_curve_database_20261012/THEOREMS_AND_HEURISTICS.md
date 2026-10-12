# Theorem, screen and heuristic register

The [theorem table](theorems.json) has ten scoped entries. Four are
published reductions or bounds (Hasse/order, anomalous trace one,
Pohlig–Hellman, and pairing transfer). The cubic cocycle entry has two
new checked Lean declarations for its algebraic core. The fourth-power
implication for the odd cubic torus has four further checked declarations,
including the exact $q=p^2$ case. Two ISO-1 entries
record the proved cubic conductor-depth restriction and exact two-branch
class-support result with their remaining Lean geometry bridges. The
finite two-branch count identities have checked Lean proofs in the
parent merged package. The dimension guard for a surjection onto a
full degree-n Weil restriction remains a standard consequence whose
abelian-variety formalization is open. The [proof map](PROOF_MAP.md)
names each precise proof boundary.
The [family proof-obligation crosswalk](proof_obligations.json) covers all
21 entries and records the next exact bridge for each.

| Rule | Scope and use | Status |
| --- | --- | --- |
| lambda^(q²+q+1)=1 for lambda=alpha^q/alpha | Cubic cyclic tower; construct a norm-one parameter | Lean algebra checked under explicit automorphism and Frobenius hypotheses |
| lambda=mu^4 when lambda^(q²+q+1)=1 and q=p² for odd p | Cubic torus; exact fourth root mu=lambda^((q²+q+2)/4) | Lean implication checked under the stated exponent equation |
| v₂(f_pi)>=2 | Necessary screen for the p>3 cubic norm-one ordinary class, both twist signs | Proved in parent report; full Lean geometry bridge open |
| Selected CM-order intersection | Exact class decision for the two named hyperelliptic genus-three branches, with exceptional-unit oracle | Parent theorem/census; Lean CM bridge open |
| Z(p)<=24p | Finite forecast for union-zero classes under the exact two-branch input law | Fits seven complete censuses; open |
| CM-mass occupancy | Statistical forecast of zero rates after freezing selected-order masses and sampling law | Open, not evaluated by this import |
| Maximal-order class number alone | Candidate occupancy covariate | Refuted as an exact classifier by p41 trace pair |
| 2p³−abs(t)<=64p for every zero row | Former boundary-only screen | Refuted in 1,184 of 2,074 exact union-zero rows |

The [heuristic table](heuristics.json) records each forecast's domain,
calibration, falsifier and link to the
[hypothesis ledger](../weak_curve_conjectures_20261012/conjectures.json).
The p41 pair is concrete: traces 137406 and 137410 have equal Frobenius
conductor 2-depth 3, equal two-splitting flag 1 and saved
maximal-order class number 414. Their cubic representative counts are
0 and 12. This rules out an exact classifier using only those three
covariates in that diagnostic. Statistical usefulness requires
held-out calibration.

These labels are mechanism-specific. A curve with a small subgroup, a
cover, a low embedding degree or an isogenous covering representative
may satisfy several entries at once. No total family count is obtained
by adding their populations.
