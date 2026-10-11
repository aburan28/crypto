# Goal: from a fitted ladder to a predictive law

**Set 2026-10-05.** The `m = 3` refutation-degree ladder ended at `ℓ = 6` because the next
rung needs a 90–300 GB matrix ([ic_dense_ladder_20261004](../ic_dense_ladder_20261004/RESULTS.md)).
Its readings so far are a fit: "`x4` refutes at `ℓ + 5`". A fit says nothing about cells
it has not seen, and nothing about the systems that matter for an attack, which cannot
be measured.

**The goal is to replace the fit with a law that predicts.** The candidate is the
generic one: if Weil-descended summation systems behave like semi-regular Boolean
systems, their refutation degree is fixed by a Hilbert series computed from the number
of unknowns and the equation degrees alone, with no measurement. On all twelve measured
cells the reading sits at `d_reg` or `d_reg + 1` of that series. That agreement was
found on the same data, so it is a hypothesis, not evidence.

**What would count as reaching it.**

1. The law is tested where it could fail. Out-of-sample cells are registered before
   measurement, chosen where the generic law and the fixed-`ℓ` fit disagree (the generic
   law says the degree falls as `n` grows; the fit says it cannot).
2. A registered verdict: the systems are generic (the law holds out of sample), or they
   are not (structure exists, in a stated direction).
3. If the law holds, the solving degree of these systems becomes computable for any
   `(n, ℓ)` with an explicit, tested heuristic. That is the input the larger-`m` exponent
   audits lack, since those systems cannot be measured at all.

**What it is not.** It is not a claim about end-to-end attack cost. At `m = 3` and
`n = 131`, index calculus is above rho by `2^70` or more however cheap each solve is
([cost_n131.py](cost_n131.py)); the law cannot change that, and the round does not claim
otherwise.

Rounds under this goal: [PREREGISTRATION.md](PREREGISTRATION.md) (round 1).
