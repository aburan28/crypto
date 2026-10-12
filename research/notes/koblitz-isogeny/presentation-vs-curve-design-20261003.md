# Design: is any index-calculus cost a property of the curve, not the setup?

Status: **design and pilot. No measurement run yet.** Pre-registered here.
The pilot only sets budgets and feasibility (§6). It produces no
hypothesis test.

## Why the earlier results are not enough

- **The Gröbner step was trivial.** The n = 31 sweeps (PRs #1209, #1258)
  used m = 2 with l = 6 or 8. Every decomposition there is ruled out after
  1.000 reductions at first fall degree 2, so "algebraically FLAT" carries
  no information about the Gröbner step.
- **The yield spread was a setup effect.** PR #1279 showed it is explained
  with no fitted constants by how the fixed subspace V = ⟨1, z, …⟩ meets
  each curve's cofactor cosets. That is a property of the setup, not of
  the curve.
- **Neither result answers the original question.** That question is
  whether the transported DLP is harder on some members of the isogeny
  class than on others, at parameters where index calculus actually does
  work.

## Instrument changes (implemented in this PR, default off)

1. **Full-group probes** (`IcCostOptions::full_group_probes`). Targets are
   R = aG + bQ + T, with T = [r]X uniform cofactor torsion. Relations are
   recorded [c]-projected, and [c]T = O, so they stay sound. Every root
   pair of S_{m+1} becomes a relation, which removes the cofactor-balance
   effect from yield by construction. The exact-probe mode stays available
   as the positive control: the effect from PR #1279 must reappear there
   and vanish here.
2. **Random factor-base subspaces** (`IcCostOptions::basis_seed`). These
   are l linearly independent random elements of F_{2^n}, so the setup can
   be varied independently of the curve.

Unit test `full_group_probes_over_a_random_subspace_still_solve_soundly`
checks that relations stay consistent and the secret is still recovered.

## Hypotheses

For each metric Y and each setting (n, m, l), fit the random-effects model
Y_{E,V} = μ + u_E + v_{E,V} + ε. Here u_E is the curve effect and v_{E,V}
the setup (subspace) effect. The quantity reported is the curve
intraclass correlation, ICC_E = σ²_u / (σ²_u + σ²_v + σ²_ε), with a 95 %
interval from a parametric bootstrap.

Metrics:
- relations per probe (full-group mode), after dividing by the measured
  |F|^m;
- µs and reductions per decomposition call;
- the first-fall-degree distribution;
- total cost to a verified solve (probes × per-call cost, plus linear
  algebra), in the units of the `ecbench` cost model. Pollard rho at the
  same r, which is class-invariant, is the reference.

| outcome | condition |
|---|---|
| **H0, folklore: cost is a property of the setup, not the curve** | ICC_E 95 % upper bound < 0.05 on every metric |
| **H1, curve structure** | ICC_E ≥ 0.10 with a 95 % lower bound > 0 on any cost metric, confirmed on the hold-out half |
| inconclusive | anything else; reported as such, never rounded to H0 |

If H1 holds, the follow-ups are regressions of u_E on volcano level,
conductor of End(E) and a₆-dependent terms of S_{m+1}. Those regressions
are exploratory and labelled as such.

## Sample

- **Full classes** at n = 17, 19 and 23. The walk already names all their
  members.
- **n = 16**, a₂ ∈ {0, 1}. These volcanoes have depth (33 rank-2 vertices
  at ℓ = 3, 5 at ℓ = 31), which n = 31 lacks entirely, since its class is
  one crater vertex over a flat floor.
- **n = 31:** the crater plus 63 floor members ranked by the existing
  seeded hash, split 32/32 into exploration and hold-out halves.
- **4 random subspaces per curve**, using the same seeds on every curve so
  that V is paired across curves.
- **Negative control:** the twist family at each n. It is not isogenous to
  the class, so it may differ only through its own r and cofactor.

## Settings

These are pinned by the pilot. The rule fixed in advance: the largest l
with m = 2 whose projected cost is at most 2 core-hours per (curve, V) cell
at the probe count needed for about 100 relations. m ≥ 3 is included only
if the pilot shows it fits that rule.

## 6. Pilot (run in this PR)

The pilot uses the Koblitz curve, a single subspace and a handful of
probes per setting. It records µs/call, reductions/call, first fall
degree, measured relations against the (2|F|)^m / (m!·#E) heuristic, and
projected core-hours per cell. Results are in
`experiments/koblitz_ic_pilot.jsonl` and are summarised in
[`presentation-vs-curve-pilot-20261003.md`](presentation-vs-curve-pilot-20261003.md).
