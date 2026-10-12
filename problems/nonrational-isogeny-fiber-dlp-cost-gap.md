# Complete DLP cost after non-rational-kernel isogeny transfer

- **Status:** open
- **Kind:** implementation and measurement gap; no mathematical
  nonexistence claim
- **Affected work:** [fiber report](../RESEARCH_ISOGENY_FIBERS_ECDLP.md),
  [2026-10-12 correction and bounded control](../research/isogeny_torsion_fibers_20261012/REPORT.md)
- **Next action owner:** unassigned
- **Last updated:** 2026-10-12

## Problem

There is no complete, matched, same-target cost comparison showing whether
transport through an explicit `F_q`-isogeny with a non-rational geometric
kernel helps a coordinate-based ECDLP solver. On the checked
`E/F₁₃: y²=x³+2 → E'/F₁₃: y²=x³+11` pair, the rational kernel is trivial
and the degree-3 map is a bijection on rational points. This validates a
structural distinction, not a speed result. The target subgroup has order
19; extrapolation from this finite control is not an acceptance criterion.

## Evidence and uncertainty

- **Observed:** [census.tsv](../research/isogeny_torsion_fibers_20261012/census.tsv)
  records 349 prime-order rows in a frozen scan of 1,014 curves; the
  [map check](../research/isogeny_torsion_fibers_20261012/map_check.txt)
  passed all 361 rational ordered-pair checks on the named `F₁₃` pair.
- **Proved:** the corrected cofactor statement depends on the order of
  the rational kernel, not the geometric degree. Its arithmetic step is
  compiled in [FiberCofactor.lean](../research/isogeny_torsion_fibers_20261012/lean/FiberCofactor.lean).
- **Suspected explanation:** any useful advantage for this type of kernel
  must arise in the destination's coordinate-based solver or its
  implementation; no such advantage is measured here.
- **Missing capability:** an audited native complete-solve comparison on a
  relevant prime-field pair with all setup, map, relation, algebra and
  recovery costs, plus matched rho and independent replay.
- **What this evidence does not establish:** any speed ratio, asymptotic
  improvement, or result for deployed curves or other fields.

## What would unblock it

- [ ] Freeze a public same-target source/codomain pair with explicit map,
  model/field/subgroup identities, a native complete-solve path, equal
  resources and stop caps, and a matched rho arm.
- [ ] Obtain verified logs for both solver arms, full phase accounting,
  failed-attempt records, and an independent replay; report the cost ratio
  and uncertainty in the repository's counted unit.

Meeting these conditions would settle only the measured pair and workload.
It would not prove a general exponent change or a result for another curve.

## Next discriminating step

- **Action:** design the above frozen paired protocol and verify the native
  destination solver can handle both isogenous models before launching the
  measured session.
- **Inputs and controls:** public prime-order challenge, explicit isogeny
  and dual or return map, same target transported in both arms, matched rho.
- **Limits:** use the repository's ecbench isolation and receipt gates;
  stop on a correctness or accounting failure and retain it.
- **Outcomes:** a verified complete-cost ratio below rho at growing sizes
  would motivate expansion; a stage-only improvement remains a stage
  diagnostic; a timeout, invalid map or missing solver support leaves the
  performance question open and identifies the exact cause.
- **Recheck trigger:** a native prime-field complete-solve implementation
  or a new source-linked paired result.

## Resolution and history

2026-10-12: opened with structural evidence and the exact remaining cost
obligation. Unresolved.
