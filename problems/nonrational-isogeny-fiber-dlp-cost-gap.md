# Complete DLP cost after non-rational-kernel isogeny transfer

- **Status:** open
- **Kind:** implementation and measurement gap; no mathematical
  nonexistence claim
- **Affected work:** [fiber report](../RESEARCH_ISOGENY_FIBERS_ECDLP.md),
  [2026-10-12 correction and bounded control](../research/isogeny_torsion_fibers_20261012/REPORT.md),
  [same-target paired round](../research/isogeny_torsion_fibers_20261012/paired_round/REPORT.md),
  [2026-10-13 isomorphism and six-class controls](../research/isogeny_torsion_fibers_20261012/paired_round/iso_control_20261013/REPORT.md)
- **Next action owner:** unassigned
- **Last updated:** 2026-10-13

## Problem

The paired round now executes a complete relation-based solve after
transport through an explicit `F_q`-isogeny with a non-rational geometric
kernel, on 16 public targets in two small prime-order pairs. It also runs
matched source-target negation rho and verifies every answer. A fully
priced, isolated, independently replayed cold-cost comparison showing
whether the transport helps ECDLP is still missing. The pairs have orders
19 and 211; extrapolation from them is not an acceptance criterion.
Both curves admit a rational order-three automorphism. A fixed six-class
rho control has now been executed and independently checked, but its
small-order restarts do not calibrate the best available reference.
Both destination models are also isomorphic to the sources over the
base field: `ψ_u(x,y)=(u²x,u³y)` with `u=2` at `p=13` and `u=14` at
`p=199`. The completed isomorphism-only solver arm showed a similar
destination relation workload at order 211. This does not assign the
recorded difference specifically to the degree-three map; a base-field
non-isomorphic pair and broader controls are still missing.

## Evidence and uncertainty

- **Observed:** [census.tsv](../research/isogeny_torsion_fibers_20261012/census.tsv)
  records 349 prime-order rows in a frozen scan of 1,014 curves; the
  [map check](../research/isogeny_torsion_fibers_20261012/map_check.txt)
  passed all 361 rational ordered-pair checks on the named `F₁₃` pair.
- **Observed in the paired round:** [raw-v4.tsv](../research/isogeny_torsion_fibers_20261012/paired_round/raw-v4.tsv)
  contains 48/48 verified arm records on 16 source targets. The degree-three
  maps passed 361 and 44,521 ordered-pair checks. At order 211, the
  source/transport counted-group totals were 5,832 / 4,803 = 1.21424;
  transport/rho was 4,803 / 690 = 6.96087. These are group-operation
  diagnostics with explicit unpriced native work, not complete-cost ratios.
  The [same-host replay](../research/isogeny_torsion_fibers_20261012/paired_round/replay-v4.tsv)
  is byte-identical. A separate [native arithmetic verifier](../research/isogeny_torsion_fibers_20261012/paired_round/verify.rs)
  checked all 48 records through its own point law; independent-host
  cost replay is absent. It also checked the two scaling isomorphisms on
  all 44,882 ordered rational-point pairs and verified `φ=ψ_u∘[m]`
  pointwise with `m=10 mod 19` and `m=29 mod 211`.
- **Proved:** the corrected cofactor statement depends on the order of
  the rational kernel, not the geometric degree. Its arithmetic step is
  compiled in [FiberCofactor.lean](../research/isogeny_torsion_fibers_20261012/lean/FiberCofactor.lean).
- **Interpretation:** the fixed pair-table relation collector used fewer
  trials after transport at order 211 (44.875 versus 61.50 on average).
  That is a coordinate-solver observation on these exact workloads.
- **Observed in the follow-on:** [raw-v3.tsv](../research/isogeny_torsion_fibers_20261012/paired_round/iso_control_20261013/raw-v3.tsv)
  has 80/80 verified arm records on the same 16 source targets. The
  [independent native checker](../research/isogeny_torsion_fibers_20261012/paired_round/iso_control_20261013/independent_v3.log)
  checked frozen targets, both destination images, all answers, phase
  sums, and six-class orbit multipliers. V1/v2 cap failures remain in
  [the report](../research/isogeny_torsion_fibers_20261012/paired_round/iso_control_20261013/REPORT.md).
  At order 211, degree-three and scaling-isomorphism destination arms
  used 359 and 358 relation trials and 4,871 and 5,046 counted group
  operations across eight targets. Six-class rho used 737 versus 690
  for negation rho in the same group-operation unit; at order 19 it
  required three charged restarts. V3 replay is byte-identical on the
  same host. These observations are bounded and do not price all native work.
- **Missing capability:** a native ecbench or equivalently audited paired
  session that prices map evaluation, factor-base lifts, table lookups,
  modular algebra and memory in a common unit, at growing sizes and with
  independent-host replay. A base-field non-isomorphic pair is needed to
  test a degree-three-specific relation effect. The current runner counts
  these activities but has no pinned conversion for all of them.
- **What this evidence does not establish:** an admitted end-to-end speed
  ratio, an asymptotic improvement, or a result for deployed curves or
  other fields.

## What would unblock it

- [x] Freeze and execute public same-target source/codomain pairs with
  explicit maps, model/field/subgroup identities, native relation solves,
  equal stop caps, and a matched rho arm, for the bounded control.
- [x] Execute scaling-isomorphism relations and a charged six-class rho
  control on the same frozen targets, with complete relation recovery,
  failure preservation, independent arithmetic verification, and
  byte-identical same-host replay. See the [follow-on](../research/isogeny_torsion_fibers_20261012/paired_round/iso_control_20261013/REPORT.md).
- [ ] Obtain verified logs for both solver arms, full phase accounting,
  **priced** native work, failed-attempt records, and independent-host replay;
  compare the degree-three arm with a scaling-isomorphism arm and a
  base-field non-isomorphic pair, calibrate both negation and six-class
  rho, then report a cold-cost ratio
  and interval in the repository's counted unit at at least four growing
  sizes. The bounded rounds have verified logs, phase counts and retained
  failed record versions, but do not meet the full pricing and replay gate.

Meeting these conditions would settle only the measured pair and workload.
It would not prove a general exponent change or a result for another curve.

## Next discriminating step

- **Completed mechanism control:** the [preregistered scaling-isomorphism
  test](../research/isogeny_torsion_fibers_20261012/paired_round/iso_control_20261013/REPORT.md)
  ran on the same 16 source targets with equal relation inputs/seeds/caps,
  rank, recovery, scalar verification, and both explicit maps. The two
  destination arms had similar order-211 relation yields, but their
  deterministic trajectories and total operation counts were not equal.
- **Completed bounded reference control:** the same follow-on ran a
  charged six-class rho arm with explicit orbit and coefficient
  canonicalization, verified all 16 scalars in v3, and preserved v1/v2
  cap failures. Stronger tuning and larger-order calibration remain open.
- **Next mechanism test:** preregister and construct at least one exact
  base-field non-isomorphic source/destination pair with an explicit
  verified map and prime-order target subgroup. Check complete relation
  solves, recovery, and both rho controls on paired public targets; this
  tests whether the present coordinate result persists when degree-one
  transport to the same model is unavailable.
- **Action:** extend native ecbench with a paired target transport path or
  an equivalently audited counted runner; preregister at least four
  increasing prime-order subgroup isogeny pairs, eight or more public
  targets per size, and five interleaved rounds. Include a base-field
  non-isomorphic pair to test a degree-three-specific effect.
- **Inputs and controls:** exact ICV1/EC1 pairs, executable isogeny,
  same source target transported to destination, source relation baseline,
  matched negation rho, six-class automorphism-aware rho, and A/A timing
  control if wall time is used.
- **Limits:** pin conversions for map, roots, lookups and modular work;
  include memory, construction and failed attempts; use ecbench's
  isolation and independent-host receipt gates. Stop on correctness or
  accounting failure and retain it.
- **Outcomes:** a verified complete-cost ratio below rho at growing sizes
  would motivate expansion; a stage-only improvement remains a stage
  diagnostic; a timeout, invalid map or missing solver support leaves the
  performance question open and identifies the exact cause.
- **Recheck trigger:** a native prime-field complete-solve implementation
  or a new source-linked paired result.

## Resolution and history

2026-10-12: opened with structural evidence and the exact remaining cost
obligation. Later the same day, the paired round discharged the bounded
transport, relation, recovery and verification control on 16 public
targets. The full-cost and independent-host gates remain open.
The boundary review identified the additional six-class rho control:
`ζ=3` and `ζ=92` satisfy `ζ³=1` over the two respective fields, so
the selected negation-only rho is an exact-target reference but does
not cover every rational automorphism of these curves.
The final structural audit found the two public scaling isomorphisms;
the isomorphism-only relation arm is now the cheapest discriminating
experiment. The v4 measurement and its frozen inputs remain unchanged.
