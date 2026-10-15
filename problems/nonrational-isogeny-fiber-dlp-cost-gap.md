# Complete DLP cost after non-rational-kernel isogeny transfer

- **Status:** open
- **Kind:** implementation and measurement gap; no mathematical
  nonexistence claim
- **Affected work:** [fiber report](../RESEARCH_ISOGENY_FIBERS_ECDLP.md),
  [2026-10-12 correction and bounded control](../research/isogeny_torsion_fibers_20261012/REPORT.md),
  [same-target paired round](../research/isogeny_torsion_fibers_20261012/paired_round/REPORT.md),
  [2026-10-13 isomorphism and six-class controls](../research/isogeny_torsion_fibers_20261012/paired_round/iso_control_20261013/REPORT.md),
  [2026-10-13 non-isomorphic rational-fiber control](../research/isogeny_torsion_fibers_20261012/nonisomorphic_pair_20261013/REPORT.md),
  [four-size cold native CPU comparison](../research/isogeny_torsion_fibers_20261012/full_cost_20261013/REPORT.md),
  [three-torsion orbit support screen](../research/isogeny_torsion_fibers_20261012/torsion_orbit_20261013/REPORT.md),
  [quadratic-extension cold DLP comparison](../research/isogeny_torsion_fibers_20261012/extension_field_20261014/REPORT.md),
  [exact pairing-transfer cost result](../research/isogeny_torsion_fibers_20261012/pairing_transfer_20261015/RESULTS.md),
  [2026-10-15 χ polynomial run](../research/isogeny_torsion_fibers_20261012/chi_polynomial_20261015/reports/run-20261015T1651Z/REPORT.md),
  [public-only direct-root repeat](../research/isogeny_torsion_fibers_20261012/chi_polynomial_20261015/reports/run-20261015T1851Z/REPORT.md),
  [different-j relation and fiber controls](../research/isogeny_torsion_fibers_20261012/different_j_20261015/reports/run-20261015T1713Z/REPORT.md),
  [different-j full-image factor-base gate](../research/isogeny_torsion_fibers_20261012/different_j_full_image_20261015/RESULTS.md),
  [different-j 60-target reuse and whole-process result](../research/isogeny_torsion_fibers_20261012/different_j_reuse_20261015/RESULTS.md),
  and [three-arm pairing batch cost report](../research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/reports/run-20261015T1740Z/REPORT.md)
- **Next action owner:** unassigned
- **Last updated:** 2026-10-15

## Problem

The paired round now executes a complete relation-based solve after
transport through an explicit `F_q`-isogeny with a non-rational geometric
kernel, on 16 public targets in two small prime-order pairs. It also runs
matched source-target negation rho and verifies every answer. A four-size
native cold-solve comparison now prices each arm directly in current-thread
CPU time, including setup, collection, recovery, transfer, descent and
verification. It has same-host exact replay and an independent relation
checker. A portable counted-unit price, L2 isolation, and independent-host
replay are still missing. The earlier pairs have orders 19 and 211;
extrapolation from them is not an acceptance criterion.
Both curves admit a rational order-three automorphism. A fixed six-class
rho control has now been executed and independently checked, but its
small-order restarts do not calibrate the best available reference.
Both destination models are also isomorphic to the sources over the
base field: `ψ_u(x,y)=(u²x,u³y)` with `u=2` at `p=13` and `u=14` at
`p=199`. The completed isomorphism-only solver arm showed a similar
destination relation workload at order 211. This does not assign the
recorded difference specifically to the degree-three map. A later
base-field non-isomorphic `F617` pair with a three-point rational kernel
has now been solved in both relation arms. Three more non-isomorphic pairs
and direct cold CPU timing have since been added. A separate isomorphic
model arm on the exact four-size targets remains an attribution control; it
cannot map to their fixed non-isomorphic destinations.

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
- **Observed on a non-isomorphic pair:** the
  [F617 report](../research/isogeny_torsion_fibers_20261012/nonisomorphic_pair_20261013/REPORT.md)
  records two order-618 models without a base-field isomorphism, an
  explicit degree-three map with a rational kernel of size three, eight
  checked target fibers, and 24/24 independently verified
  source-relation, destination-relation and matched rho rows. Source and
  destination relation totals were 3,312 and 3,794 counted group operations
  on the same targets; rho used 664. The destination/source ratio is
  1.145531 in that unit. Clean v2 execution and replay were byte-identical;
  v1's post-child timing-wrapper failure remains preserved. This is one
  small subgroup order on a shared host, with native costs unpriced.
- **Observed at four sizes:** the [native run](../research/isogeny_torsion_fibers_20261012/full_cost_20261013/REPORT.md)
  used subgroups 103, 439, 1,613 and 6,427 on explicit base-field
  non-isomorphic degree-three pairs. Eight public targets per size, five
  measured rounds and one warmup produced 640/640 verified measured arm
  records. An independent Rust checker reconstructed all 384 relation
  records per version including warmups; same-host replay matched 768
  deterministic answer/counter rows. Complete cold current-thread CPU
  transport/rho ratios were 3.449, 9.789, 25.459 and 81.820;
  transport/source ratios were 1.043, 1.028, 1.058 and 0.975. All four
  paired transport/source intervals contain one, within residual 8–19%
  rho A/A CPU imbalance. V1 wall and V2 CPU control failures are preserved.
  These are L0 direct-time observations, without admitted isolated wall.
- **Observed in the torsion-orbit screen:** the [exact native counts](../research/isogeny_torsion_fibers_20261012/torsion_orbit_20261013/raw-v1.tsv)
  found two-sum subgroup support 64, 74, 75 and 78 for both twelve
  subgroup points and all 36 of their rational three-torsion lifts.
  Thirty-six distinct subgroup points reached 103, 351, 559 and 642 at
  the same four subgroup sizes. A complete scalar index agreed with group
  addition. The conditional support identity compiled in Lean without
  axioms. This is a relation-hit diagnostic, with no new full-solve timing.
- **Constructed equivariant coordinate:** on the same source family,
  `chi=(y+1)/(y−1)` obeys `chi(P+(0,1))=1/(1−chi(P))` and takes distinct
  values `1, infinity, 0` on the three kernel points. The
  [native verifier](../research/isogeny_torsion_fibers_20261012/torsion_orbit_20261013/equivariant-v1.tsv)
  passed both translation identities on all 51,492 rational points across
  the four fields and the order-three projective action on all 51,492
  projective-line points. Denominator-cleared algebraic steps have
  [`IDC1h41adc1f4e1c5522c`](../research/isogeny_torsion_fibers_20261012/torsion_orbit_20261013/group-law.identity.certificate.json)
  and [`IDC1ha0c3d83775e32754`](../research/isogeny_torsion_fibers_20261012/torsion_orbit_20261013/mobius-conjugacy.identity.certificate.json)
  certificates; the report gives the curve-law derivation.
  Lean checks the induced linear action's order but lacks the curve-level
  bridge. The bounded χ polynomial decomposition system has since been
  built and measured; see the next bullet.
- **Observed in the native χ polynomial run:** the
  [exact census](../research/isogeny_torsion_fibers_20261012/chi_polynomial_20261015/reports/run-20261015T1635Z/REPORT.md)
  matched independent pair enumeration in 22/22 cases. The
  [paired run](../research/isogeny_torsion_fibers_20261012/chi_polynomial_20261015/reports/run-20261015T1651Z/REPORT.md)
  recovered and independently replayed every χ arm scalar in 12/12
  rounds, while projected-pair and χ had identical relation counts,
  rank, and query trials in each round. The native χ implementation scans
  polynomial candidates and filters them by exact group addition. A
  subsequent [direct S3 root and quotient exactness run](../research/isogeny_torsion_fibers_20261012/chi_polynomial_20261015/reports/run-20261015T1827Z/REPORT.md)
  matched 22/22 independent decomposition sets; its
  [paired relation run](../research/isogeny_torsion_fibers_20261012/chi_polynomial_20261015/reports/run-20261015T1828Z/REPORT.md)
  recovered each new arm's scalar in 12/12 rounds with unchanged
  relation/query/rank outcomes. The independent checker replayed 123/132
  arm-round scalars; the other nine were raw-pair cap exhaustions. The
  direct χ-coordinate resultant remains unimplemented. Local measured
  pipeline CPU for each polynomial arm exceeds rho on the
  four frozen cells, but the formal counted `S` for the polynomial
  primitive is still incomplete.
- **Observed with public-only solver input:** the
  [V3 custody gate](../research/isogeny_torsion_fibers_20261012/chi_polynomial_20261015/reports/run-20261015T1850Z/REPORT.md)
  rejected an answer-bearing solver input with exit 101, then matched
  independent signed-pair enumeration in 22/22 target sets with 144
  decompositions and 2,992 pair/map checks. The
  [V3 full relation run](../research/isogeny_torsion_fibers_20261012/chi_polynomial_20261015/reports/run-20261015T1851Z/REPORT.md)
  used a framework path without a planted answer. Both direct root arms
  solved 12/12 rounds; independent replay checked 123/132 recovered
  scalar rows across all arms, and nine raw-pair rows hit their inherited
  cap. All 132 normalized relation, rank, trial, recovery and status rows
  matched V2. Paired source-root/rho pipeline CPU ratio medians were
  13.237, 8.254, 13.427 and 21.844 on the four frozen cells. The
  redaction and checker processes read the private fixture; the solver
  did not. The frozen protocol's narrower custody sentence is corrected
  in the [versioned note](../research/isogeny_torsion_fibers_20261012/chi_polynomial_20261015/PUBLIC_ONLY_PROTOCOL_V3_CORRECTION.md).
- **Observed on a different-j extension-field pair:** the
  [explicit degree-three map](../research/isogeny_torsion_fibers_20261012/different_j_20261015/reports/run-20261015T1713Z/REPORT.md)
  connects `j=7` and `j=11` curves over `F_(67²)`, with order-61 target
  subgroup. Source, image and destination-direct relation arms all
  reached rank nine and recovered the same scalar in three rounds; the
  independent checker verified 15/15 arm-round rows. The
  [exhaustive fiber audit](../research/isogeny_torsion_fibers_20261012/different_j_20261015/reports/run-20261015T1738Z/REPORT.md)
  found 4575 source points, 1525 image points, exactly three preimages
  per image and 61 distinct target-subgroup images, with 3721/3721
  subgroup homomorphism checks. Thus the rational fiber multiplicity
  did not add relation yield on this projected subgroup workload.
- **Observed with a full rational-image destination base:** the
  [separately frozen V2 run](../research/isogeny_torsion_fibers_20261012/different_j_full_image_20261015/RESULTS.md)
  enumerated all 1,525 rational image points, selected eight projected
  columns from 21 deterministic candidates, and replayed 18/18 scalar
  rows across four relation and two rho arms. All four bases had 56/60
  nonzero two-sum support; full-image and projected-image support shared
  52 sums and had four unique sums each. The three-round paired
  full-image/rho local CPU ratio median was 70.281, with a 98.532%
  median full-image setup share. The custom timer omits common parsing,
  model construction and process startup; portable counted-`S` remains
  open. The follow-on distinct-target reuse gate is now executed below.
  The checked Lean finite-support
  proof uses frozen scalar-coordinate lists; a faithful curve-level
  point-to-scalar and Vélu bridge remains open.
- **Observed over every distinct target in the frozen subgroup:** the
  [V3/V4 reuse study](../research/isogeny_torsion_fibers_20261012/different_j_reuse_20261015/RESULTS.md)
  solved all 60 nonidentity order-61 public targets in seven arms and
  three rounds per version. Each independent native replay checked
  1,260/1,260 scalars and 10,960 relation witnesses. At `N=60`, the
  V3 full-image cold native timer median was 11,049.667 µs against
  643.208 µs for reusable source BSGS; V4 whole-child CPU medians were
  15,351 and 6,148 µs. The paired V4 full-image/BSGS ratio median was
  2.483. V3 prefix timing omits startup and final stdout; V4 includes
  those costs for the full batch only. Shared-host L0 timing and this
  one group do not establish a portable or growing-size crossover.
- **Observed in the pairing batch:** the
  [V6 three-arm run](../research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/reports/run-20261015T1740Z/REPORT.md)
  replayed 594/594 scalar outputs on 66 distinct order-67 `F_401`
  targets. At `N=64`, pairing/direct-prime-rho whole CPU ratios were
  2.523, 1.620 and 3.000. Direct prime-field rho beat pairing in all
  three L0 rounds, even though pairing was below the weaker extension-rho
  arm in two V6 rounds and all three earlier V5 rounds. The
  [V6B same-input repeat](../research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/reports/run-20261015T1752Z/REPORT.md)
  started after the release tests exited, replayed another 594/594
  scalars, and again found direct rho faster in all three rounds; paired
  pairing/direct ratios were 1.944, 2.470 and 1.992. A portable
  counted-unit comparison remains open.
- **Observed against reusable source BSGS:** the
  [V7 cold three-arm run](../research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/reports/run-20261015T1805Z/REPORT.md)
  rebuilt a nine-entry direct `F_401` baby table each round and replayed
  594/594 pairing, direct-rho and direct-BSGS scalars. At `N=64`, paired
  pairing/direct-BSGS complete CPU ratios were 5.190, 4.178 and 5.749;
  source BSGS was the lower-cost arm at every measured prefix in all
  three shared-host L0 rounds. This rules out attributing a multi-target
  CPU advantage to pairing against the available direct methods on this
  fixed cell, while leaving counted-unit and growing-size questions open.
- **Missing capability:** a native ecbench or equivalently audited paired
  session that prices map evaluation, factor-base lifts, table lookups,
  modular algebra and memory in a portable counted unit, with L2 isolation
  and independent-host replay. Direct CPU includes those costs locally but
  does not pin their conversion to ecbench's counted `S` unit. The four
  pairs remove a public base-field scaling isomorphism; a same-target
  isomorphism-only model arm is still needed for attribution.
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
- [x] Construct and verify a base-field non-isomorphic degree-three pair
  with rational kernel size three, exact target fibers, subgroup projection,
  complete paired relation recovery, original-source scalar checks,
  matched negation rho, and independent native verification. See the
  [F617 report](../research/isogeny_torsion_fibers_20261012/nonisomorphic_pair_20261013/REPORT.md).
- [x] Execute four increasing non-isomorphic pairs with paired relation
  recovery and matched rho, independent arithmetic replay, retained
  failures, full phase accounting, and complete local cold CPU pricing.
  See the [four-size report](../research/isogeny_torsion_fibers_20261012/full_cost_20261013/REPORT.md).
- [ ] On those same targets, add a separate base-field isomorphic model arm
  and a duplicate relation A/A control; admit an L2 isolated run and a
  second-host exact replay. Pin all native-work conversions for a portable
  complete `S` price and interval, including construction and memory with
  exact setup reuse. The [next protocol](../research/isogeny_torsion_fibers_20261012/full_cost_20261013/NEXT_GATE.md)
  gives fixed inputs, stop conditions and receipts.
- [ ] Complete a Lean theorem for `chi(P+T)=1/(1−chi(P))` as a projective
  identity on `E:y²=x³+1`, including `O,T,2T`, with a faithful elliptic-curve
  group-law definition and checked curve/subgroup assumptions. The current
  [proof map](../research/isogeny_torsion_fibers_20261012/torsion_orbit_20261013/PROOF_MAP.md)
  has conditional algebra and polynomial-step certificates, but no
  curve-level bridge. Completing it would establish the coordinate identity
  under the theorem's field assumptions; it would not establish an oracle
  or a DLP speed advantage.
- [x] Price the implemented quadratic-extension projected-pair, χ orbit-index,
  degree-three map, public scaling-isomorphism map and matched rho controls
  cold on four frozen groups, with exact support census and independent
  native arithmetic replay. The
  [2026-10-14 report](../research/isogeny_torsion_fibers_20261012/extension_field_20261014/REPORT.md)
  records 87 verified V4 arm rows and nine capped raw-pair rows; all four
  projected arms had the same relation outcomes and the χ index had higher
  median CPU cost than the plain projected table at all sizes. The V3
  seven-arm result, with 75 verified rows, remains separately retained.
- [x] Compare a deterministically selected full rational-image base with
  source, projected-image and destination-direct bases on the frozen
  different-j pair. The [V2 report](../research/isogeny_torsion_fibers_20261012/different_j_full_image_20261015/RESULTS.md)
  records exact support membership, 18/18 independently replayed scalar
  rows and a charged local timer interval over three paired rounds.
- [x] Price that full-image construction over all 60 distinct reusable
  targets against direct source BSGS and rho in native CPU. V4 includes
  input/model validation, process startup, relation recovery, final
  verification and output; memory is recorded separately as peak RSS.
  The [V3/V4 result](../research/isogeny_torsion_fibers_20261012/different_j_reuse_20261015/RESULTS.md)
  preserves exact receipts and independent scalar replay. A portable
  counted-unit price, larger-order pair, isolation and second-host
  performance replay remain open; the native CPU result does not
  discharge those separate acceptance conditions.
- [x] Construct and cold-price an exact pairing transfer on the first two
  frozen quadratic-extension targets, with a certified evaluation field,
  selected torsion input, separated H search/Miller/final exponent/field-DLP
  costs, same-target negation rho, peak RSS, and separate native scalar replay.
  The [V2 report](../research/isogeny_torsion_fibers_20261012/pairing_transfer_20261015/RESULTS.md)
  records six checked pairing solves, 18 public-multiple pairing checks,
  and six rho plus six rho A/A controls. Its cold solve plus replay median
  paired CPU/rho ratios were 2,666.48 and 660.71 at orders 19 and 61.
  This discharges the bounded two-cell construction/cost subtrack; growing
  sizes, a best finite-field solver, L2 isolation, and independent-host
  performance replay remain open.
- [x] Construct and price the bounded nonfactorable χ polynomial-pair
  decomposition oracle and complete relation recovery, with exhaustive
  pair-enumeration and independent scalar replay. See the
  [run report](../research/isogeny_torsion_fibers_20261012/chi_polynomial_20261015/reports/run-20261015T1651Z/REPORT.md).
- [x] Implement and exhaustively check direct Semaev roots and the
  degree-three quotient polynomial control, then complete the frozen
  paired relation solves and independent replay. See the
  [V2 exactness](../research/isogeny_torsion_fibers_20261012/chi_polynomial_20261015/reports/run-20261015T1827Z/REPORT.md)
  and [solve](../research/isogeny_torsion_fibers_20261012/chi_polynomial_20261015/reports/run-20261015T1828Z/REPORT.md) reports.
- [x] Repeat the direct-root and quotient controls with public-only
  solver input, reject a private answer-bearing input, and replay recovered
  scalars in a separate process. See the
  [V3 custody](../research/isogeny_torsion_fibers_20261012/chi_polynomial_20261015/reports/run-20261015T1850Z/REPORT.md)
  and [solve](../research/isogeny_torsion_fibers_20261012/chi_polynomial_20261015/reports/run-20261015T1851Z/REPORT.md) reports.
- [ ] Implement a direct resultant/root solver in the χ coordinate,
  derive a portable counted cost for each polynomial primitive, and
  repeat on isolated and independent hosts. The earlier χ finite scan
  and the direct `S3` abscissa root are distinct implementations; do not
  transfer either one's cost to the missing χ resultant.

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
- **Completed non-isomorphic mechanism test:** the
  [F617 control](../research/isogeny_torsion_fibers_20261012/nonisomorphic_pair_20261013/REPORT.md)
  used an explicit rational-kernel degree-three map, a prime-order target
  subgroup, both relation solves, recovery, and matched negation rho on
  eight public targets. No base-field degree-one scaling exists for this
  pair. Destination relation work was higher in the aggregate and lower
  on three of eight individual targets. This observation is exact for the
  frozen workload.
- **Action:** execute the [frozen next gate](../research/isogeny_torsion_fibers_20261012/full_cost_20261013/NEXT_GATE.md)
  on an available isolated host: same four pairs and targets, the public
  scaling model, duplicate relation and rho controls, portable native
  conversion, and independent-host replay.
- **Separate torsion lead:** the free three-lift relation-yield multiplier
  failed its exact pair-table test. The subsequent
  [quadratic-extension run](../research/isogeny_torsion_fibers_20261012/extension_field_20261014/REPORT.md)
  completely priced a χ orbit *index* on frozen targets and found unchanged
  projected relation outcomes with higher local CPU cost. A finite census
  found two extra raw support points from torsion shifts at `p=787`, while
  projected support stayed fixed. The native χ polynomial scan then
  reproduced projected relation outcomes in all twelve rounds at higher
  local CPU cost. The direct source `S3` and quotient-root controls have
  since been run on the same inputs and memory cap, also with unchanged
  projected relation outcomes. A χ-coordinate resultant solver and
  portable counted primitive pricing remain open.
- **Quadratic-extension map attribution:** all four exact ordinary fields
  admit `u²=−3` in the base field, giving a checked public scaling map
  `ψ_u:E→E'` to the degree-three destination. The V4 paired run included
  that map arm and independently replayed all 87 recovered scalars. All
  four projected arms had identical relation/query/rank outcomes on the
  twelve paired cells. The new `p=67`, different-j pair supplies the
  separately requested non-isomorphic extension-field destination and
  paired solve; its full rational fiber audit is linked above.
- **Pairing route, bounded result:** [V2 raw and phase records](../research/isogeny_torsion_fibers_20261012/pairing_transfer_20261015/RESULTS.md)
  establish an exact finite-field transfer and independently checked
  scalars at orders 19 and 61. H search used 26,469/27,169 and
  39,115/39,673 selection extension multiplications. The one-target cold
  route remained above rho. The 66-target order-67 batch now charges
  reusable H and the finite-field table against actual distinct targets.
  Its direct-prime-rho arm beat pairing in three V6 shared-host rounds at
  `N=64`; a later cold direct source-BSGS arm was faster still in V7.
  Portable counted pricing, isolation and larger-size scaling remain the
  next cost gates.
- **Inputs and controls:** exact ICV1/EC1 pairs, executable isogeny,
  same source target transported to destination, source relation baseline,
  matched negation rho, six-class automorphism-aware rho where rationally
  available, and A/A timing control if wall time is used.
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

2026-10-13: the isomorphism-only arm and charged six-class rho control
completed on the frozen F13/F199 targets. A separate preregistered F617
pair then supplied a non-isomorphic destination, rational kernel of size
three, and eight explicitly checked three-point target fibers. All 24
same-target relation/rho rows verified, with destination/source counted
group operations 3,794/3,312 and same-host byte-identical replay. Full
native pricing, growing sizes and independent-host replay remain open.

2026-10-13: four growing non-isomorphic pairs received a native complete
cold current-thread CPU comparison with 640 independently verified
measured arm records and preserved V1/V2 timing-control failures. At
orders 103, 439, 1,613 and 6,427 the destination/rho CPU ratios were
3.449, 9.789, 25.459 and 81.820. The remaining gap is portable counted
pricing, an iso-only same-target attribution control, L2 admission and
independent-host replay; the problem stays open.

2026-10-13 local: a frozen exact screen checked rational three-torsion
orbits on the same four pairs. Orbit-expanded bases had identical two-sum
support to their twelve underlying subgroup points. The existing research
note's unconditional factor-three collection interpretation was corrected;
the extension-field polynomial-oracle mechanism remains open. A further
frozen exact check established the nonfactorable degree-three coordinate
`chi=(y+1)/(y−1)` and its order-three torsion action on every rational
point of the four curves. This is an explicit oracle construction lead;
it carries no measured DLP cost yet.

2026-10-14: a frozen four-size quadratic-extension DLP run priced an exact
cofactor-projected pair table, χ orbit index and explicit degree-three
isogeny control against same-target negation, A/A and six-class rho. A
separate Rust process replayed 75 recovered scalars; nine raw-pair rows
reached their cap. The χ index matched projected relation outcomes and
cost more CPU on all four sizes. Exact torsion shifts added two raw support
points on one frozen base, with no extra projected support. The direct
Semaev and χ polynomial oracles, pairing transfer, peak-memory equality,
L2 isolation and independent-host cost remain open.

2026-10-14 V4: a preregistered scaling check found explicit public
isomorphisms on all four quadratic-extension pairs (`u=27,57,28,264`). The
eight-arm rerun retained the same four targets, seeds and caps and placed
that scaling arm next to the degree-three arm. Independent arithmetic
checked 87 recovered scalars across 96 slots; nine raw-pair rows capped.
A native paired audit found equal relation, query, rank and answer outcomes
for the projected, χ-index, isogeny and scaling arms in all twelve cells.
The fixed destination therefore cannot isolate a degree-three-specific
support effect. For any odd `p≠3`, every nonzero `F_p` element is a square
in `F_(p²)`, so this fixed `x³+1 → x³−27` pair always has such a scaling
over the quadratic extension. The χ polynomial-system oracle, a different
non-isomorphic quadratic-extension destination, portable memory/cost pricing,
isolated host and pairing route remained unresolved at that date.

2026-10-15: the separately preregistered native pairing route completed
on the first two frozen quadratic-extension cells. The solver received no
answer scalar, and six separate native arithmetic replays checked the
returned values. V1 retained a combined selection/pairing phase; V2
separated candidate search, Miller, and final exponentiation without
changing deterministic primitive totals. The V2 fully validated paired
CPU/rho ratios were 2,781.81 and 698.15. Shared-host L0 noise and only
two subgroup sizes leave scaling, reusable setup, and isolated timing open.

2026-10-15 direct-root V2: both source and quotient `S3` root arms
matched independent enumeration on all 22 smallest-cell targets and
recovered the scalar in all 12 paired rounds each. Their relation
counts, useful rank and trials matched the projected-pair control in
every round. The source-root/negation-rho paired CPU ratio medians were
15.474, 6.014, 14.669 and 30.631 across the four frozen cells. These
shared-host ratios include the measured pipeline but exclude common
input parsing and calibration. The [proof map](../research/isogeny_torsion_fibers_20261012/chi_polynomial_20261015/PROOF_MAP_ROOTS_V2.md)
records a checked conditional relation-transfer theorem and the open
concrete polynomial/curve bridge. The destination full-image factor-base
and χ-coordinate resultant questions remain open.

2026-10-15 full-image V2: the separately frozen different-j destination
base selected eight projected columns after enumerating all 1,525
rational image points. Four relation arms solved the same order-61 target
in three rounds each, and two rho arms did likewise; a separate native
checker replayed all 18 scalar rows. The source and projected-image
two-sum support sets coincided, while the full-image set exchanged four
of 56 nonzero sums with the projected-image set. The full-image cold
local timer interval had a paired rho ratio median of 70.281, dominated
by construction. Reusing that construction over distinct targets,
portable counted pricing and the direct χ-coordinate resultant remain
open obligations.
