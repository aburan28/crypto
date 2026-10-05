# F6-IC against inherited F4 and strong rho on an E_0 size ladder

Status: **preregistered, no outcomes measured.** This protocol, its four
`SPEC-n*.json` files and the decision rules below are committed before any
point in them is generated or any measured run begins.

**Disclosure.** While the adapter was being written, one implementation
smoke test ran on m = 13 with d = 5 and scratch seeds that are not in any
spec here (target seed 7, two targets). All six runs verified. The F4 and F6
arms issued identical query streams (36 and 17 queries, 16 and 10
relations). F6-IC cut Boolean reductions by about 2.4×, cut word XORs by
about 1.3×, and spent 17,784–55,488 charged geometric additions doing it.
Those runs are not evidence and are not committed. The decision rules
below were written before the smoke test and are unchanged by it. No curve, seed,
dimension rule, budget, cap or decision rule may change after an outcome is
seen. An amendment is a new protocol with new runs; every original row is
kept.

## Why this round

PR #1333 measured F6-IC, which closes branches of an inherited Boolean F4
search with exact factor-base geometry, against inherited F4 on one curve,
`icv1-f2m17-tm101-00378d4e`. On eight prepared targets it was 1.291–1.583×
faster online in wall time on an unisolated Mac. Its preregistered 2× gate
failed on every one. It never ran rho, never charged its work in ecbench's
counted unit, and measured one size, so it cannot say whether the gap is a
constant or grows with the field.

A constant factor cannot move the method against the floor (AGENTS.md §3):
inherited F4 costs 9,104.845× same-point rho online on that curve
(`research/ic_solver_online_20261003/confirmation/RESULT.md`), and 1.5× of
that is still about 6,000×. That figure is an extrapolation from those two
measurements, not a measurement. The only result that could matter is an
**exponent change**: geometric closure pruning a fraction of the Boolean
search that grows with the field degree `m`. This round tests that and
nothing else.

## Family and sizes

The ECC2K-130 family `E_0: y² + xy = x³ + 1` (`a = 0`) over `GF(2^m)`,
AGENTS.md §8b, at four prime degrees with a recorded prime-order subgroup.
None has a proper intermediate subfield over `GF(2)`.

| m | ICV1 slug | r | log₂ r | h | ord_m(2) | d = ⌈m/3⌉ | S_floor = √(π/4m) |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 13 | `icv1-f2m13-t181-515ee569` | 2,003 | 10.97 | 4 | 12 | 5 | 0.24579512472436318 |
| 19 | `icv1-f2m19-t797-b6cf2467` | 130,873 | 17.00 | 4 | 18 | 7 | 0.20331440047859611 |
| 23 | `icv1-f2m23-t5197-69e76b73` | 2,095,853 | 21.00 | 4 | 11 | 8 | 0.18479108808238470 |
| 31 | `icv1-f2m31-tm90707-c95f16f5` | 1,439,393 | 20.46 | 1,492 | 5 | 11 | 0.15917105461020273 |

Disclosures, per §8b:

- **Cyclotomic block.** At m = 31 the order of 2 mod 31 is 5, so the
  nontrivial cyclotomic block of `x^31 − 1` splits. At 13 and 19 it is
  irreducible. At 23, `ord_23(2) = 11` and it splits into two degree-11
  factors. F6-IC does not use invariant subspaces, but the factor base is a
  standard linear subspace, so the split is stated here.
- **Subgroup size.** `r` is not monotone in `m`: the m = 31 subgroup is
  smaller than the m = 23 one, and its cofactor is 1,492. Exponents of solver
  work are therefore fitted against `m`, which sets the descent system's
  shape. `S` is normalised by `√r` as always.
- **Not the gate size.** None of these is m = 83. Nothing here discharges the
  §8a gate, and no claim transfers to ECC2K-130 (m = 131).

The `a = 1` curves of #1333 (n9, n17) are not in this ladder. They are cited
only as the source of the F6 implementation and its prior numbers.

## Method under test and its control

Both IC arms run the native cold `ic.pipeline` on the same factor base,
query stream, linear algebra and verification. They differ in one place:
the PDP engine inside the decomposition oracle.

- **Factor base `koblitz-standard-subspace:dimension=d`.** This is #1333's
  `build_standard_subspace_factor_base(kc, d)`, the base its F4/F5/F6
  workers used. It has every point whose abscissa lies in the standard
  `d`-dimensional linear subspace. `d = ⌈m/3⌉` is fixed in advance. Three
  summands of `2^d` abscissae then cover about `2^{3d}/(6·2^m)` of the
  `x`-line: 0.67 at m = 13, 19 and 31, and 0.33 at m = 23.
- **Oracle `pdp3-koblitz`, three summands, Semaev `S_4` Weil descent,
  `max_degree = 3`, `node_budget = 8192`.**
  - Control arm: `engine=inherited-f4` calls #1333's `groebner_decompose`
    with `SolverEngine::InheritedF4 { max_degree: 3 }`.
  - Candidate arm: `engine=f6-ic` calls #1333's `groebner_decompose_f6_ic`
    with the same engine.
  - Every #1333 constant stays as committed. In particular, the
    single-fixed-summand closure applies only when the base has at most
    256 points. The predicted base sizes are about 32, 128, 256 and 2,048
    points, so **at m = 31 F6-IC runs without that closure**, with support
    pruning only. Whether the closure is active at m = 23 depends on the
    actual point count, which is recorded. Changing that cap would create a
    different candidate, not a tuning of this one.
- **Reference: `rho.signed_frobenius_strong`** (32 lanes, `dp_bits = 4`,
  step-cap factor 2000) on the same public point, with no table shared
  across targets.

Both IC arms are exact solvers on identical systems, with the same seed and
the same base. Their query streams are therefore **identical until an
outcome differs**. The adapter must show this: on every target the two arms
emit the same sequence of PDP outcomes and the same relations, unless one
arm exhausts its node budget on a query the other completes. Each such
divergence is recorded and listed. This pairing makes the per-query
reduction counts directly comparable.

## Workloads

Each size has its own spec, `SPEC-n{13,19,23,31}.json`. Each spec draws
eight public points per curve with `hash_to_subgroup_v1` from its own seed:

| m | target seed | measurement seed |
| ---: | ---: | ---: |
| 13 | 202610051301 | 202610051302 |
| 19 | 202610051901 | 202610051902 |
| 23 | 202610052301 | 202610052302 |
| 31 | 202610053101 | 202610053102 |

Every point is a separate one-target workload. No relation, log table or
collision table is shared across points or arms. No scalar is given to any
solver.

- **Rounds.** For m = 13, 19 and 23: one warm-up and three measured rounds,
  in alternating arm order. For m = 31: one measured round, no warm-up.
  Counts are deterministic, so repeats only serve wall time, and wall time
  at m = 31 is not used.
- **Deterministic caps.** `max_trials = 200000` relation trials and the
  node budget above. Exhausting either is a recorded failure, never a
  dropped row.
- **Wall-clock cap.** 3,600 s per run. A run that hits it is a recorded
  timeout.

## Accounting

The unit is ecbench's: `S = charged group-addition equivalents / √r`
(`docs/ecbench/README.md` §7).

- **Charged.** Base construction, every relation trial, every point
  addition the oracle performs, relation verification, linear algebra,
  target descent, and final scalar replay. That includes F6-IC's
  `geometric_group_additions`, which are affine point additions, so the
  adapter charges them on the `GroupOps` ledger.
- **Counted, uncharged.** The Boolean solver's word XORs, as counted by
  `koblitz_groebner::f4_word_ops_thread` around each call, under the unit
  string the repository's other inherited-F4 adapters use: "word XORs
  (elimination, specialisation and linear elimination only)". Also
  reductions, splits, propagations, geometric support checks, and residual
  lookups (exact and fast). These appear as `*_uncharged` counters. Every IC
  `S` here is therefore a **lower bound**, and its record says so
  (`cost.lower_bound = true`).

This asymmetry is the central accounting fact of the round. F6-IC's extra
work is charged, and the F4 reductions it removes are not. So a lower `S`
for F6 is a real saving in charged work, but a higher `S` does not show that
F6 costs more in total. The primary F4-versus-F6 comparison therefore uses
both quantities, never `S` alone.

**Sensitivity row, at m = 13 and 31 only.** `docs/ic/calibration.json`
pins `ns_per_word_xor / ns_per_add` for this family only at those two
degrees (its `koblitz` entries for `icv1-f2m13-t181-515ee569` and
`icv1-f2m31-tm90707-c95f16f5`). There,
and only there, the result also reports `S` with the solver's word XORs
priced at the pinned ratio. That ratio was calibrated on another host, so
the row is labelled as such, carries no fit, and never enters a headline.
No ratio is invented for m = 19 or 23.

## Hypothesis, statistic and decision

**H1 (exponent).** Geometric closure removes a fraction of inherited F4's
Boolean work that grows with `m`.

For target `t` at size `m`, let:

- `W4(t)` = total Boolean solver word XORs in the F4 arm, summed over all
  PDP calls, ordinary and target. A word XOR is the same unit at every
  `m`; a reduction is not, because its matrix grows with the field.
  Reductions are reported beside it;
- `W6(t)` = the same quantity in the F6 arm;
- `G6(t)` = F6-IC's geometric group additions.

The statistic is `y(t) = log₂(W4(t)/W6(t))`. Fit
`y = α + β·m` by ordinary least squares over all completed targets at all
four sizes. Get a 95% interval for `β` from 10,000 bootstrap resamples of
targets within each size, with resampling seed 202610050001.

- **Supported, a stage-level exponent lead.** `β`'s interval lies wholly
  above 0, **and** the fitted slope of `log₂ G6` per relation against `m` is
  no steeper than that of `log₂ W4` per relation. The first condition says
  F6 removes a growing share of the reductions. The second says the charged
  geometry it adds grows no faster than the reductions it replaces. The
  class is still a stage diagnostic: reductions are unpriced, so no
  whole-method advance can be claimed until a pinned reduction price
  exists. The next step would be pricing that unit, then the m = 83 gate.
- **Rejected, engineering.** `β`'s interval contains 0 or lies below it.
  The gain is a constant factor, F6-IC is classed **engineering**, and the
  F6 thread closes with a negative result.
- **Not decidable.** Fewer than six targets complete at any size. The
  ladder then falls back to m ∈ {7, 13, 19, 23}, using
  `icv1-f2m7-t13-616700dd` (r = 29) with d = 3 and target seed
  202610050701, and is labelled small-size-only. If that also fails, the
  round reports the failure and stops.

**Against the boundaries, at every size, for both IC arms.** Report:

- `S_IC,lower / S_floor`;
- `S_IC,lower / S_rho`, the median over targets, with the per-target
  minimum and maximum;
- the phase split of `S`;
- fitted exponents of `S_IC,lower`, `W4` and `W6` against `m`, next to
  rho's `S ≈ const`.

`S_IC` is a lower bound, so a size where `S_IC,lower ≥ S_rho` proves IC is
slower than rho there. No lower bound can prove the opposite.

**Abandon early** if, at m = 19, every target has `W6 ≥ W4`. Then F6 does
not reduce the work it targets on this family, and m = 23 and 31 are not
run.

**Inadmissible**:

- changing `d`, the node budget, the closure cap or the engines;
- choosing, dropping or replacing a target;
- pooling divergent query streams as if paired;
- treating a timeout, failure or unverified scalar as anything but a loss;
- quoting the reduction ratio as an end-to-end speedup;
- moving the boundary or the unit.

## Correctness and evidence

- Every completed run verifies `[k]G = Q` in the runner's own process.
- `ecbench verify --replay 2` passes for every run before any figure is
  cited.
- Before measurement, the adapter is cross-checked against #1333's code on
  every PDP query of a full run: for each query, the adapter's outcome and
  stats equal a direct call of the same #1333 function.
- Sessions are committed under `sessions/`, with the spec, records,
  receipts and audit.
- A `vs_rho` claim also needs an audit receipt from another host class
  (AGENTS.md §12). If none is available, that claim fails on that ground and
  the result says so.

## Host and wall time

The host is the session's Linux x86-64 cloud container. Its CPU model,
features, core count, memory, `rustc` version and commit go in the host
manifest. Counts are the metric. Wall time is a practicality note, used
only from uncontended runs at the isolation level the host achieves, and
read against an A/A interval. No wall-clock speedup is claimed. The result
names the hardware class it covers (x86-64 Linux VM) and the ones it does
not cover (Arm64, GPU).

## Deliverables

1. **Adapter PR.** The `koblitz-standard-subspace` base and `pdp3-koblitz`
   oracle, behind the Cargo feature `f6-ic-oracle`. Some frozen-source
   replay workflows rebuild the crate with pre-F6 snapshots of `src/lib.rs`,
   `Cargo.toml` and `koblitz_index_calculus.rs`, which do not declare the
   feature, so those builds exclude the adapter. The PR carries unit tests
   recovering known logarithms on a Koblitz curve with both engines, the
   cross-check above, and the README §7 charge row.
2. **Result PR.** Sessions, audit, `RESULT.md`, the table, and the
   dashboard updates §7/§7a require: a scoreboard panel, progress-chart
   points for each IC/rho ratio, the leaderboard, and the lab browser.
