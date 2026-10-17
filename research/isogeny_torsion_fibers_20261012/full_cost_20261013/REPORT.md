# Fully timed native isogeny-fiber relation comparison at four sizes

2026-10-13, Conductor T-2486. This report covers the preregistered
[protocol](PROTOCOL.md), the retained [V2](REMEDIATION_V2.md) and
[V3](REMEDIATION_V3.md) timing remediations, and the exact raw records.
The degree-three formula is the `A=0` specialization of Vélu's map in
[Pizzo, Pomerance and Voight, equation (2.7)](https://math.dartmouth.edu/~carlp/3isog-070319.pdf).
It narrows the existing [complete-cost and independent-host problem](../../../problems/nonrational-isogeny-fiber-dlp-cost-gap.md)
but does not discharge its isolation, ecbench pricing, or second-host gates.

## Result and its unit

The native runner solved the same eight public source targets in four
prime-order subgroups through cold source relations, transported destination
relations, and matched source negation rho. Five measured rounds and one
warmup round produced 640 measured and 128 warmup arm records. Every record
recovered a scalar; the independent Rust checker reproduced relation
collection, row rank, factor-base log recovery, target descent, original
source scalar checks and destination checks. The exact identities, including
ICV1 and EC1, are in [identities.tsv](identities.tsv).

**Primary unit:** current-thread CPU milliseconds per complete cold solve.
The clock starts before each arm's own setup and ends after scalar checks.
It directly includes factor-base roots/projection, pair-table construction
and queries, modular rows, map evaluations, rho class work, allocation,
descent and verification. It removes time the process was descheduled, but
remains a shared-host L0 measurement. Phase wall clocks and native counts
are retained as an audit vector; their wall ratios are not used for a speed
claim. The eight targets are independent resampling units: each bootstrap
draw carries that target's five rounds together (10,000 fixed-seed draws).
The intervals below are paired percentile intervals, descriptive at L0.

| subgroup `r` | arm | verified | mean CPU ms | CPU / rho | paired 95% interval | CPU / source |
| ---: | :--- | ---: | ---: | ---: | :--- | ---: |
| 103 | negation rho | 40/40 | 0.007055 | 1.000 | — | 0.302 |
| 103 | source relations | 40/40 | 0.023328 | 3.307 | [2.318, 4.953] | 1.000 |
| 103 | transport relations | 40/40 | 0.024330 | 3.449 | [2.452, 5.014] | 1.043 |
| 103 | rho A/A | 40/40 | 0.005716 | 0.810 | [0.531, 1.300] | 0.245 |
| 439 | negation rho | 40/40 | 0.007010 | 1.000 | — | 0.105 |
| 439 | source relations | 40/40 | 0.066758 | 9.524 | [8.087, 11.382] | 1.000 |
| 439 | transport relations | 40/40 | 0.068613 | 9.789 | [7.470, 12.272] | 1.028 |
| 439 | rho A/A | 40/40 | 0.006435 | 0.918 | [0.787, 1.040] | 0.096 |
| 1,613 | negation rho | 40/40 | 0.014116 | 1.000 | — | 0.042 |
| 1,613 | source relations | 40/40 | 0.339743 | 24.068 | [20.222, 28.456] | 1.000 |
| 1,613 | transport relations | 40/40 | 0.359370 | 25.459 | [22.261, 28.871] | 1.058 |
| 1,613 | rho A/A | 40/40 | 0.011562 | 0.819 | [0.687, 0.974] | 0.034 |
| 6,427 | negation rho | 40/40 | 0.020155 | 1.000 | — | 0.012 |
| 6,427 | source relations | 40/40 | 1.691627 | 83.930 | [67.667, 104.252] | 1.000 |
| 6,427 | transport relations | 40/40 | 1.649107 | 81.820 | [69.169, 96.622] | 0.975 |
| 6,427 | rho A/A | 40/40 | 0.016618 | 0.824 | [0.732, 0.925] | 0.010 |

All four transport/source CPU intervals include 1: `[0.954,1.141]`,
`[0.911,1.136]`, `[0.903,1.249]`, `[0.875,1.081]` in ascending size.
The destination/source differences of −2.5% to +5.8% are smaller than the
8–19% aggregate CPU A/A imbalance. The transport/rho ratios are much
larger than this noise control at every measured size. This is a bounded
complete-cost observation for these exact algorithms and fields, not an
asymptotic conclusion or a result for deployed curves.

## Structure, boundary, and operation diagnostic

The selected `(p,r)` pairs are `(617,103)`, `(2633,439)`, `(9677,1613)`,
and `(38561,6427)`. On each `E/F_p: y²=x³+1`, the map
`φ(x,y)=((x³+4)/x², y(x³−8)/x³)` sends the three-point rational kernel
`{O,(0,1),(0,−1)}` to `E'/F_p: y²=x³−27`. Native preflight enumerated
both curves, checked order `p+1=6r`, all rational fibers, and absence of
an `F_p` scaling isomorphism. It checked all 381,924 ordered rational-point
pairs at `p=617`, 256 frozen sampled pairs at each larger size, and every
generator subgroup multiple. Each destination generator and target is the
image of its source counterpart. Because `p≡2 mod 3`, these fields have no
nontrivial rational cube root of unity for a six-class rho control.

The two-class generic walk-step floor is `sqrt(π/4)=0.886227` in
`S=(group additions+doublings)/sqrt(r)`. The table below is a **counted
group-operation diagnostic**: map field work, roots, lookups, modular work,
hashes and memory are outside `S` but inside the direct CPU totals above.

| `r` | arm | mean group ops | `S` | `S` / floor | group ops / rho |
| ---: | :--- | ---: | ---: | ---: | ---: |
| 103 | rho | 74.475 | 7.338 | 8.280 | 1.000 |
| 103 | source | 384.525 | 37.888 | 42.752 | 5.163 |
| 103 | transport | 399.050 | 39.320 | 44.367 | 5.358 |
| 439 | rho | 97.925 | 4.674 | 5.274 | 1.000 |
| 439 | source | 1,331.550 | 63.551 | 71.710 | 13.598 |
| 439 | transport | 1,395.400 | 66.599 | 75.149 | 14.250 |
| 1,613 | rho | 133.525 | 3.325 | 3.751 | 1.000 |
| 1,613 | source | 5,769.925 | 143.666 | 162.109 | 43.212 |
| 1,613 | transport | 6,002.875 | 149.466 | 168.654 | 44.957 |
| 6,427 | rho | 188.750 | 2.354 | 2.657 | 1.000 |
| 6,427 | source | 29,129.750 | 363.356 | 410.004 | 154.330 |
| 6,427 | transport | 28,816.200 | 359.445 | 405.590 | 152.669 |

The reported direct CPU unit prices native work that `S` leaves out. For
example, every destination target evaluates the map twice; its semantic
counter is `10` multiplications, `4` squarings and `2` inversions. The
actual evaluator computes inverses by exponentiation, whose instructions
are included in CPU time. [summary-v3.tsv](summary-v3.tsv) gives every
phase wall clock, native lift/lookup/row count, rho steps and process peak
RSS. The peak is process-wide, 25.9–27.2 MiB, not an arm-specific allocation.
Each relation arm rebuilds a 12-point factor base and its 78 pair entries
per target; setup is never divided among targets. The map and target are
public inputs to the destination method. Structural validation/preflight
was timed separately: 2.858 s, 0.0014 s, 0.0073 s, and 0.711 s in size
order. It includes research checks and target generation, and is excluded
from the primary single-target solve. Adding the entire verification
preflight to one destination solve is a conservative verified-from-scratch
upper bound, not a measured algorithmic map-construction cost. Candidate
band search was outside that preflight clock and remains separately unpriced.

## Retained attempts, replay, and validation

- **V1:** [raw](raw-v1.tsv), [summary](summary-v1.tsv), [run log](run-v1.log).
  The 640 measured answers and independent relation replays passed. Same-seed
  rho A/A operation counts matched exactly, but its wall ratios were
  `420.987`, `0.771`, `2.071`, `5.690`; V1 wall comparisons are censored by
  scheduling noise. The exact [V1 source snapshot](runner-v1.rs) matches the
  SHA-256 recorded before V2 edits.
- **V2:** [raw](raw-v2.tsv), [summary](summary-v2.tsv), [run log](run-v2.log).
  The 640 answers and independent relation replays passed. Thread-CPU A/A
  ratios were `0.827`, `0.924`, `0.754`, `0.817`; the third misses its frozen
  `[0.8,1.25]` gate. [V2 source](runner-v2.rs) was reconstructed from V3's
  sole arm-order change; recompilation reproduced all deterministic V2 rows.
- **V3:** [raw](raw-v3.tsv), [summary](summary-v3.tsv), [run log](run-v3.log).
  The balanced arm order passed its frozen CPU A/A gate at all sizes
  (`0.810`, `0.918`, `0.819`, `0.824`). The independent
  [checker receipt](independent-v3.log) reconstructs all 384 source and
  destination relation rows including warmups. The [same-host replay](replay-check-v3.log)
  matches all 768 deterministic answer/counter rows; clocks and RSS are
  excluded from that equality. V1 and V2 have corresponding source-snapshot
  [replay receipts](replay-check-v1.log) and [V2 receipt](replay-check-v2.log).

The three versions preserve their separate raw outputs and frozen protocols.
See [ARTIFACTS.sha256](ARTIFACTS.sha256) and the [validation record](validation/README.md)
for source/binary/input hashes and exact commands. The host was shared
Darwin arm64 (`25.6.0`), rustc `1.93.1`, L0; there was no configured
isolab client or second-host replay in this environment. V3 passes the
local complete-CPU and arithmetic gates. L2 wall admission, independent
host reproduction, and a native ecbench counted-price method are still
open. A growing-size scaling or speed claim remains unsupported by this
record.

| requirement | status | evidence or exact remaining gap |
| :--- | :--- | :--- |
| Four growing non-isomorphic pairs and exact identities | verified complete | `identities.tsv`; structural checks in `run-v3.log` |
| Paired source/destination relations and original-source recovery | verified complete | 640/640 measured rows; independent collection, rank, log and descent replay |
| Cold native work in one direct unit, including maps and setup | verified local L0 | full current-thread CPU nanoseconds per arm in `raw-v3.tsv`; process memory and wall phase vector retained |
| Same-host deterministic replay | verified complete | 768 matching answer/counter records in `replay-check-v3.log` |
| A/A CPU gate | verified complete at frozen threshold | V3 ratios 0.810–0.918; residual imbalance limits small-difference interpretation |
| L2 wall admission and independent-host replay | blocked by environment | no configured isolab client or second host in this session |
| Portable counted-unit price and iso-only attribution | not implemented in this round | unpriced native GAE primitives and scaling-model arm in [next gate](NEXT_GATE.md) |

## Proof map

The written non-isomorphism argument uses the classification of short
Weierstrass models with nonzero `b`: an `F_p` isomorphism between these two
models would require `u^6 = -27` in `F_p`. For the four chosen primes the
ratio is a quadratic nonresidue, as checked by exhaustive native residue
enumeration. The structural preflight also checked the degree-three map,
kernel, subgroup transport, and fiber multiplicity with the scope above.

The pinned [Lean module](lean/NonisomorphicFour.lean) formalizes the finite
sixth-power exclusion for `p=617,2633` as `ratio617,ratio2633`, and the
coprimality of each `r` with the kernel order three as
`coprime103,coprime439,coprime1613,coprime6427`. These are conditional
arithmetic components: the Weierstrass-isomorphism classification, the two
larger finite exclusions, and the elliptic-curve map/fiber theorem have no
Lean declaration in this round. The first full-module `decide` build hit a
recursion-depth limit; the raised-limit build was interrupted after a
bounded wait. Both attempts and their source versions are retained under
`lean/`. The final partial-module build and axiom inspection are recorded
in [validation](validation/README.md). None of these partial components is
presented as a Lean proof of the complete structural assertion.

## Visuals and decision

The [comparison graph](graph-comparison.svg) shows the complete V3 CPU
ratios; the [control graph](graph-controls.svg) shows why V1 and V2 were
retained but not used for the primary table. Both are embedded in the
[results PDF](report.pdf). The measured transport mechanism has not
reduced total cost against matched rho on these four fields. Its small
source-versus-destination differences do not clear the local A/A variation.
The earlier F617 counted-group graph remains unchanged: it is a frozen
single-size run in a different unit, and this round has its own source-linked
four-size CPU graph. The canonical ecbench frontier and boundary ledger
remain unchanged because this L0 CPU table has no portable counted-native
price or admitted isolation receipt.
The next discriminating run requires an L2 isolated host and independent
host replay of this frozen Rust workload, plus an isomorphism-only
destination control if attribution among coordinate models is sought.
Its concrete, unexecuted acceptance protocol is [NEXT_GATE.md](NEXT_GATE.md).
