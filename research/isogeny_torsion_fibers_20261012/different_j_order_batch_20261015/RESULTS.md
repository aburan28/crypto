# Complete cold cost of the order-229 different-\(j\) fiber base

**Status:** V1 execution complete and independently checked on the
frozen 64-target `r=229` instance; PR review and CI remain separate
gates. The [descriptive run report](reports/run-20261015T2114Z/REPORT.md),
[PDF](reports/run-20261015T2114Z/REPORT.pdf), [native
summary](SUMMARY_V1.json), [graph](GRAPH_V1.svg), and [raw
records](RUN_V1.jsonl) carry the measurements. The source and
destination are the distinct-`j` `F_(211²)` models from the
[full rational-fiber audit](../different_j_order_screen_20261015/FIBER_REPORT.md).

## Request and measured boundary

The user requested explicit isogeny and torsion-fiber transport,
relation collection, recovery and same-target verification, fully
priced against generic ECDLP controls. The generic negation-orbit
collision floor for a single order-229 target is
`S_floor = sqrt(pi/(2·2)) = sqrt(pi/4) ≈ 0.8862`, with
`S = total equivalent group operations/sqrt(229)`. The two measured
complete direct references are reusable source BSGS and negation rho
on the exact same public target prefixes. The method's portable
counted `S` is **bounded/unpriced** because construction, hash/table
work, field inversions, linear algebra, startup and output do not yet
have an audited common primitive conversion. The native
[summary](SUMMARY_V1.json) gives an explicitly local CPU-to-addition
conversion using a median 197.50791 CPU ns per repeated source add;
that is a host-specific diagnostic and is not compared to the generic
floor as a portable operation count.

The one table reports **complete cold child CPU in microseconds**,
medians of three rotated rounds. These totals include map/model
construction, every scan and failed query, table setup, relation
collection, elimination, recovery, scalar verification, process
startup and output. Each target-count column is a fresh child with
the same frozen prefix; none inherits a preceding child's setup.

| Arm | N=1 | N=4 | N=16 | N=64 |
| --- | ---: | ---: | ---: | ---: |
| Source projected relations | 4,771 | 4,759 | 9,535 | 20,774 |
| Isogeny-image relations | 5,999 | 5,778 | 8,575 | 20,042 |
| Destination-direct relations | 5,039 | 6,229 | 11,332 | 22,054 |
| Full rational-image relations | 56,482 | 64,045 | 53,560 | 72,774 |
| Reusable source BSGS | 4,529 | 5,330 | 6,328 | 6,383 |
| Source negation rho | 5,030 | 6,333 | 4,797 | 7,653 |
| Rho A/A repeat | 5,248 | 4,980 | 6,088 | 7,357 |

At `N=64`, the **paired** full-image/BSGS child-CPU ratios were
11.037, 11.139 and 11.401 (median **11.139**); paired
full-image/rho ratios were 8.156, 9.509 and 10.177 (median **9.509**).
At all four prefixes and in all three rounds, both direct controls
used less whole-child CPU than the full-image arm. These are exact
ratios of the frozen L0 child meters, not isolated wall-time or
growing-size results. The A/A repeat and min–max whiskers in the
[graph](GRAPH_V1.svg) show the shared-host variation.

## Fiber and relation signal

The independently checked map has 44,655 rational source points,
14,885 distinct rational images and exactly three rational preimages
per image. The selected full-image factor base uses eight signed
subgroup columns and public source preimage witnesses; its
construction re-enumerates and verifies the full rational image in
each cold child. At `N=64`, that setup alone used **3,667,277**
base-field multiplications and **89,636** inversions per child. The
full-image median complete counters were 5,138,577 multiplications
and 143,637 inversions, against 43,313/1,641 for BSGS and
160,515/6,345 for rho. Maximum full-image child RSS was 14,794,752
bytes, below the common 512 MiB acceptance ceiling.

The full-image destination two-sum table covered **98** point keys
and the destination-direct table **90**. Their intersection was 40;
full-only and direct-only support were 58 and 50. The changed support
is an exact finite-instance observation, not a rank or cost gain by
itself. Over three 64-target rounds, full-image collected 2,329
relations in 5,389 trials, versus destination-direct 2,628 in 6,776
trials. Source and transported-image bases had identical 2,433
relations in 5,763 trials. All relation arms reached rank nine and
recovered every scalar. The separate Rust checker verified **84/84**
process records, **1,785/1,785** scalar outputs and **13,141**
relation witnesses, using independent field, point and row formulas.

## Status and interpretation

| Requirement | Status | Evidence or remaining gap |
| --- | --- | --- |
| Different-`j` explicit pair and full rational fibers | Verified complete on `p=211`, `r=229` | [Full-fiber audit](../different_j_order_screen_20261015/FIBER_REPORT.md): 44,655/14,885 points, exact three-point fibers and 52,441 subgroup identities. |
| Exact same-target transport, relations, recovery | Verified complete on four frozen prefixes | [Public input](PUBLIC_V1.json), [raw runs](RUN_V1.jsonl), [replay](REPLAY_V1.json). |
| Complete local cold cost against source BSGS and rho | Verified complete on this L0 host | [84 child meters](METERS_V1.jsonl), [summary](SUMMARY_V1.json), [run PDF](reports/run-20261015T2114Z/REPORT.pdf). |
| Portable counted `S`, isolated performance, independent-host replay | Open | Construction/table/hash/linear-algebra weights and second-host L2/L3 run are missing. The current CPU ratios have the exact shared-host scope above. |
| Growing-order or exponent comparison | Open | This run has one subgroup order. It cannot establish a scaling crossover. |
| Direct χ-coordinate resultant | Open | Earlier [χ run](../chi_polynomial_20261015/reports/run-20261015T1851Z/REPORT.md) executes finite-scan and direct S3/quotient polynomial controls, but not this resultant. |
| Faithful curve-level Lean map proof | Partial | The [screen proof map](../different_j_order_screen_20261015/FIBER_REPORT.md) has compiled finite arithmetic; the curve/Vélu bridge remains open. |

The complete cold comparison on this instance places full-image
above both direct ECDLP controls in measured CPU, while its
destination relation support differs. This is a bounded cost and
support result. It does not identify a portable operation-count
improvement or rule out a cheaper construction of useful fiber
columns. The next discriminating experiment is to implement and
exhaustively check the direct χ-coordinate resultant against
independent pair enumeration, then run it on frozen public targets
with the same complete cost and replay gates. A separate larger-order
different-`j` panel needs a fresh bounded screen, exact map/fiber
audit, and matched direct controls before any scaling assessment.

The [validation receipt](VALIDATION.md) identifies the passing
release suite, touched examples, formatting, native replay, and the
retained failed launcher/lint attempts. Historical records are not
rewritten. The canonical ECDLP frontier and scoreboard graphs were
checked and left unchanged: the L0 single-order panel supplies no
qualified portable `S` point for them.
