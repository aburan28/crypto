# Four-order complete-solve panel: counted arithmetic and relation support

**Status:** the V2 prime-to-characteristic panel completed on 2026-10-15.
Its [frozen target/run spec](BATCH_RUN_SPEC_V1.json), [source and binary
hashes](BATCH_SOURCE_FREEZE_V1.json), [three full-fiber audit
report](reports/run-20261015T2335Z-fibers/REPORT.md), and [descriptive
run report](reports/run-20261015T2341Z-batch/REPORT.md) separate the
construction, complete solve, and analysis gates. The new solver processed
252/252 cold child processes and returned 5,355/5,355 correct scalar
outputs. Independent Rust formulas replayed 37,061 relation witnesses.
Together with the checked order-229 reference, the
[native trend derivation](TREND_V1.json) covers 336 children and 7,140
scalar outputs across four orders. PR review, isolated wall timing and
second-host performance replay are separate gates.

## Boundary and one counted unit

The generic negation-orbit collision floor is `sqrt(pi/4)=0.886226925`
per target. Direct source negation rho on each identical public target
prefix is the applicable generic reference; reusable direct source BSGS
is the memory control. `S_group = (all additions + all doublings) /
(N sqrt(r))`, cold from common checks through setup, 64 unknown-scalar
solves and final verification, follows the repository's group-operation
count. It remains a **lower bound on complete portable `S`**: field work
outside group operations, square roots, hash/table and row operations,
allocation, startup and I/O lack a pinned common conversion. The separate
Fp multiplication/inversion axes in `TREND_V1.json` overlap counted
group work and must not be added to this table.

The table uses the *same unit throughout*: median `S_group` per target for
the 64-target cold child, three rotated rounds per cell. `r=229` is the
independently checked V1 reference; the other cells have `r != p` and
different source/destination `j` values.

| Arm | `p=211,r=229` | `p=457,r=421` | `p=409,r=433` | `p=421,r=457` |
| --- | ---: | ---: | ---: | ---: |
| Source projected relations | 58.28 | 80.04 | 78.84 | 64.15 |
| Isogeny-image relations | 64.56 | 86.17 | 85.39 | 69.73 |
| Destination-direct relations | 69.56 | 81.83 | 68.04 | 72.28 |
| Full rational-image relations | 148.49 | 405.89 | 332.16 | 336.69 |
| Reusable source BSGS | 1.76 | 1.60 | 1.52 | 1.54 |
| Source negation rho | 6.64 | 5.59 | 5.14 | 5.31 |
| Rho A/A repeat | 6.64 | 5.59 | 5.14 | 5.31 |

![Four-order counted cost trend](TREND_V1.svg)

The paired median **full-image/rho `S_group` ratios** are 22.355, 72.401,
64.275 and 63.389 in the table's order. Paired full-image/BSGS ratios
are 84.147, 253.690, 218.661 and 218.524. In all four cells, every
prefix and every round, the full-image arm exceeded both direct controls
in each recorded arithmetic axis: group additions/doublings, Fp
multiplications, and Fp inversions. This is a componentwise result for
the fixed implementation, regardless of positive weights assigned to
those three counted axes. It does **not** establish a complete portable
cost ordering when native categories are still unpriced. The largest
full-image N64 peak RSS in the new cells was 35,110,912 bytes, well below
the common 512 MiB ceiling.

## Relation geometry and complete recovery

The full-image two-sum table changed support relative to destination-direct
at every cell, but it did not produce a consistent support increase:

| `r` | full/direct support keys | full relations/trials | direct relations/trials |
| ---: | ---: | ---: | ---: |
| 229 | 98 / 90 | 2,329 / 5,389 | 2,628 / 6,776 |
| 421 | 108 / 112 | 2,457 / 9,475 | 2,348 / 8,853 |
| 433 | 106 / 120 | 2,393 / 9,595 | 2,234 / 8,113 |
| 457 | 112 / 122 | 2,376 / 9,556 | 2,339 / 8,830 |

All four relation arms reached rank nine in the complete children,
recovered each requested scalar and passed source/image replay. The
separate checker rebuilt every signed base, map, pair-table witness,
relation row and scalar from its own field/curve formulas. The three
new replay receipts are
[`409/433`](BATCH_409_433_V1.replay.json),
[`421/457`](BATCH_421_457_V1.replay.json), and
[`457/421`](BATCH_457_421_V1.replay.json). Their [run
records](BATCH_409_433_V1.run.jsonl),
[`421/457`](BATCH_421_457_V1.run.jsonl), and
[`457/421`](BATCH_457_421_V1.run.jsonl) retain every capped attempt,
count and output. All 252 child exits were zero; none required a replacement
target or cap change. The initial target-freeze shell error and its
corrected continuation remain in [TARGET_FREEZE_LAUNCHER_V1.md](TARGET_FREEZE_LAUNCHER_V1.md).

## Interpretation and next discriminating work

The observed `S_group` rises from 148.49 at order 229 to 332–406 in the
larger cells, while direct rho remains about 5–7 and BSGS about 1.5–1.8.
The full rational-image base is rebuilt by enumerating `#E(F_(p²))`
source points in every cold child; its setup therefore grows with the
field and dominates this implementation. The four orders span only a
factor of two and the curves have different cofactors, so a fitted
exponent would be unstable and is deliberately unset. A broader screen
with larger prime-to-characteristic orders and identical memory/cost
controls would be needed for scaling inference.

The most specific structural lead is to construct useful destination
columns from a *direct chi-coordinate resultant/root oracle* without
enumerating the entire rational source. Its exact decomposition set must
first match independent pair enumeration, and then its whole cold cost
must be compared with direct source rho and BSGS. The pairing-transfer
track remains separate: compare complete torsion-input search, Miller,
final exponentiation and finite-field DLP cost with direct ECDLP for an
exact target count. Isogeny transport, relation support, finite-field
transfer and a practical advantage remain distinct obligations.

For an ECDLP frontier claim, the missing gates are a fully priced portable
unit with every native category and memory, an L2-or-better interleaved
wall session with A/A admitted, and a different-host performance replay
of the exact frozen inputs and binary. The current child CPU timings are
shared-host L0 descriptions only. No canonical frontier or scoreboard
graph is changed by this bounded result; this report and its visuals carry
the four-order evidence. The direct chi resultant and the curve-level Lean
Vélu bridge remain open work, not conclusions from this panel.

An [isolated-host access check](ISOLATION_ACCESS_V1.md) found no configured
isolab and a zero Standard On-Demand EC2 vCPU quota in five checked regions;
the requested bare-metal session therefore has no available placement from
this environment. A [frozen independent Linux replay
protocol](INDEPENDENT_HOST_PROTOCOL_V1.md) and CI workflow are prepared for
the exact public inputs. Its hosted VM is a different host class, but its
timing remains descriptive and does not satisfy the L2 gate.
