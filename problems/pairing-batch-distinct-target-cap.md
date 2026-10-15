# Distinct pairing-target batch exceeds the frozen subgroup sizes

Status: `resolved by a separately frozen order-67 cell`; the counting
obstruction on the two original groups remains true. Kind: finite counting
obstruction to the original per-group workload. Last updated: 2026-10-15.

## Problem

The [next-test specification](../research/isogeny_torsion_fibers_20261012/pairing_transfer_20261015/NEXT_TESTS.md)
asks for `N in {1,4,16,64}` **distinct** targets per frozen pairing group.
The two groups have prime orders `r=19` and `r=61`, so they contain only
`r-1=18` and `60` nonidentity targets respectively. Thus `N=64` cannot be
executed on either frozen group under the stated distinct nonidentity-target
contract. Reusing a target would measure a different workload. This is a
finite counting obstruction, independent of implementation or budget.

## Evidence and uncertainty

- **Observed:** The [frozen public fixtures](../research/isogeny_torsion_fibers_20261012/extension_field_20261014/FIXTURES.jsonl) and [pairing protocol](../research/isogeny_torsion_fibers_20261012/pairing_transfer_20261015/NEXT_TESTS.md) specify `r=19,61` and 64 distinct targets per group.
- **Proved:** A group of order `r` has exactly `r-1` nonidentity elements, hence at most 18 or 60 distinct nonidentity targets here. This says nothing about the cost of setup reuse for feasible batch sizes.
- **Suspected explanation:** None needed for the counting limit.
- **Implemented and checked:** The [order-67 protocol](../research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/CELL_401_PROTOCOL.md) selected `E/F_401:y²=x³+1`, `r=67`, with a degree-two pairing field and 66 distinct nonidentity public targets. The [V4 batch](../research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/RUN_401_V4.jsonl) recovered all 66 scalars in each of three rounds, with [198 independent replays](../research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/REPLAY_401_V4.json). [V6](../research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/COMPACT_PRIME_401_V6.jsonl) and the [post-test V6B repeat](../research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/COMPACT_PRIME_401_V6B.jsonl) price the 64-target prefix against both extension-field and stronger prime-field rho; their separate replay receipts each check 594 scalar outputs.
- **What this evidence does not establish:** It does not rule out amortized pairing benefit on a larger group or decide the cold one-target comparison.

## What would unblock it

The distinct-target workload gap is closed by the additional order-67 cell.
The original order-19 and order-61 groups still cannot supply 64 distinct
nonidentity targets. The new cell leaves the cold one-target and scaling
questions separate.

- [x] Freeze and verify an explicit order-67 pairing cell with 66 distinct nonidentity targets and a degree-two evaluation field.
- [x] Recover and independently replay the 64-target prefix in three rounds, charge setup, every target and failure, and compare with both extension-curve and direct prime-field rho.

Meeting these conditions resolves the batch-size obstruction only. It does
not establish a single-target speedup or a scaling exponent.

## Remaining research step

The new `N=64` workload is complete. The V5 extension-curve rho control
showed a lower CPU cost for pairing in three shared-host rounds, while the
V6 direct prime-field rho control was lower cost in all three. The V6B
repeat after the release tests exited kept the same ordering in all three
further shared-host rounds. A V7 cold direct source-group BSGS control
was lower cost than both pairing and rho at each measured prefix. Thus
the V5 pairing/extension-rho contrast is a representation diagnostic on
this cell, and the multi-target comparison must keep the direct reuse
route. Its qualified performance question requires isolated matched native
measurement, complete counted-unit conversion and growing-order cells.
Track the remaining work
in the study result report, not by reopening this finite-cardinality problem.

## Resolution and history

2026-10-15: identified before the reusable-input run. The new order-67 cell
passed field, pairing-order, bilinearity, distinct-target and independent
replay checks. V6 supplied the stronger prime-field rho reference. Resolved
for the requested `N=64` workload on an added cell; impossible on the two
original frozen groups as stated.
