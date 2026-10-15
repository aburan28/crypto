# Distinct pairing-target batch exceeds the frozen subgroup sizes

Status: `open`. Kind: proved finite counting obstruction to the specified
workload, plus an open implementation and measurement gap for a replacement
cell. Next action owner: this research task. Last updated: 2026-10-15.

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
- **Missing capability:** A separately frozen, practical pairing cell with `r>=67` and 64 distinct nonidentity public targets, plus native batch pricing and rho controls.
- **What this evidence does not establish:** It does not rule out amortized pairing benefit on a larger group or decide the cold one-target comparison.

## What would unblock it

The existing one-target costs and the possibility of setup reuse remain open.
A checkable exit is a separately frozen pairing cell with at least 65 group
elements and a practical evaluated field, with 64 distinct nonidentity
public targets, each independently replayed, plus matched 64-target rho.
Until then, run `N=1,4,16` on `r=19`, `N=1,4,16` and at most `N=60` on
`r=61`; label these as partial to the original four-count request. A new
`r>=67` cell may supply `N=64` without changing the old inputs.

- [ ] Freeze and verify an explicit pairing cell of prime order at least 67 with at least 64 distinct nonidentity targets and an evaluated pairing field.
- [ ] Recover and independently replay all 64 logarithms, charge setup plus every target and failure, and compare with 64 matched rho solves.

Meeting these conditions resolves the batch-size obstruction only. It does
not establish a single-target speedup or a scaling exponent.

## Next discriminating step

The next discriminating test is a bounded native screen of small explicit
pairing cells for `r>=67` and low embedding degree. Freeze the selector and
resource cap before the screen. A passing cell proceeds to distinct-target
batch pricing; an exhausted screen leaves this exit open and records every
candidate and reason. Neither outcome changes the cold single-target result.

Inputs and controls: the two original groups remain fixed; a new cell must
carry its own field/curve/subgroup identity, frozen seed, match to rho, and
source-answer exclusion. The bounded screen and its cap belong in a committed
protocol before execution. A passing cell proceeds to the batch; an exhausted
screen leaves this record open. Recheck on an accepted new cell or a corrected
target-count specification.

## Resolution and history

2026-10-15: identified before the reusable-input run. Unresolved.
