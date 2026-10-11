# Inherited-F4 remap round: protocol

Status: **preregistered**, written before any confirmatory run. This is
the second engineering round on the inherited-F4 engine, after
[`../inherited_f4_overhead_20261005`](../inherited_f4_overhead_20261005/RESULT.md)
(#1386). The exploratory steps that motivated it are disclosed below with
their numbers.

## Question

After round 1, the largest single loop in a whole `ecbench exec` of the
inherited-F4 arm is `rewrite`'s column remap. It moves a row into a new
layout one set bit at a time, through a lookup in the composed column map,
including the bits the map deletes. On `L19` seq 25 it was about 13% of
the process's instructions. Can it move a 64-column word at a time where
the map allows it, without changing a single decision, count or answer?

## What changes, and what must not

Only `src/cryptanalysis/inherited_f4.rs` changes:

1. **`ColMap`.** Every layout map now carries, for each 64-column word of
   its source layout, `live` (the columns it keeps) and `out` (the first
   target when the kept columns land on consecutive targets in order).
   - A deletion step (`v := 0`) builds both in the loop that lays out the
     child. It is consecutive on every word.
   - A fold (`v := 1`) and a layout extension delete nothing, so they get
     all-ones masks and are marked scattered.
   - A composed map takes its live masks in the same pass as the gather
     and is marked scattered.
2. **`remap_row`.** It skips the deleted bits with one AND per word.
   - On a CPU with BMI2 (detected at run time), a consecutive word moves
     with one `pext` plus a shift.
   - Every other word goes bit by bit, as before.
   - Both paths write the same bits and report the same lowest word, so
     the specialisation charge, which counts words from that one, is
     unchanged.

Nothing algebraic changes: same rows, same order, same layout, same
pivots, same word-XOR charges. The gate below enforces this run by run.

## Boundaries (AGENTS.md §1)

These are unchanged from round 1.

- **Floor:** `√(π/4m)`.
- **Reference:** same-target strong rho, from the
  [F6-IC ladder](../f6_ic_ecbench_ladder_20261005/RESULT.md).
- **Class:** in the ecbench unit this round cannot move `S`, so it is
  **engineering** at best.

## Frozen inputs

These are the round-1 suite plus two fresh holdout sets whose seeds were
fixed before any run of this round.

| set | spec | sequences | role |
| --- | --- | --- | --- |
| L13 | `../f6_ic_ecbench_ladder_20261005/SPEC-n13.json` | round 1, all 16 IC runs | ladder, m = 13 |
| L19 | `../f6_ic_ecbench_ladder_20261005/SPEC-n19.json` | round 1, all 16 IC runs | ladder, m = 19 |
| L23 | `../f6_ic_ecbench_ladder_20261005/SPEC-n23.json` | seq 24, 25, 27, 28 | ladder, m = 23 |
| H13 | `../inherited_f4_overhead_20261005/SPEC-holdout-n13.json` | all 8 | round-1 holdout, m = 13 |
| H19 | `../inherited_f4_overhead_20261005/SPEC-holdout-n19.json` | all 4 | round-1 holdout, m = 19 |
| H13b | `SPEC-holdout-n13b.json` (target seed 202610071013) | all 8 | fresh holdout, m = 13 |
| H19b | `SPEC-holdout-n19b.json` (target seed 202610071019) | all 4 | fresh holdout, m = 19 |

Exploration in this round used only L13 seq 25, 28 and 31 and L19 seq 25.

## Binaries

- **Baseline:** `ecbench` from `main` at
  `a3652eddafb949073bf0f96b210a75b671cfb346`, which includes round 1.
- **Candidate:** the same, plus this branch's `inherited_f4.rs`.
- Both are release builds of the default profile, built on the measuring
  host with the same `rustc`. Their digests are recorded with the result.

## Measurements

Each runs with `run.sh` and is the same as round 1:

1. **Gate.** Byte-identical stripped child output on every input in every
   set (`strip.jq`).
2. **Primary: Callgrind total `Ir`, whole process.** One run per binary and
   input, on L13, H13, H13b, L19 seq 24, 25, 27 and 28, H19 and H19b.
3. **Secondary: pinned CPU time.** Five ABAB rounds with a byte-identical
   A/A copy, on L19 seq 25 (F4), L19 seq 24 (F6) and L23 seq 25 (F4).

## Decision rule

- **Abort** if any input fails the gate.
- **Success (engineering)** if the gate passes on all inputs and both of
  these hold:
  - the pooled candidate/baseline `Ir` over the confirmatory m = 19 inputs
    (L19 seq 24, 27 and 28, H19 and H19b) is at most **0.97**;
  - no `Ir` input regresses by more than **2%**.
- **Null** if the gate passes but either test fails. That is reported, and
  nothing is merged as a speedup.
- Wall time never decides.

The 0.97 is set knowing the exploratory figure on L19 seq 25, which the
last exploratory step put at 0.936 against a baseline built from round 1's
head. It is preregistered against the other inputs.

## Host-class scope, stated in advance

The `pext` path runs only where the CPU reports BMI2, and this host has it.
On x86-64 without BMI2 and on Arm64, the masked bit-by-bit path runs; a
unit test pins it to the same output, but its speed is not measured here.
Some AMD cores before Zen 3 implement `pext` in microcode, slowly. They
are not measured either, and the result does not claim a gain on them.

## Exploratory record (before this protocol)

All on L19 seq 25, Callgrind total `Ir`, with every candidate's output
identical to the baseline's. The baseline here is round 1's head
(`5ea04af9`), at 52,095,031,593.

| step | change | `Ir` | vs round-2 baseline |
| --- | --- | --: | --: |
| r2c1 | `ColMap` masks read off every map; `pext` remap | 50,434,351,564 | 0.968 |
| r2c2 | masks built by the deletion step itself; folds and extensions all-live | 49,334,070,255 | 0.947 |
| r2c3 | branch-free generic mask pass | 50,044,613,069 | 0.961 (rejected) |
| r2c4 | masks fused into the composition's gather | 49,430,580,065 | 0.949 |
| r2c5 | composed maps: live mask only, marked scattered | 48,779,426,384 | **0.936** |
| r2c6 | drafts held as unique `Rc` slices, reduced in place | 49,198,373,797 | 0.944 (rejected) |

Rejected:

- **r2c3:** conditional moves cost more instructions than the branches
  they replaced.
- **r2c6:** `Rc::get_mut` in the reduction loop cost more than the
  allocation it saved.
- **Consecutive testing of composed words (r2c2/r2c4):** it cost more than
  `pext` saved on them.

## Commands

```sh
OUT=$PWD/results/run-1 ./run.sh build
OUT=$PWD/results/run-1 ./run.sh inputs
OUT=$PWD/results/run-1 ./run.sh gate
OUT=$PWD/results/run-1 ./run.sh ir 3
OUT=$PWD/results/run-1 ./run.sh wall 5     # alone on the host
```
