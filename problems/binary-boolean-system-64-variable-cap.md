# Binary decomposition systems are capped at 64 Boolean unknowns

- **Status:** `open`
- **Kind:** implementation gap
- **Affected work:** [`research/ic_m3_algebraic_vs_mitm_20261014/`](../research/ic_m3_algebraic_vs_mitm_20261014/README.md)
  (the `m = 4` algebraic arms at n ≥ 23); any Gröbner/SAT decomposition of a
  binary summation-polynomial system with more than 64 Weil-descended unknowns
- **Next action owner:** unassigned
- **Last updated:** 2026-10-14

## Problem

We cannot measure the repository's Gröbner arms (`groebner_decompose`,
`groebner_decompose_f6_ic`) or its SAT arm (`sat_decompose_with`) on a binary
decomposition system with more than 64 Boolean unknowns. All three build the
system with `koblitz_groebner::build_decomposition_system`, which refuses
`n_vars > MAX_VARS = 64` (`src/cryptanalysis/koblitz_groebner.rs:140`). The
`pdp_split` systems pack assignments into a `u64` and have the same limit.

The S₃-chain system for `m` summands over a dimension-`d` subspace of
`F_{2^n}` has `md + (m − 2)n` unknowns. The cap therefore removes:
- `m = 4` at `d = ⌈n/4⌉` from n = 23 on (70 unknowns at n = 23);
- `m = 3` beyond n = 31 at `d = ⌈n/3⌉`.

The study wanted an `m = 4` slope over n = 13…31. It got four cells
(n = 11–19) instead of seven.

## Evidence and uncertainty

- **Observed:** `ic_algebraic_floor solve --m 4 --arm f4-chain --cell
  1:23:6` returns `unsupported` immediately. The same holds for `f6-ic`, and
  for n = 29 and 31 (`research/ic_m3_algebraic_vs_mitm_20261014/runs/m4/`,
  `runs-m4.log`). The floor ran on all seven cells.
- **Derived or proved:** the counts follow from the system's construction
  (`build_decomposition_system`). Nothing about cost is derived.
- **Suspected explanation:** monomials are `u64` bitmasks throughout the
  Boolean polynomial and F4 code. This is a representation choice, not a
  mathematical limit.
- **Missing capability:** a Boolean polynomial and monomial type for more
  than 64 variables (`u128`, or a `[u64; k]` bitset) in the chain builder,
  inherited F4 and the SAT encoder.
- **What this evidence does not establish:** what the `m = 4` algebraic
  arms cost at n ≥ 23. Extrapolating the n ≤ 19 slope is not a measurement.

## What would unblock it

- [ ] `build_decomposition_system` and `matrix_f4_f2` / inherited F4 accept at
  least 128 unknowns, and pass the existing decomposition tests. They also
  agree with the floor's ground truth on 12 targets at `m = 4`, n = 23,
  d = 6 (`ic_algebraic_floor run --m 4 --cells 1:23:6 --arms f4-chain`).
- [ ] The same holds for `sat_decompose_with` (chain encoding).

Meeting these lets the frozen `m = 4` grid run to n = 31. It does not by
itself say anything about cost.

## Next discriminating step

- **Action:** widen the monomial type behind `MAX_VARS` (a 128-bit monomial
  is enough for the `m = 4`, n ≤ 31 cells: at most 94 unknowns).
- **Inputs and controls:** the study's frozen `m = 4` cells and seed. Use
  the floor as ground truth, and the n ≤ 19 cells as a regression check
  (verdicts and word XORs unchanged).
- **Execution limits:** the study protocol's 120 s and 4 GiB per target.
- **Outcomes:**
  - The arms finish: refit the `m = 4` slope on seven cells.
  - They censor: the cells become lower bounds under the same rule.
- **Recheck trigger:** a change to `koblitz_groebner::MAX_VARS` or to the
  Boolean polynomial type.

## Resolution and history

- 2026-10-14: recorded from the `m = 4` grid of
  `research/ic_m3_algebraic_vs_mitm_20261014/`. Unresolved.
