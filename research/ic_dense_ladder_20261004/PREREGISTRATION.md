# The `ℓ = 6` rung with a dense in-tree engine: does the direct `S₄` descent reach degree 11?

Registered before any registered cell runs. §2 lists everything that ran before
registration and what was seen.

## 1. Why

- **The rung the last round could not reach.**
  [ic_gb_ladder_20261003](../ic_gb_ladder_20261003/RESULTS.md) confirmed with an
  external engine that the refutation degree of the direct `S₄` descent (`x4`) reads 5, 8,
  9, 10 at `ℓ = 2…5` on `K₁/2¹⁷` and the symmetric norm form (`rr`) 4, 5, 6, 8, and left
  `ℓ = 6` censored: Singular's truncated `slimgb` ran out of 12 GB at degree 8 (`rr`) and 9
  (`x4`) on both `n = 19` curves. Its verdict was *inconclusive at `ℓ ≥ 6`*, and the
  registration named the next step an engine question. This round is that engine.
- **The engine.** `examples/macaulay_dense.rs` builds the affine Boolean Macaulay matrix
  at one degree from the same exported systems (every product `m·f`, `deg m ≤ D − deg f`,
  reduced modulo `x_i² = x_i`, over the multilinear monomials of degree `≤ D`) and reduces
  its rows, generated in batches, against a bit-packed echelon basis with the lowest set
  column as pivot and the constant monomial as the last column, so that a row reducing to
  that single bit is the polynomial `1`. It is the textbook Macaulay/XL observable, the same
  one the in-tree sparse scan and the external engine measure, from a third piece of code;
  it has no floor at the system degree and no Gröbner machinery. The basis is at most
  `n_cols × n_cols` bits: 3.6 GB at `rr ℓ = 6` degree 8 (169k columns), 5.0 GB at `x4 ℓ = 6`
  degree 10 (199k) and 6.7 GB at degree 11 (231k), inside this host's 15 GB where Singular's
  representation was not.
- **What is being asked.** At `ℓ = 6` on the two `n = 19` curves: does `x4` refute at 11
  (the law `ℓ + 5`, four above the sharp first-fall bound 7) and `rr` at 8 (`ℓ + 2`)? A
  reading of 11 would be the first measured excess of 4 over the Kousidis–Wiemers bound,
  on the canonical object, still rising; a reading at or below 10 for `x4` on both curves,
  or at or below 8 for `rr` with 7 on the other curve, would be the first flattening.
- **The prior is growth,** as before; the asymmetric payoff is unchanged.

## 2. What ran before registration (disclosed; not evidence)

On this container, 2026-10-04, after the engine was written:

- **Calibration, draw 0 of `K₁/2¹⁷ ℓ = 2, 3, 4`, `rr` and `x4`, degrees `D − 1` and `D`
  around every known reading.** Every reading equals the external engine's and (above its
  floor) the in-tree scan's: `rr` 4, 5, 6; `x4` 5, 8, 9; the trace-constant draw refuted at
  1. All in under 50 ms.
- **Reach, draw 0 of `K₁/2¹⁷ ℓ = 5`:** `rr` refuted at 8 in 31 s (39,203 columns,
  111,911 rows, rank 29,931 at refutation) and `x4` at 10 in 24 s (30,827 columns, 32,997
  rows), four threads, well under 1 GB; the external engine needed 60 s and 430 s with more
  than 4.5 GB for the same two readings.
- **Plumbing** (`SMOKE`): the `K₁/2¹⁷ ℓ = 2` cell through `run.sh` into
  `gb_ladder_analyze`, readings as above.

Nothing above is cited as a result; the registered run re-runs the calibration cells.

## 3. Object

- **Systems.** The exported systems of the earlier round, byte-identical (their sha256 is
  recorded in [ic_gb_ladder_20261003/runs/registered/dump.sha256](../ic_gb_ladder_20261003/runs/registered/dump.sha256);
  the files were regenerated from the recorded binary and seed and re-verified against it).
  Same draws, same rootless rule.
- **Observable.** The least `D` at which the dense Macaulay matrix contains `1`, scanned
  upward one process per `D` from 1 to a cap of 9 (`rr`), 11
  (`x4`) and 7 (`ctrl`). Satisfiable controls cannot be resolved by this engine (it does
  not detect pinning); they are read as "not refuted by 7", which is all P4 uses.
- **Cells.** `K₁/2¹⁷ ℓ = 2, 3, 4, 5` (calibration against both earlier engines), then
  `K₁/2¹⁹ ℓ = 6` and `K₀/2¹⁹ ℓ = 6`, `rr` then `x4` then controls, decisive readings first.
- **Budget and censoring.** One process at a time, four threads, 86,400 CPU-s (`ulimit -t`,
  about six hours wall) and 13.5 GB of address space per process, with the engine refusing
  before allocation any basis above 12 GB. A killed or refused process censors the draw at
  `≥ D`; censoring is never negative evidence. No retries.
- **Isolation.** Degrees are not timings; wall seconds are advisory.

## 4. Analysis

`examples/gb_ladder_analyze.rs`, unchanged: readings, medians (same rule, `triv`
excluded), slopes, pairs, excess over 7, and calibration against the in-tree ladder; the
external engine's readings are compared by hand from its readout.

## 5. Predictions (pass/fail)

- **P1, calibration.** On every draw of the four `K₁/2¹⁷` cells, this engine's reading
  equals the external engine's (`rr` 4, 5, 6, 8; `x4` 5, 8, 9, 10; the two trivial draws).
- **P2, `x4` at `ℓ = 6`.** The median is exact and equals 11 on each `n = 19` curve
  (excess 4).
- **P3, `rr` at `ℓ = 6`.** The median is exact and equals 8 on each `n = 19` curve.
- **P4, controls.** No control is refuted at or below 7 at any `ℓ`.
- **P5, pairs.** `rr` is below `x4` on every draw exact on both arms, by 2 to 4.

## 6. Decision rule, fixed now

- **Growing**: P2 and P3 hold on at least one curve each and no reading on the other is
  below them. The law holds one rung further, the first-fall assumption is four degrees
  behind the data on the canonical object, and the stop decision's item 3 stays closed.
- **Saturation signal**: an exact `x4` median `≤ 10` at `ℓ = 6` on both curves, or an exact
  `rr` median `≤ 8` on one curve with `7` on the other, with P4 holding. Then the algebraic
  route's exponent audit is reopened as a Coordinator decision, and the next round is the
  same cells with more draws before anything else.
- **Inconclusive**: `x4` at `ℓ = 6` resolves on neither curve, or `rr` on neither. Then
  the reach is reported with its censoring; the step after is `ℓ = 7`, which needs a
  larger host (1.2 M columns at degree 9).

## 7. Scope

Two Koblitz curves, `n` 17–19, `ℓ` 2–6, `m = 3`, four rootless draws per cell, one
in-tree engine. A stage diagnostic (AGENTS.md §2, §5): no end-to-end cost, nothing at
`n ≈ 83` or 131, nothing asymptotic; at `m = 3` nothing that could tie rho even with a free
oracle.
