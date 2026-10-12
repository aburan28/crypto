# The external-engine `m = 3` ladder: results

Written after the run. [PREREGISTRATION.md](PREREGISTRATION.md) is unchanged since its
registration commit `5f7f088ee`, apart from its five dated, additive amendments (the retry
phase's CPU budget; one memory death per cell and arm; the controls' retries, their degree
cap, and none at `ℓ ≥ 7`). The readout
is [runs/registered/readout.txt](runs/registered/readout.txt), printed by
`examples/gb_ladder_analyze`.

## Run

| | |
|:--|:--|
| engine | Singular 4.3.2 (Ubuntu `1:4.3.2-p10+ds-1.1build1`), `slimgb` on the homogenisation with `degBound = D`, driver [refute.sing](refute.sing) |
| systems | exported by `examples/rr_degree_ladder.rs --dump-dir` from the unchanged builders, seed 20260930, four rootless draws per cell; `runs/registered/dump/` holds every system and the draw catalogue |
| phase 1 | 2026-10-03 21:57Z – 2026-10-04 00:07Z, three lanes on CPUs 1–3, 3,600 CPU-s and 4.5 GB per process |
| phase 2–3 | 2026-10-04 00:08Z – 09:27Z, one process at a time at 12 GB (Amendments 1–5) |
| engine record | [runs/registered/engine.txt](runs/registered/engine.txt): Singular version, sha256 of the driver and of the two binaries |
| systems | every exported system is hashed in [runs/registered/dump.sha256](runs/registered/dump.sha256); the `.sing` files themselves (280 MB) are not committed and regenerate byte-identically from `rr_degree_ladder --dump-dir` at the recorded binary |
| cells | nine registered cells, all run; no cell skipped |

## Registered verdict: inconclusive at `ℓ ≥ 6`, growth confirmed at `ℓ = 5`

By §6's rule the verdict is **inconclusive**: one new rung resolved for each arm (`ℓ = 5`),
not the two (`x4`) and three (`rr`) the "growing" reading requires. The rung that would
have separated the readings, `ℓ = 6`, is beyond this engine on this host: both `n = 19`
curves die at the 12 GB limit at degree 8 for `rr` (after 44 and 55 minutes) and at degree 9
for `x4`. What did resolve is exact, on every draw, and on the side of growth:

| `ℓ` | `rr` (3ℓ + 1 unknowns) | `x4` (3ℓ unknowns) | `x4 − 7` | control |
|--:|:--|:--|:--|:--|
| 2 | 4 4 4 (one `triv`) | 5 5 5 6 | −2 | pinned at 5 (three), ≥13 |
| 3 | 5 5 5 5 | 8 8 8 8 | 1 | pinned at 6 (two), ≥13 |
| 4 | 6 6 6 6 | 9 9 9 9 | 2 | pinned at 7 (one), ≥13 |
| 5 | 8 8 8 (one `triv`) | **10 10 10 10** | **3** | ≥ 8 ≥ 8 ≥ 8, ≥ 13 |
| 6 (`n = 19`, both curves) | ≥ 8 (one draw each), ≥ 7 | ≥ 9 | ≥ 2 | ≥ 7 (one each), ≥ 6 |
| 7 (`n = 23`) | ≥ 7 | ≥ 9 (one draw), ≥ 8 | ≥ 1 | ≥ 6 (not retried) |
| 8 (`n = 29, 31`) | ≥ 7 (one draw each), ≥ 6 | ≥ 8 | ≥ 1 | ≥ 6 (not retried) |

Slopes of the exact medians on `K₁/2¹⁷`: `rr` 1.30 over `ℓ = 2…5`, `x4` 1.60 over
`ℓ = 2…5`. Paired on the 14 draws exact on both arms, `rr − x4` is −2.43 (range −3 to −1).

## Predictions

1. **P1, calibration: holds.** On the 20 shared draws where both readings are exact and the
   in-tree reading is above its floor, the two engines agree on every one. `x4` at `ℓ = 2`
   reads 5 on three of four draws (the in-tree floor of 6 hid it).
2. **P2, growth of `x4`: holds at `ℓ = 5`, unresolved at `ℓ = 6`.** 10 on all four draws;
   excess 3 over the sharp first-fall bound. At `ℓ = 6` the one draw retried at 12 GB on
   each curve dies at degree 9 after 30 minutes: `≥ 9`, consistent with the law's 11 and
   with nothing below 9.
3. **P3, growth of `rr`: holds at `ℓ = 5`, unresolved at `ℓ = 6`.** 8 on every genuine
   draw; at `ℓ = 6`, `≥ 8` on the one retried draw of each curve, consistent with 8 or 9
   and with nothing lower.
4. **P4, control: holds.** No control is refuted at all. At `ℓ = 2, 3, 4` the controls
   resolve by pinning (5, 6, 7) or scan to the cap unresolved; at `ℓ = 5` all four pass
   degree 7 unrefuted (`≥ 8`, one scanned to the cap with 15 of 16 variables pinned),
   so none is below `rr`'s 8. At `ℓ ≥ 6` the controls are censored at or above `rr`'s
   own bound.
5. **P5, pairs: holds.** `rr` is below `x4` on all 14 paired draws, by 1 to 3.

## What the engine reached, and what it cost

Advisory wall seconds per process (Singular's own clock agrees within 3%):

| cell, arm, degree | wall | memory |
|:--|--:|:--|
| `K₁/2¹⁷ ℓ = 4`, `rr` 6 / `x4` 9 | 0.4 / 9 | < 1 GB |
| `K₁/2¹⁷ ℓ = 5`, `rr` 8 | 60 | < 1 GB |
| `K₁/2¹⁷ ℓ = 5`, `x4` 9 (unrefuted) / 10 (refuted) | 150 / 430 | > 4.5 GB at 10, < 12 GB |
| `n = 19 ℓ = 6`, `rr` 7 (unrefuted) / 8 | 500 / dies | > 4.5 GB at 7, > 12 GB at 8 |
| `n = 19 ℓ = 6`, `x4` 9 | dies in 50 s at 4.5 GB, in 30 min at 12 GB | > 12 GB |
| `n = 23 ℓ = 7`, `rr` 7 / `x4` 8 (unrefuted) / `x4` 9 | dies in 120 s / 310 / dies in 305 s | > 12 GB / < 12 GB / > 12 GB |
| `n = 29, 31 ℓ = 8`, `rr` 6 (unrefuted) / 7 | 370–500 / dies | < 12 GB / > 12 GB |

Against the in-tree sparse Macaulay scan, which did not resolve `ℓ = 5` `rr` in 6,000 s and
never scanned `x4` past 9, this engine resolved `ℓ = 5` on both arms in minutes. It is
memory that stops it one rung later: the truncated basis at degree 8 in 19 unknowns, and at
degree 9–10 in 18, does not fit in 12 GB. Two trials of Singular's `std` in place of `slimgb`
during the smoke runs were 10–30× slower; neither was tried for memory and no other engine
was registered.

## What this says

- **Measured, on the canonical object.** The Weil descent of `S₄` with `x_i ∈ V`
  refutes at degree 5, 8, 9, 10 for `ℓ = 2, 3, 4, 5` at `n = 17`: three degrees above the
  sharp first-fall bound 7 at `ℓ = 5` and still rising one per rung from `ℓ = 3`. The
  first-fall-degree assumption, `D_reg = D_ff + o(1)`, is behind these readings by 3 at the
  largest size resolved, with an independent engine confirming every in-tree reading above
  its floor. This is the strongest statement the round can make; it is a statement at
  `n = 17`, `ℓ ≤ 5`, and nothing asymptotic.
- **The symmetric norm form tracks it two to three degrees lower** and rises at the same
  rate (8 at `ℓ = 5`; `≥ 8` at `ℓ = 6`), as the earlier ladder read: a constant, not a
  slope.
- **No saturation signal.** Nothing read equal or lower across a rung; the `ℓ = 6` bounds
  are consistent with the law and inconsistent with any reading below 8 (`rr`) or 9 (`x4`).
- **The next rung is an engine question.** By §6, the step after an inconclusive verdict is
  reach, not mathematics: degree 8 in 19 Boolean unknowns needs more than 12 GB in this
  engine's truncated-basis representation. A dense Macaulay elimination at that size is
  about 650k × 262k bits, 21 GB, so a 32–64 GB host would resolve `ℓ = 6` for `rr` with
  either method; `x4` at `ℓ = 6` (degree 11 in 18 unknowns) is larger still. Nothing here
  reopens the stop decision's item 3; the route's `m ≥ 4` exponent stays as the audits read
  it.

## Erratum carried

The in-tree ladder's floor and its two trace-constant draws are recorded as Erratum 1 of
[ic_rr_norm_ladder_20260930/RESULTS.md](../ic_rr_norm_ladder_20260930/RESULTS.md); this
round's `triv` readings are the same draws.

## Scope

Two Koblitz curves, `n` 17–31, `ℓ` 2–8, `m = 3`, four rootless draws per cell, one external
engine. The observable is the Macaulay refutation degree of one presentation; a stage
diagnostic (AGENTS.md §2, §5): no end-to-end cost, no yield, nothing at `n ≈ 83` or 131,
nothing asymptotic, and at `m = 3` nothing that could tie rho even with a free oracle.
