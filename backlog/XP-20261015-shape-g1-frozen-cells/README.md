# XP-20261015-shape-g1-frozen-cells: goal G1 on its frozen cells

**Status: blocked** on PR #1892 (open; the card is ready only on a merged
commit) and on memory. Funding 250,000.

## Question

Goal G1 was frozen (commit `110117d8c`) to ask whether, on class-0 targets at
`m = 3`, the `gp` and `fp` factor-base shapes have a mean refutation degree at
least 0.5 below random at two of `(n, l) = (13, 5), (17, 6), (19, 7)`.

## Why it matters

The G1–G6 study (PR #1892, `research/shape_search_g1_g6_20261015/`) killed G1
on its *scaled* cells (`n = 5, 7, 11`: gap 0.000) and recorded the frozen
cells as "not testable (resource)". So G1's verdict on its own terms is open,
and the prior here is a null.

- A gap ≥ 0.5 at two frozen cells revives the shape line and opens G6
  (`n = 31` with the real F4 engine, then the `m = 83` gate of AGENTS.md §8a).
- No gap closes G1 on its frozen cells and retires the shape builders.
- A censored result (both arms at the degree bound) closes nothing, and the
  report says so.

## Prerequisites

- PR #1892 merged; pin its merge commit (the card pins the branch head
  `f12cb92` until then).
- A host per cell: degree-7 matrices need 47–50 GB at `(13, 5)` and more at
  `(17, 6)` and `(19, 7)` (35 and 40 unknowns; `C(35, 6) ≈ 1.6·10⁶`,
  `C(40, 6) ≈ 3.8·10⁶` degree-6 columns). Size each from a timed probe.
- `PROTOCOL.md` here, frozen before ready: instances and targets per cell,
  the dense budget, the probe, and the one-sided rule for `(13, 5)` that the
  study's protocol used.

## Expected cost

Not estimable from committed data. One random draw at `(13, 5)`, `d_max = 6`,
did not finish in 60 minutes under a 2.5 GB cap (study README, "The frozen
cells"). The probe in `PROTOCOL.md` produces the estimate.

## What "done" means

All three cells scored by `factor_base_shape_ladder score`, with
`selftest` passing, landed in `research/shape_search_g1_frozen_cells_<date>/`
by a PR that updates the shape-search note's G1 row. Execution: three cairn
lab tasks; see [backlog/README.md](../README.md).
