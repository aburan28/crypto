# Pre-registration: does the MutantXL closure stay at degree 4 past 32 unknowns?

Registered 2026-10-16, before the closure was measured on any draw below. It
is the successor named in `research/dreg_closure_mutantxl_20261016/RESULTS.md`.

## Why

The closure study found that MutantXL's fixed-degree closure refutes all 92
committed unsatisfiable chained-`S₃` draws at degree 3–4. That held for
11–32 unknowns, ℓ = 1–5 and m = 3, 4. XL needs 5–7 on the same draws.

Held at 4 to `n = 131`, the descent cost model's best `m = 4` row falls from
`2^173.13` to `2^68.43`
(`research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md` §6). That
row assumes the degree stays at 4 from `N = 32` to `N ≈ 400`, and nothing
measures that. This study pushes `N` as far as the dense engine allows at
degree 4.

## Instrument

- **The driver.** `examples/dreg_closure` and `closure_refutation_scan` are
  unchanged in what they compute. Two changes were made, both before any
  measurement here:
  - `closure_refutation_scan` drops the plain XL basis before it builds the
    closure, so the two dense matrices are never held together. Every value
    it returns is the same.
  - `dreg_closure --to D` scans to degree `D` whatever the committed
    outcome, and `--draw K` replays one draw.
- **XL identity.** The check still applies below a committed `≥D` bound.
  Above it, XL is recorded as a new reading.
- **Fresh draws.** They come from `examples/dreg_ladder` with
  `--unsat 4 --ffd-max 3`, `d_max = 3` and seed `20261016`.
  - Each draw's solutions are counted exactly before it is kept. So every
    measured draw is unsatisfiable, and XL is known only not to refute below
    degree 4.
  - The generator output is committed under `draws/` before the closure
    runs.
- **Run order.** Every draw runs in its own process, scanning degrees 3..4
  with `--max-cols 130000`:
  `prlimit --as=15000000000 timeout 7200 dreg_closure --to 4 --max-cols 130000 --draw K FILE`.

## Cells (frozen)

`C(N, ≤4)` is the column count of the degree-4 matrix.

| m | (n, ℓ) | N | S | C(N, ≤4) | source |
|--:|---|--:|--:|--:|---|
| 3 | (11, 5) | 26 | −4 | 17,902 | committed: `dreg_surplus_control_20260925/runs/cell-11-5-6.u{0..3}` (XL ≥7) |
| 3 | (13, 5) | 28 | −2 | 24,158 | committed: `dreg_fixed_surplus_20260923/runs/cell-13-5-6.u{0..3}` (XL ≥7) |
| 3 | (15, 5) | 30 | 0 | 31,931 | committed: `dreg_fixed_surplus_20260923/runs/cell-15-5-6.u{0..3}` (XL ≥7) |
| 3 | (16, 6) | 34 | −2 | 52,956 | fresh |
| 3 | (18, 6) | 36 | 0 | 66,712 | fresh |
| 3 | (19, 7) | 40 | −2 | 102,091 | fresh |
| 3 | (21, 7) | 42 | 0 | 124,314 | fresh |
| 4 | (12, 3) | 36 | 0 | 66,712 | fresh |
| 4 | (14, 3) | 40 | +2 | 102,091 | fresh |
| 4 | (13, 4) | 42 | −3 | 124,314 | fresh |

- **Draws.** Four unsatisfiable draws per cell, 40 in total.
- **Soundness control.** The first two satisfiable draws the generator
  records in each fresh cell with `N ≤ 36` are replayed with
  `--satisfiable --to 4`.
  - If a cell records none, the report says so.
  - Neither XL nor the closure may refute a satisfiable draw. Any such
    refutation is an instrument failure, and it voids the study.

## Predictions

- **Q4: flat at 4.** In every measured draw the closure refutes at degree
  ≤ 4.
  - Confidence: moderate, about 55%.
  - Reason: in the earlier study the closure degree did not move over 11–32
    unknowns, ℓ = 1–5 and m = 3, 4. But ℓ = 6, 7 and m = 4 at ℓ = 4 are new,
    and XL rose about one degree per ℓ there.
- **Q5: XL does not refute at degree 4 in any draw.** Every committed
  XL value at these ℓ is ≥ 6, and the fits give ≥ 7. Confidence: high.
- **Q6, descriptive.** At the refuting degree, log₂(closure word ops) ÷
  log₂(columns) lies in [2.0, 2.4], as the earlier study's cells did.

## Verdict rules

Per draw:

- **closure ≤ 4:** the closure refutes at degree 3 or 4.
- **above 4:** degree 4 is built in full and saturated, and the closure has
  not refuted.
- **unmeasured:** the draw was killed by the memory or time limit. That is a
  resource limit, never evidence.

Per cell:

- **flat:** every measured draw refutes at ≤ 4, and at least 3 of the 4 are
  measured.
- **rises:** any measured draw is above 4.
- **unmeasured:** otherwise.

Study:

- **"flat to N = X":** every cell is flat or unmeasured. `X` is the largest
  flat `N`.
- **"rises at N = Y":** some cell rises. `Y` is the smallest `N` where one
  does. Every cell is reported either way.

## Cost accounting, descriptive

- Recorded per draw: columns, closure products, counted word operations,
  wall seconds and peak RSS.
- No wall-clock claim. The container is the 4-core, 16 GB, non-isolated one
  (`AGENTS.md` §10).

## Stop conditions

- Each draw has a 2-hour limit and a 15 GB address-space limit. A kill
  marks that draw unmeasured.
- An XL-identity failure on a committed draw aborts that draw, and it is
  reported as an instrument failure.

## Scope and class

- **Class: stage diagnostic.** No `S`, no rho ratio, no speed claim.
- **Systems.** The same chained `S₃` systems as before: `b = 1`, random
  subspaces, `n ≤ 21`. Field subfields are as the generator draws them
  (`AGENTS.md` §8b).
- **Against ECC2K-130.** The descent note's verdict is degree-independent,
  and a flat result cannot reverse it: the free-oracle floor leaves
  `2^5.71` per target.
  - A flat result extends the measured range of the degree-4 row's
    assumption from `N ≤ 32` to `N ≤ 42`, still about 10× short of the
    priced `N`.
  - A rise bounds that row.
- **Conductor.** No Conductor task was taken: the CLI is not installed in
  this container, and its hooks fail open.
