# Pre-registration: what sets the MutantXL closure degree — ℓ, N or m·ℓ?

Registered 2026-10-16, before the closure was measured on any draw below. It
is the successor named in `research/dreg_closure_scale_20261016/RESULTS.md`.

## Why

The scale study measured the closure at degree 4 on 40 draws.

| | outcome at degree 4 |
|---|---|
| m = 3, ℓ = 5, 6 (N ≤ 36) | refutes on every draw |
| m = 3, ℓ = 7 (N = 40, 42) | does not refute on any draw |
| m = 4, ℓ = 3, 4 (N up to 42) | refutes on every draw |

Three explanations fit those cells. Each one decides how the degree-4 and
degree-5 rows of the ECC2K-130 cost model extrapolate
(`RESEARCH_ECC2K130_DESCENT_DEGREE.md` §7).

- **N:** the degree rises once the unknown count passes about 38, whatever
  `ℓ` is.
  - This is already in tension with m = 4 holding at N = 40 and 42. A
    different `m` could still have a different threshold.
- **ℓ:** the degree rises at ℓ = 7, whatever `m` and `N` are.
- **m·ℓ:** the degree rises once the summand unknowns `m·ℓ` pass a value
  between 18 and 21.

Four cells separate them:

| m | (n, ℓ) | N | S | m·ℓ | C(N, ≤4) | N says | ℓ says | m·ℓ says |
|--:|---|--:|--:|--:|--:|---|---|---|
| 3 | (22, 6) | 40 | +4 | 18 | 102,091 | rises | holds | holds |
| 3 | (24, 6) | 42 | +6 | 18 | 124,314 | rises | holds | holds |
| 4 | (11, 5) | 42 | −9 | 20 | 124,314 | (holds at m = 4 so far) | holds | rises if the threshold is ≤ 20 |
| 4 | (12, 5) | 44 | −8 | 20 | 149,986 | rises | holds | rises if the threshold is ≤ 20 |

- **The ℓ = 6 cells against N.** They hold ℓ fixed at a value that refuted
  before and move N to where ℓ = 7 rose.
- **The m = 4, ℓ = 5 cells.** They place `m·ℓ = 20` between the largest
  refuting value (18) and the smallest rising value (21).

## Instrument

- **Unchanged from the scale study.**
  - `examples/dreg_closure --to 4 --draw K`, one process per draw, under
    `prlimit --as=15000000000 timeout 7200`.
  - `--max-cols 155000`, raised from 130,000 only so that `(12, 5)` fits.
  - Draws come from `examples/dreg_ladder`: exact solution count, XL to
    degree 3, `--ffd-max 3`, seed `20261017`.
  - The draws are committed before any closure measurement.
- **Feasibility probe, disclosed.** Before registering, I ran
  `dreg_ladder` on m = 4 `(11, 5)`, `(12, 5)` and m = 3 `(13, 7)`,
  `(14, 7)` with a different seed (1). It ran to degree 3 only and measured
  no closure.
  - Its purpose was to check that unsatisfiable draws exist at these
    surpluses.
  - m = 3 at ℓ = 7 with small N had almost none (0 in 243 at N = 35), so
    no such cell is registered.
- **Memory.** The scale study's peak was 7.1 GB at 124k columns. `(12, 5)`
  at 150k is extrapolated to about 10–11 GB.

## Predictions

- **Q7.** Both m = 3 ℓ = 6 cells refute at ≤ 4 on every draw.
  - This rejects the N explanation and holds under the ℓ and m·ℓ
    explanations.
  - Confidence: about 70%.
- **Q8.** At m = 4 ℓ = 5, at least one of the two cells does not refute
  at 4.
  - The m·ℓ explanation predicts this with a threshold ≤ 20. The ℓ
    explanation predicts that both refute.
  - Confidence: about 50%. I hold no strong prior between the two.

## Verdict rules

Per draw and per cell, the rules are exactly those of the scale study:

- **closure ≤ 4:** the closure refutes at 3 or 4;
- **above 4:** degree 4 is built in full and saturated without refuting;
- **unmeasured:** the draw was killed. That is never evidence.
- **Cells:** **flat** when ≥ 3 of 4 draws are measured and all refute at
  ≤ 4. **rises** when any measured draw is above 4. Otherwise
  **unmeasured**.

Which explanation survives:

| outcome | reading |
|---|---|
| ℓ = 6 cells flat; m = 4 ℓ = 5 cells flat | the **ℓ explanation** survives; the m·ℓ one does not (with a threshold ≤ 20) |
| ℓ = 6 cells flat; any m = 4 ℓ = 5 cell rises | the **m·ℓ explanation** survives (threshold in 19–20 at m = 4) |
| any ℓ = 6 cell rises | the **N explanation**, or an `N`-and-`ℓ` interaction, is live; neither the ℓ nor the m·ℓ explanation survives alone |
| otherwise (cells unmeasured) | undecided, reported as such |

## Soundness control

The first two satisfiable draws the generator records in each fresh cell are
replayed with `--satisfiable --to 4`.

- Neither XL nor the closure may refute them.
- If a cell records none, the report says so.
- The m = 4 cells, with surplus −8 and −9, will record many satisfiable
  draws. The m = 3 cells, at +4 and +6, may record none.

## Scope and class

- **Class:** stage diagnostic. No `S`, no rho ratio, no speed claim.
- **Container:** 4-core, 16 GB, not isolated (`AGENTS.md` §10).
- **Against ECC2K-130:** the result selects which extrapolation the degree
  rows use. The descent route's verdict against rho is degree-independent
  (the free-oracle floor), and no outcome here reverses it.
- **Conductor:** no task was taken, because the CLI is not installed in
  this container.
