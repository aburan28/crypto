# Results: degree 7 at `ℓ = 5` (`m = 3`)

Scored by `score.py` against `PREREGISTRATION.md` and its amendments 1 and 2.
Every value comes from a committed run record in `runs/`.

## The measurement

| cell | `N` | `S` | `ℓ = 3` | `ℓ = 4` | **`ℓ = 5` at degree 7** | reading | semi-regular `D_reg` | verdict |
|---|--:|--:|---|---|---|---|--:|---|
| `(10, 5)` | 25 | −5 | 5 5 5 5 | 7 7 6 7 | **U3_CELL** | **7** | 7 | **Q6: one degree above 6** |
| `(13, 5)` | 28 | −2 | 6 6 6 6 | — | not run (about 47–50 GB) | not measured | 7 | Q5: not testable |
| `(11, 5)` | 26 | −4 | 6 6 6 6 | 6 7 6 7 | not run (about 18–19 GB) | not measured | 7 | secondary, not decided on |

The ℓ = 3 and ℓ = 4 columns are the committed values at the same surplus,
from Results 5 and 6 of `RESEARCH_DREG_MEASUREMENT.md`.

- **ℓ = 5 refutes at exactly 7: one degree above 6.**
  - This is the first exact ℓ = 5 value. Every earlier ℓ = 5 draw was a ≥7
    bound.
  - The degree-7 row space contains `1` on U3_COUNT four draws. The FFD is 3 on
    every draw.
  - Three draws already fix the registered median at 7, whatever the last
    one reads.
- **The predictions held.**
  - Q6 was predicted at 7, with low confidence: ≥8 was named as the likelier
    miss, and it did not happen.
  - Q5, `(13, 5)`, is not testable. It needs about 47–50 GB, and this
    machine has 16 GB (amendments 1 and 2).
- **The semi-regular reference held.** It predicted 7 for `(10, 5)`, and
  that is the measured value. Across the 14 cells with exact values it now
  equals the measured value in 10, and is within one in all 14.
- **The size of the growth, read down the `S = −5` column:**
  - ℓ = 3 → 5, ℓ = 4 → 7 7 6 7, ℓ = 5 → 7.
  - So there are two degrees between ℓ = 3 and ℓ = 5 at this surplus.
  - From ℓ = 4 to ℓ = 5 the degree does not rise, except on the one ℓ = 4
    draw that read 6.
- **Against the ℓ = 3–4 plateau of 6 at `S ≥ −2`:** exactly one degree more.
  - Whether `S ≥ −2` also reads 7 at ℓ = 5 is the unrun `(13, 5)`.
  - The reference predicts 7 there, and 8 at `(15, 5)`.
- **What it does not separate.** 7 at ℓ = 5 fits both a degree that has
  levelled off at 7 and one that follows the semi-regular reference, which
  is 7 at 25 unknowns and grows about linearly in them. The cell that tells
  them apart is `(15, 5)`, where the reference is 8. It needs a larger
  machine.

## The runs

| draw | outcome | FFD | wall time | peak memory | binary, budget |
|---|---|--:|--:|--:|---|
| `(10, 5)` u0, draw 22 | **7**, refuted | 3 | 9,376 s | 10.99 GB | `7eb043b0`, 11,000 MB, F5 rows |
| `(10, 5)` u1, draw 34 | **7**, refuted | 3 | 9,055 s | 10.96 GB | `7eb043b0`, 11,000 MB, F5 rows |
| `(10, 5)` u2, draw 35 | **7**, refuted | 3 | 9,085 s | 11.00 GB | `7eb043b0`, 11,000 MB, F5 rows |
| `(10, 5)` u3, draw 57 | U3_ROW |

- **Draws.** Each row replays its committed degree-6 draw from `--d-min 7`.
  `score.py` checks the draw index, subspace, target and the `≥7`
  premise of every row, and voids none.
- **Engine.** Sparse elimination of band 7, then a dense M4RI finish over
  the remaining 245,506 columns. It uses the F5-criterion rows of
  amendment 1: the same row space, so the same outcome.
- **Wall times are practicality notes** on a four-core, 16 GB container,
  single-threaded in the sparse phase and multi-threaded in the dense one.
  No speed claim is made.

## What did not go to plan, in order

1. **Band 7 is rank-deficient.**
   - The pre-registration sized the block after band 7 assuming full rank:
     356k rows, 10.2 GiB. More rows survived.
   - The first `(10, 5)` u0 run therefore did not switch at 11,000 MB.
     It went on through band 6 sparsely, at under one column a second.
   - It was stopped by hand after 4 h 26 min with no outcome
     (`runs/stopped/NOTE.md`). That is a resource limit, not evidence.
2. **F5-criterion rows (amendment 1).**
   - They drop 171,003 of the 836,820 degree-7 rows on that draw, and as
     many band-7 survivors.
   - Identity check: 332 committed rows, 0 mismatches
     (`runs/identity-check-f5/`).
   - The rerun switched to dense right after band 7:
     - band 7 took about 1 h 40 min;
     - the dense block used about 10.9 GB of the 11,000 MB budget.
3. **`(13, 5)` and `(11, 5)` are out of reach here.**
   - `(10, 5)`'s switch within budget bounds its band-7 rank at 60–65% of
     the band.
   - At that ratio the block after band 7 is about 47–50 GB at `(13, 5)`
     and 18–19 GB at `(11, 5)`.
   - Both are therefore not measured (memory), under amendments 1 and 2.
     Amendment 1's "≥20.7 GB" for `(13, 5)` assumed full rank.
   - They need hosts with about 64 GB and 24 GB. The runner was restarted
     without them. u1 had run a minute, and restarted from the beginning.

## Engineering record

- **Commits.**
  - `d6c1ee51`: the budgeted multi-band switch and `--d-min`.
  - `0227cd61`: frees the pivot rows.
  - `7eb043b0`: F5 rows.
- **Identity checks.** Each passed on 332 committed rows (41 measured draws)
  with 0 mismatches, before the binary it checked ran any `ℓ = 5` draw:
  `runs/identity-check/` (`0227cd61`) and `runs/identity-check-f5/`
  (`7eb043b0`).
- **Class.** Engineering: no outcome moved. The `(8, 4)` degree-7 draw took
  284 s at 1.31 GB without the F5 rows and 231 s at 0.85 GB with them. That
  is one draw on a shared container, a practicality note and not a speed
  claim.

## Scope

- **System.** `m = 3`, the chained `S₃` system, `b = 1`, random subspaces.
  - `(10, 5)` lives over GF(2^10), which has the proper intermediate
    subfields GF(2²) and GF(2⁵). `(11, 5)` lives over GF(2^11), which has
    none. No subfield structure is chosen.
  - Four unsatisfiable draws a cell.
- **Degree.** The refutation degree of bounded Macaulay linear algebra: the
  constant `1` in the row space, or every variable pinned. It is not an
  F4/F5 solving degree, and not a certified first fall degree.
- **Class.** Stage diagnostic: no `S`, no rho ratio. What this means for
  ECC2K-130 is in
  `research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md`, where
  every `n = 131` figure is marked as extrapolation.
