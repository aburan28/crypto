# Results: solving degree at fixed surplus (`m = 3`), scored as far as it ran

**Verdict, by the registered rule: inconclusive.** The primary pair's large
cell `(13, 5)` has no finished draw. The two secondary pairs that can be
scored both read **grows**. They share a confound the design did not
anticipate, stated below; it cuts against reading them as field-size growth.

Scored on 2026-09-24, at the user's direction, from the committed runs only
(`score.py`, `score-output.txt`). The pre-registration and its four addenda
are in `PREREGISTRATION.md`.

## What ran

Every measured draw was verified unsatisfiable by exact count before it was
measured. **Every measured draw resolved by refutation.** None was pinned,
none hit the caps, and none needs a lower bound.

| cell `(n, ℓ)` | `S` | unknowns, equations | `d_max` | resolving degree (4 unsat draws) | FFD | control | time a draw |
|---|--:|---|--:|---|---|---|--:|
| `(5, 2)` | `−1` | 11, 10 | 7 | 5 5 5 5 | 3 2 3 3 | not resolved by 7 | < 0.1 s |
| `(7, 2)` | `+1` | 13, 14 | 7 | 5 5 5 5 | 3 3 3 3 | 6, pinned | 0.1 s |
| `(7, 3)` | `−2` | 16, 14 | 7 | 6 6 6 6 | 2 3 3 3 | 7, pinned | 3–8 s |
| `(9, 3)` | `0` | 18, 18 | 7 | 6 6 6 6 | 3 3 3 3 | stopped unfinished (user's direction) | 67–80 s |
| `(11, 4)` | `−1` | 23, 22 | 6 | 6 6 6 6 | 3 3 3 3 | not resolved by 6 | 21–27 min |
| `(13, 4)` | `+1` | 25, 26 | 6 | 6 6 6 6 | 3 3 3 3 | lost to a container reboot | 1.7–2.0 h |
| `(13, 5)` | `−2` | 28, 26 | 6 | **none finished** | — | — | > 4.5 h |
| `(15, 5)` | `0` | 30, 30 | 6 | **never started** | — | — | — |

Times are wall-clock on the four-core container: a practicality note, never
the metric.

**What stopped `(13, 5)` and `(15, 5)`.** The container restarted four times
between about 05:17 and 14:14 UTC on 2026-09-24. That includes an outage
from about 05:38 to 13:42, and from 13:42 it restarted on every wake from
idle.

- `(13, 5)`'s first draw was killed about 4.5 hours in. Three later
  relaunches died within minutes to hours, and the fourth was stopped when
  scoring began.
- A draw that needs 4.5 uninterrupted hours cannot finish in an environment
  that restarts on every idle wake.
- These are **resource limits**. They are recorded as such, and they are not
  evidence about any degree.

## Scored by the registered rule

| pair | `S` | small | median | large | median | verdict |
|---|--:|---|--:|---|--:|---|
| **primary** | `−2` | `(7, 3)`: 6 6 6 6 | 6 | `(13, 5)`: — | — | **not testable**: no committed draws |
| | `−1` | `(5, 2)`: 5 5 5 5 | 5 | `(11, 4)`: 6 6 6 6 | 6 | **grows** |
| | `0` | `(9, 3)`: 6 6 6 6 | 6 | `(15, 5)`: — | — | **not testable**: no committed draws |
| | `+1` | `(7, 2)`: 5 5 5 5 | 5 | `(13, 4)`: 6 6 6 6 | 6 | **grows** |

**Overall: inconclusive.** The rule makes that the verdict whenever the
primary pair cannot be tested. Every testable pair points the predicted way,
and none reads flat or falls.

## The confound, found after the run

This was not registered, it does not change the verdict, and it is the most
important thing in this note for whoever runs `(13, 5)` next.

- **The two pairs that grow are exactly the two whose small cell has `ℓ = 2`.**
  An `ℓ = 2` subspace has four elements, so each summand carries two unknowns.
  Both `ℓ = 2` cells resolve at **5**.
- **Every cell with `ℓ ≥ 3` resolves at 6**, from 16 to 25 unknowns and
  across surpluses `−2, 0, −1, +1`: `(7, 3)`, `(9, 3)`, `(11, 4)` and
  `(13, 4)`, all 6 6 6 6.

So the same data fit two readings equally well:

1. **At fixed surplus the solving degree grows with the field.** This is the
   registered reading of the two scored pairs.
2. **The degree is 5 at `ℓ = 2` and 6 at every `ℓ ≥ 3` measured, and does
   not move between 16 and 25 unknowns.** On this reading, "grows" is an
   `ℓ = 2` floor effect, and the solving degree was flat across the range.

The design cannot separate these, because holding `S = n − 3ℓ` fixed while
`n` grows forces `ℓ` to grow too. The primary pair is the one that can:
`(7, 3) → (13, 5)` has `ℓ ≥ 3` at both ends. **It is the decisive
follow-up, and it is the pair this environment could not finish.**

The `ℓ ≥ 3` cells alone are not a registered comparison, since no two share a
surplus. Read descriptively, they are the first observation in this
repository of the solving degree **not** moving over a range of unknowns:
6 at 16, 18, 23 and 25. That is consistent with Petit–Quisquater's
constant-gap assumption over this range, and it establishes nothing about
larger sizes.

## Secondary, not decided on

- **FFD is 3 on 22 of 24 draws, and 2 on the other two.** The gap to the
  solving degree is 2 at the `ℓ = 2` cells and 3 at every `ℓ ≥ 3` cell. This
  matches the DREG note's gap of 3 at `n = 5` and `n = 7`.
- **Wherever a control finished, the Semaev structure resolves at least one
  degree lower:**

  | cell | Semaev draws | random control |
  |---|--:|---|
  | `(5, 2)` | 5 | not by 7 |
  | `(7, 2)` | 5 | 6, pinned |
  | `(7, 3)` | 6 | 7, pinned |
  | `(11, 4)` | 6 | not by 6 |

- **Random subspaces reproduce the invariant-subspace result.** The
  random-subspace `(7, 3)` cell resolves at 6 on all four draws, as the
  invariant-subspace Result 3 did.
- **Satisfiable fractions** were 8/12 at `S = −2`, 3/7 and 2/6 at `S = −1`,
  0/4 at `S = 0`, and 1/5 and 1/5 at `S = +1`. The yield falls with surplus,
  as `λ = 2^{−S}/3!` says it should. The samples are too small for more.

## Scope

- **`m = 3`, the chained `S₃` system, `b = 1`, random subspaces, `n ≤ 13`
  measured, four draws a cell.**
- **The solving degree is that of Macaulay-matrix linear algebra** (sparse
  elimination).
- **Class: stage diagnostic.** It computes no `S` and no rho ratio, so it
  owes the scoreboard no row, as registered.
- It says nothing about `n = 131` beyond what is stated above, at this scale.

## To finish it

Build commit `968b6cf`, or any later commit on this branch. Then run:

```sh
python3 research/dreg_fixed_surplus_20260923/run_queue.py <dreg_ladder binary> <workdir>
```

- It runs the primary cell's four draws first, then `(15, 5)`'s, then the
  controls, one process per draw.
- It resumes after an interruption, skipping every finished job.
- It needs about 10 GB and two cores. Each `(13, 5)` draw needs **more than
  4.5 uninterrupted hours** on a four-core container.
- The binary is the one checked against the frozen `f03dc02` rows (0
  mismatches in 18; `runs/identity-check-*.jsonl`).

When `(13, 5)` has three or more finished draws, `score.py` scores the
primary pair by the registered rule with no other change.

## Addendum, 2026-09-25: the cost of finishing, re-estimated

Additive; nothing above is changed. "More than 4.5 uninterrupted hours" a
`(13, 5)` draw is true as a lower bound and badly understates the cost.

- **The model.** Time is fitted as `(rows × cols)^α` against this study's
  own measured draw times: `(11, 4)`, `(13, 4)` and the degree-5 probe.
- **The estimate.** One `(13, 5)` draw takes about **11–24 h** on the
  four-core container, and one `(15, 5)` draw about **1.7–5.6 days**.
- **It is an extrapolation.** It is used for scheduling only. The script
  and its output are `research/dreg_ell_grid_20260925/cost_model.py` and
  `cost-model-output.txt`.

The primary pair still needs a large machine that stays up for days.
Meanwhile, `research/dreg_ell_grid_20260925/` pre-registers a cheaper design
that separates the two readings in "The confound" above. It compares
`ℓ = 2` with `ℓ = 3` at matched unknown counts, and adds an `ℓ = 5` cell,
`(10, 5)`, at `(13, 4)`'s 25 unknowns.

## Addendum, 2026-09-26: the primary pair, scored — grows at fixed surplus

Additive; nothing above is changed. `(13, 5)` now has four finished draws,
so `score.py` scores the primary pair by the registered rule, with no other
change. The output is in `score-output-20260926.txt`. The verdict above,
"inconclusive", stands for the state it described. **The registered verdict
with the primary pair in is "grows at fixed surplus".**

| pair | `S` | small | large | verdict |
|---|--:|---|---|---|
| **primary** | −2 | `(7, 3)`: 6 6 6 6 | **`(13, 5)`: ≥7 ≥7 ≥7 ≥7** | **grows** |
| | −1 | `(5, 2)`: 5 5 5 5 | `(11, 4)`: 6 6 6 6 | grows |
| | 0 | `(9, 3)`: 6 6 6 6 | `(15, 5)`: running | not testable |
| | +1 | `(7, 2)`: 5 5 5 5 | `(13, 4)`: 6 6 6 6 | grows |

- **How it was measured.** The four `(13, 5)` draws were measured on the
  dense-finish path under the pre-registration's last addendum. That path
  was identity-checked on 376 committed rows with 0 mismatches.
  - The draws were draw indices 0, 2, 4 and 7, so 4 of 8 were
    unsatisfiable.
  - Each took 31–34 min on the four-core container.
  - FFD is 3 on all four draws.
  - Every draw is a mathematical lower bound: the degree-6 Macaulay matrix
    does not contain `1`. So "grows" here means at least one degree.
- **The confound "The confound" above named is answered.** The primary pair
  has `ℓ ≥ 3` at both ends and is held at one surplus. Going from 16 to 28
  unknowns at `S = −2`, the degree goes from 6 on every draw to above 6 on
  every draw.
  - So the growth is not an `ℓ = 2` floor effect.
  - The surplus control (`research/dreg_surplus_control_20260925/`) showed
    that the surplus can move the degree at `ℓ = 4`. It is fixed within
    this pair, so that confound does not apply here.
- **The prediction was "grows", and it held for the primary pair.**
- **What this does not say.**
  - How fast the degree grows: the `(13, 5)` values are lower bounds.
  - Anything at `S = 0`, which is `(15, 5)`, running now.
  - Anything about `n = 131`.
- **Class:** stage diagnostic, as registered, so no scoreboard row.

## Addendum, 2026-09-26 (later): the `S = 0` pair — every pair grows

Additive; nothing above is changed. `(15, 5)` now has four finished draws,
measured on the same dense-finish path. `score.py`, unchanged, scores every
pair; the output is in `score-output-20260926b.txt`.

| pair | `S` | small | large | verdict |
|---|--:|---|---|---|
| **primary** | −2 | `(7, 3)`: 6 6 6 6 | `(13, 5)`: ≥7 ≥7 ≥7 ≥7 | **grows** |
| | −1 | `(5, 2)`: 5 5 5 5 | `(11, 4)`: 6 6 6 6 | grows |
| | **0** | `(9, 3)`: 6 6 6 6 | **`(15, 5)`: ≥7 ≥7 ≥7 ≥7** | **grows** |
| | +1 | `(7, 2)`: 5 5 5 5 | `(13, 4)`: 6 6 6 6 | grows |

**Overall, by the registered rule: grows at fixed surplus, with every pair
testable and every pair growing.**

- **The `(15, 5)` draws.** They are draw indices 0–3, so every draw was
  unsatisfiable. That fits `S = 0`, where the expected yield is `1/3!`.
  - Each draw took 96–104 min on the four-core container, and FFD is 3 on
    all four.
  - Each value is a mathematical lower bound: the degree-6 Macaulay matrix
    over 30 unknowns and 30 equations contains no `1`.
- **Why `S = 0` matters.** It is the regime index calculus works in, with
  about one decomposition in six points. From 18 unknowns to 30, the
  refutation degree goes from 6 on every draw to above 6 on every draw.
- **Execution, disclosed.**
  - The first attempt at draw 0 was lost to a container restart while
    idle. That is a resource limit, not evidence, and the draw was rerun
    from the same seed.
  - Draws 2 and 3 were launched by hand with `run_queue.py`'s exact command,
    staggered so their dense phases, of about 6.7 GB each, did not overlap.
  - `runs/queue-dense.log` records every launch. Outcomes depend only on
    the seed.
- **The three controls** (`(13, 5)`, `(15, 5)` and `(13, 4)`) are
  secondary. They are running now and will be recorded if they finish.
- **Scope.**
  - The large-cell values are lower bounds, so the size of the growth is
    unmeasured.
  - It covers `m = 3` and `n ≤ 15`, and says nothing about `n = 131`.
  - Class: stage diagnostic, so no scoreboard row.
