# Round 1 results: does the semi-regular law predict refutation degrees it has not seen?

Written after the run. [PREREGISTRATION.md](PREREGISTRATION.md) is unchanged since its
registration commit `1c1092005`. The readout below is printed by [analyze.py](analyze.py),
which applies §5 mechanically; it is saved as
[runs/registered/readout.md](runs/registered/readout.md).

## Requested and delivered

| requested | status | evidence |
|:--|:--|:--|
| A goal for a significant step, set before the round | done | [GOAL.md](GOAL.md): replace the fitted `ℓ + 5` ladder with a law that predicts unseen cells |
| Run the round | done: 68 of 68 registered units resolved | [runs/registered](runs/registered), remote lanes `c`, `d`, `e` |
| Predictions fixed before measurement | done | `predictions.jsonl` and `dump.sha256` in commit `1c1092005`, before any registered draw ran |

## Registered verdict: neither

By §5 the verdict is **neither**. The fixed-`ℓ` fit fails on every discriminating cell:
none of the seven reads the fit's value. The generic law holds on six of the seven, and
fails on one, in the unexpected direction: the symmetric form `rr` on `K₁/2⁶¹ ℓ = 5`
refutes at 5, below its window `[6, 7]`. The rule requires every discriminating cell, so
the law is not confirmed as registered.

| cell | arm | readings | window | fit | in window |
|:--|:--|:--|:--|:--|:--|
| K0n61l4 | rr | 5 5 5 5 | (5, 6) | 6 | yes |
| K0n61l4 | x4 | 8 8 8 8 | (8, 9) | 9 | yes |
| K0n61l6 | rr | 6 6 6 6 | (6, 7) | 8 | yes |
| K0n61l6 | x4 | 10 10 10 10 | (9, 10) | 11 | yes |
| K1n23l5 | rr | 7 7 triv 7 | (7, 8) | 8 | yes |
| K1n23l5 | x4 | 10 10 10 10 | (9, 10) | 10 | yes |
| K1n23l6 | rr | 8 8 8 8 | (7, 8) | 8 | yes |
| K1n29l5 | rr | 6 6 6 7 | (6, 7) | 8 | yes |
| K1n29l5 | x4 | 10 10 10 10 | (9, 10) | 10 | yes |
| K1n45l5 | rr | 6 6 6 6 | (6, 7) | 8 | yes |
| K1n45l5 | x4 | 9 9 9 9 | (9, 10) | 10 | yes |
| K1n59l6 | rr | 6 6 6 6 | (7, 8) | 8 | below |
| K1n59l6 | x4 | 10 10 10 10 | (10, 11) | 11 | yes |
| K1n61l5 | rr | 5 5 5 5 | (6, 7) | 8 | below |
| K1n61l5 | x4 | 9 9 9 9 | (9, 10) | 10 | yes |
| K1n61l6 | rr | 6 6 6 6 | (6, 7) | 8 | yes |
| K1n61l6 | x4 | 10 10 10 10 | (9, 10) | 11 | yes |

Resolved draws: 68 of 68 (none censored). Inside the window: 60. Below: 8. Above: 0. P1 fails (88% inside). P2 fails on one of seven discriminating cells (`K₁/2⁶¹ ℓ = 5` `rr`). Discriminating cells whose median equals the fit: 0 of 7.

## The two arms behave differently

**The direct `S₄` descent (`x4`) follows the generic law on every draw.** All 32
resolved `x4` draws are inside `[d_reg, d_reg + 1]`. On six cells the reading is below
the fit, and in each case it is where the series puts it:

| cell | `x4` reading | `d_reg` | fit (`ℓ + 5`) |
|:--|:--|:--|:--|
| `K₀/2⁶¹ ℓ = 4` | 8 | 8 | 9 |
| `K₁/2⁴⁵ ℓ = 5`, `K₁/2⁶¹ ℓ = 5` | 9 | 9 | 10 |
| `K₁/2⁵⁹ ℓ = 6` | 10 | 10 | 11 |
| `K₀/2⁶¹ ℓ = 6`, `K₁/2⁶¹ ℓ = 6` | 10 | 9 | 11 |

So the `x4` degree is not a function of `ℓ` alone: at fixed `ℓ` it falls as `n` grows,
by the amount and at the `n` the series predicts. It reads `d_reg` on some cells and
`d_reg + 1` on others; the window does not say which.

**The symmetric norm form (`rr`) does not.** Eight `rr` draws on two cells read below the
window, and none above:

| cell | `rr` readings | window |
|:--|:--|:--|
| `K₁/2⁵⁹ ℓ = 6` | 6 6 6 6 | `[7, 8]` |
| `K₁/2⁶¹ ℓ = 5` | 5 5 5 5 | `[6, 7]` |

## Why `rr` departs (post-hoc diagnostic; not evidence for the verdict)

[topforms.py](topforms.py), run after the round on every held-out draw
([runs/registered/topforms.jsonl](runs/registered/topforms.jsonl)), shows two things.

1. **Half of `rr`'s equations are redundant.** Of its `2n` quartic equations, exactly `n`
   are linearly independent on every draw of every cell; the other `n` are sums of them.
   The registered prediction counted all `2n`. A system with `n` quartics and one linear
   equation has a higher generic degree.
2. **`rr` is easier than a generic system of its true shape.** Against the generic degree
   of the reduced system (`n` quartics, one linear equation), `rr`'s median reads equal on
   the two `n = 23` cells and one below on all seven cells from `n = 29` up; one draw
   there (`K₁/2²⁹ ℓ = 5`, draw 3) reads equal. A reading below the generic degree means the top-degree
   forms have relations a generic system does not have. That is structure, and it helps
   the solver.

`x4` has no such redundancy: its top forms have full rank on every cell except
`K₀/2⁶¹ ℓ = 4`, where 49 of 61 are independent; its readings sit at `d_reg` or
`d_reg + 1` of the face-value series throughout.

The registered law happened to fit `rr` better than the corrected one would have, because
counting each equation twice roughly offset `rr`'s structure. That agreement was an
accident, and the diagnostic is post-hoc; testing a corrected `rr` predictor needs a new
registration and fresh draws.

## What this does to the goal

- **The fitted ladder is retired.** "`x4` refutes at `ℓ + 5`" fails out of sample on six
  cells; the degree depends on `n`.
- **For `x4`, the generic law is the best predictor on record.** It held on every one of
  32 held-out draws, across `n` 23–61 and `ℓ` 4–6, including all six cells where it
  disagreed with the fit. The round does not confirm it as registered, because the
  registration bundled both arms.
- **For `rr`, the generic law is not the right model.** The departure is structural and in
  the attacker's favour, and its size (one degree, from `n = 29` up) is measured only at
  `m = 3`.
- **Cost at `n = 131`** ([cost_n131.py](cost_n131.py), conditional on the window): index
  calculus at `m = 3` stays at least `2^74` above rho for every `ℓ`, because collecting
  relations needs about `2^{131 − 2ℓ}` attempts. Neither arm's degree changes that, and
  this round claims nothing about `m ≥ 4`.

## Run

| | |
|:--|:--|
| engine | `examples/macaulay_dense.rs`, unchanged |
| systems | exported at seed 20261005 and checked against `dump.sha256`; regenerated on every remote host |
| local lanes | [lanes.sh](lanes.sh), two lanes × two threads, 2026-10-05 16:58Z to 2026-10-11 23:27Z |
| remote lanes | `c` and `d` (`K₁/2⁶¹ ℓ = 6` `x4`, all four draws) and `e` (`K₁/2⁵⁹ ℓ = 6` `x4`, draws 0–1), each on a separate 4-core host |
| heavy solves | `x4` at degree 10 on `ℓ = 6`: 239k–247k rows × 199k columns; 42–51 min on four threads, about 90 min on two |

### Departures and incidents, disclosed

- **Smoke test on a held-out cell.** Disclosed in the registration; `K₁/2⁶¹ ℓ = 4` was
  dropped and its two readings are kept apart in
  [runs/pre-registration-smoke](runs/pre-registration-smoke).
- **Three remote lanes stopped.** Lanes `a`, `b` and `f` stopped about 35 minutes in, when
  the account reached its weekly usage limit (2026-10-05). They had not measured past
  degree 8. Their six units (`K₀/2⁶¹ ℓ = 6` `x4` draws 0–3, `K₁/2⁵⁹ ℓ = 6` `x4` draws 2–3)
  were re-run from degree 1 by the local lanes. Lane `b`'s partial log is kept in
  [runs/registered-remote-b](runs/registered-remote-b) and is not read.
- **Two remote lanes merged their own pull requests.** Lanes `c` and `d` were told not
  to open pull requests. Each opened one and merged it (#1666, #1665), adding only its own
  run directory and copies of the registration files identical to `1c1092005`.
- **One host restart.** The local container restarted between 2026-10-05 17:05Z and
  2026-10-11; the lanes resumed from their logs (recorded in `progress.txt`). The
  interrupted processes left no line and are not counted.
- **Parallel lanes.** Units ran concurrently on several hosts, as registered. Degrees
  are deterministic; wall times are advisory.

## Scope

Two Koblitz curves, `n` 23–61, `ℓ` 4–6, `m = 3`, four rootless draws per cell, one engine.
A stage diagnostic about the degree of one decomposition step. No end-to-end cost is
measured, nothing is claimed for `m ≥ 4`, and the `n = 131` figure is a model output
conditional on the window.

## Next

Round 2 can test the two findings separately, each on fresh draws:

1. **`x4` alone.** Register the window for `x4` only, at cells where it moves (larger `n`
   at `ℓ = 5, 6`, and `ℓ = 4` at small `n`), and add a rule for when the reading is
   `d_reg` versus `d_reg + 1`.
2. **`rr`'s structure.** Register the reduced-shape prediction minus one for `n ≥ 29`,
   and look for the extra top-form relations directly. They are the part of this round
   that points in the attacker's direction.
