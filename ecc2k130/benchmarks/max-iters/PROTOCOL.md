# Protocol: what the step cap discards, on a scale model of the σ walk

Frozen before any run. Results go in [README.md](README.md); nothing here is
edited after the first measurement.

## Question

`aws/campaign.json` restarts a walk that has gone `maxIters = 2^30` steps
without a report and discards its trail. At `dpWeight = 32` a trail averages
`2^28.41` steps (`../dp-interval/`), so the cap is `r = 2^1.59 = 3.01` mean
trail lengths. [WALK-CONSTANT.md](../../WALK-CONSTANT.md) §6 prices that with
a model: the σ walk has no fruitless cycles, so the guard cuts only honest
trails, the longest `e^{-r} = 4.9%`, and they hold 15.6% of all trail steps.
At `2^32` (`r = 12.04`) the model gives 0.007%. Nothing has measured it.

## Hypothesis

1. **Mechanism.** The guard reports a trail iff its length `T ≤ c`, stops a
   longer one after exactly `c` steps, and changes nothing else: the next
   trail of that walk is the next seed either way, and a seed's trail is the
   same with or without a cap.
2. **Distribution.** `T` is geometric on `{0, 1, …}` with
   `θ = Σ_{even k ≤ w} C(131, k) / 2^130`, the even-weight tail of the
   subgroup's `x` coordinates. At `w = 32` this is `2^-28.409`, the fleet's
   measured interval.
3. **Consequence.** The discarded fraction of trail steps is
   `D(c) = c·q / (E[T; T ≤ c] + c·q)`, `q = (1 − θ)^(c+1)`: 15.60% at
   `2^30`, 1.47% at `2^31`, 0.0071% at `2^32` for `θ = 2^-28.41`.

## Scale model

The production walk, curve and client, with the three lengths that set the
loss shrunk together so that trails are short enough to count:

| length | production | scale model |
|---|---|---|
| DP weight → mean trail `μ = 1/θ` | 32 → `2^28.41` | 44 → `θ = 1.447744e-4`, `μ = 6907.3` |
| guard check interval `ECC_GUARD_PERIOD` | 4096 (`2^-18` of the cap) | 64 (`2^-8.3` of the smallest cap) |
| launch (`--steps`) | 1024 | 1024 |

`ECC_GUARD_PERIOD` is a compile-time macro (`include/kernel.h`); the scale
model is `make cpu` from this tree with `-DECC_GUARD_PERIOD=64` added and
nothing else changed. With caps and trail starts both multiples of 64, a cut
lands exactly at `c`, as it does in production to within `2^-18`. A walk that
reports or is cut waits for the end of its launch before it restarts; those
idle steps are counted separately and are not trail steps. Production idles
`1024 / 2^28.41 = 2^-18.4` of a trail; the scale model about 7%.

## Arms and frozen inputs

Caps are `64 · round(r · μ / 64)` at production's three ratios
`r = 2^e · 2^-28.41`:

| arm | production cap | `r` | scale cap `c` | `c / μ` |
|---|---|---:|---:|---:|
| `ref` (reference) | none | ∞ | `--max-iters 0` | ∞ |
| `c30` (baseline) | `2^30` | 3.0105 | 20800 | 3.0113 |
| `c31` | `2^31` | 6.0210 | 41600 | 6.0226 |
| `c32` (candidate) | `2^32` | 12.0420 | 83200 | 12.0452 |

Every arm runs

```
ecc2k130-cpu --curve 131 --threads 1 --steps 1024 --launches 108 \
    --dp-weight 44 --dp-cap 65536 --verify 4 --run-id R \
    --max-iters C --dp-file ARM.bin --checkpoint ARM.ck
```

with `OMP_NUM_THREADS=1`: 16,384 walks (1 thread × 32 slots × 512 lanes on an
AVX-512 host; the lane width is read back from the checkpoint), `S = 110,592`
steps per walk (`16.0 μ`, longer than the largest cap). Run ids are fixed
now: **30001** (development) and **30002, 30003** (holdouts). All twelve
runs are kept, whatever they show.

## Accounting

Unit: walk steps, one step being one `R ← R + σ^j(R)`. Every step of every
walk is assigned to exactly one of: a **reported** trail (its `iters`), a
**cut** trail (`c`), **idle** (`1024·(⌊t/1024⌋ + 1) − t` after a trail that
ended at length `t`), or the trail **in flight** at the end
(`S − startIter`, from the final checkpoint). Per walk, `ended = seed &
0xffff` from the checkpoint, `cut = ended − reported`, and

```
startIter = Σ_reported (T + idle(T)) + cut · (c + idle(c))       exactly.
```

The measured loss is `D̂ = C / (R + C)`, with `R` the steps on reported
trails and `C = cut · c`, over ended trails; the loss factor is
`(R + C) / R`. Wall time is recorded as a practicality note only.

## Checks

Any failure is a defect in the guard or the harness and stops the thread.

1. The accounting identity above holds for every walk of every run.
2. Every trail reported by a capped arm has `T ≤ c`.
3. Same run id, per seed: a seed reported by two arms has the same length
   and orbit; a seed the reference reports with `T > c` is never reported by
   the    arm with cap `c`; a seed it reports with `T ≤ c` is reported by that
   arm whenever that arm finished it (reported or cut it before the end).
4. Every record's eight witness counts sum to its length, and the four
   `--verify` rewalks per run match.
5. Every run exits 0: no report overflow, no exhausted seed counter.

## Comparison with the model

The prediction is fixed now, at `θ = 1.447744e-4` (not a fitted value): a
Monte Carlo of the exact process (16,384 walks, `S`, 1024-step launches, the
arm's cap, the idle rule, geometric trails), 400 replicas, NumPy seed
20260925, gives the mean and standard deviation of `D̂` and of the cut
fraction `q̂ = cut / ended` for each arm. The steady-state `D(c)` is reported
beside it. The pooled censored estimate
`θ̂ = reported / (Σ_reported (T + 1) + Σ_in-flight F)` over the three
reference runs is reported against `θ`.

## Success and stop

- **Success:** checks 1–5 hold on all twelve runs; for every capped arm and
  every run id, `D̂` and `q̂` lie within 3 standard deviations of the
  prediction; `c32` discards under 0.1% of trail steps.
- **Stop:** a check failure means the guard is not what the model assumes;
  the campaign cap is then not raised on the model's authority. A holdout
  outside 3 standard deviations means the geometric model is wrong for this
  walk; the production figures below are then withdrawn, not refitted.

Inadmissible: changing the weight, caps, horizon, run ids, `θ` or the
tolerance after a run; dropping or rerunning a run; counting idle steps as
trail steps or as savings.

## Production projection

A model row, marked as one: `D(c)` at `θ = 2^-28.41` for `c = 2^30, 2^31,
2^32`, with the guard overshoot (at most 3072 steps at 1024-step launches)
and launch idle (at most 1024 steps per trail) left out as below `2^-18` of a
trail. It is admissible only if the scale model succeeds. Expected work per
solve of the σ walk scales by the loss factor; `2^60.91–60.92` on completed
trails (WALK-CONSTANT.md §6) becomes `2^61.15–61.17` at `2^30`.

Expected class under AGENTS.md §3: **engineering**. The walk, its constant
and the floor do not move; a self-imposed loss is removed.

## Amendment 1, before any run

Made after checking the rule above against the rare arm and before the first
measurement. `c32` cuts a trail only if it exceeds `12.05 μ` and starts in
the first `4 μ` of the horizon: about 0.45 cuts per run, a Poisson count, for
which a 3-SD z-score flags 1.1% of honest runs. For `c32` the model
comparison is therefore the exact one-sided tail: a run fails if
`P(X ≥ cuts)` for `X ~ Poisson(λ)`, `λ` the Monte Carlo mean cut count, is
below 0.00135, the normal 3-SD tail. The `c30` and `c31` rule, the 0.1%
bound and everything else stand.

The sizing runs that chose weight 44 and 1024-step launches measured
throughput and report counts only, on run id 9 with no cap; they are not
measurements of the cap and are not used.
