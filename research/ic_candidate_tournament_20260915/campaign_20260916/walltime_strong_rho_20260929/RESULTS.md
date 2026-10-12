# Small-cell re-timing against main's strongest rho: results

Written after the run. [PREREGISTRATION.md](PREREGISTRATION.md) is unchanged since commit
`50158845`. The raw records are [runs/registered/runs.jsonl](runs/registered/runs.jsonl),
and the readout is [readout.txt](runs/registered/readout.txt).

## Answer

**The small-cell lead survives the strongest rho on main.**

- **Against the strongest rho.** In-process, the IC arm is 0.74–0.84 of whichever rho arm is
  fastest at each of the four small cells. Every upper 95% limit is below 1, so **P2 holds**.
- **At the control cell.** At `n23a0` the two are at parity: 1.03 [0.96, 1.11].
- **The batched rho does not help here.** At these sizes it is no faster than the lean rho,
  and 32 lanes are no faster than 1.

This contradicts what I said I expected before the run, namely that the lead would vanish
against a batched rho. The pre-registration itself made no prediction on that point.

## Conditions

- **Runs.** 3,000 timed runs, with 90 (3.0%) contended and excluded, so **P4 holds**.
- **Answers.** No answer failures.
- **A/A.** `strong32` against itself reads 0.979–0.999 in-process, with every interval
  inside [0.955, 1.018], so **P3 holds**.

## In-process time

Geometric mean over 12 cases, in ms:

| cell | IC | lean rho | strong, 32 lanes | strong, 1 lane | strongest | IC / strongest |
|:--|--:|--:|--:|--:|:--|--:|
| `n13a0` | 0.281 | 0.380 | **0.376** | 0.378 | strong32 | **0.748** [0.731, 0.764] |
| `n17a1` | 0.309 | **0.415** | 0.442 | 0.443 | lean | **0.744** [0.705, 0.785] |
| `n19a0` | 0.343 | **0.432** | 0.631 | 0.625 | lean | **0.793** [0.766, 0.824] |
| `n19a1` | 0.380 | **0.451** | 0.616 | 0.630 | lean | **0.842** [0.807, 0.883] |
| `n23a0` | 0.507 | **0.492** | 0.759 | 0.809 | lean | 1.031 [0.964, 1.107] |

**Whole-process wall time, reported, not primary:**

- **Between the worker's modes.** IC/lean rho is 0.875, 0.883, 0.900 and 0.911 at the four
  small cells, and 1.004 at `n23a0`, in line with the first isolated re-timing.
- **The strong rho.** It takes 3.5–4.0 ms per process against about 1 ms for the worker,
  because it is a standard `std` binary. That is process start-up, not the walk.

## What it means

- **The batching cure does not apply at these sizes.**
  - The strong rho's advantages come from 32 lanes sharing one inversion and from the
    normal-basis representative. Both are measured to pay at `n = 37–53` and `L = 1,024`
    targets (#943, #966).
  - At `n ≤ 23` with one target, its setup costs more than they save, and lockstep lanes buy
    nothing: `strong1` ≈ `strong32` everywhere.
  - So this run does not support the latency-bound-walk hypothesis in the form a batched
    rho would fix. The IC arm's higher instructions per second stays unexplained.
- **What IC has is a constant-factor lead in time**, 16–26% in-process, at the four smallest
  cells.
  - It shrinks with `n` and is gone by `n = 23`.
  - Rho wins in instructions at every cell, and in time at every larger cell of round 0025.
  - This is **engineering, not an advance** (AGENTS.md §3). It does not move any exponent.

## Scope

- **Where.** One host class: a 4-vCPU Intel Xeon @ 2.80 GHz cloud container. Five cells at
  `n ≤ 23`.
- **Binaries.** One IC binary, and main's strong rho at `KIC_RHO_DP_BITS=2`, because its
  default setting aborts at `n = 13, 17`.
- **Targets.** The strong-rho arms draw their own targets, so they are paired with the
  worker's cases by index, not by target.
