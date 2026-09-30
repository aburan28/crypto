# A rho with its cache sized to the walk, against the IC arm: results

This was written after the run. [PREREGISTRATION.md](PREREGISTRATION.md) is unchanged since its
registration commit `e64e2d22`. The readout is [runs/registered/readout.txt](runs/registered/readout.txt);
the regression is [runs/regression.txt](runs/regression.txt).

## Run

| | |
|:--|:--|
| binaries | round 0025 `worker` sha256 `1ca7f379…0eda`; patched worker `45f41719…efef2` ([runs/binaries.sha256](runs/binaries.sha256)) |
| wall | 2026-09-30 02:20:15Z to 02:20:53Z, pinned to CPU 3 under the benchmark lock |
| timed runs | 6,000 (5 cells × 12 cases × 5 arms × 20 rounds), 154 contended and excluded (2.6%), 0 answer mismatches |
| first attempt | refused by the preflight at 02:19:57Z and produced no run: the agent's own process used 0.32 CPU-s in 3 s against a limit of 0.30. Its output directory is kept as [runs/attempt1-refused-by-preflight](runs/attempt1-refused-by-preflight). The second attempt is the registered run. (Its `readout.txt` is only the analysis script failing on an empty run list.) |

## Registered verdict: every prediction holds

| | prediction | result |
|:--|:--|:--|
| P1 | `sized/rho`, `worker` timer, upper limit below 1 at all four small cells | **holds**: 0.755, 0.794, 0.856, 0.843, upper limits 0.777–0.879 |
| P2 | IC lead gone at ≥ 3 of 4 small cells | **holds** at `n13a0`, `n17a1`, `n19a1`; **not** at `n19a0` (see below) |
| P3 | A/A inside [0.92, 1.09] | holds: 0.997–1.029 |
| P4 | ≤ 5% contended | holds: 2.6% |
| P5 | layout control inside [0.92, 1.09] | holds: 0.963–1.008 |

## Wall-clock time, geometric mean over 12 cases, 95% bootstrap

| cell | `sized/rho` | `incumbent/rho` | **`incumbent/sized`** | A/A `sized_aa/sized` | layout `inc_new/incumbent` |
|:--|:--|:--|:--|:--|:--|
| `n13a0` | 0.875 [0.856, 0.893] | 0.862 [0.843, 0.881] | **0.986 [0.959, 1.012]** | 1.027 [1.009, 1.048] | 0.982 [0.949, 1.017] |
| `n17a1` | 0.897 [0.880, 0.912] | 0.879 [0.854, 0.903] | **0.980 [0.959, 1.000]** | 1.027 [1.004, 1.052] | 1.005 [0.979, 1.036] |
| `n19a0` | 0.913 [0.893, 0.933] | 0.895 [0.878, 0.911] | **0.980 [0.963, 0.995]** | 1.000 [0.984, 1.019] | 0.992 [0.973, 1.013] |
| `n19a1` | 0.892 [0.875, 0.913] | 0.918 [0.899, 0.941] | **1.030 [0.996, 1.060]** | 1.029 [1.003, 1.052] | 0.963 [0.949, 0.978] |
| `n23a0` | 0.944 [0.930, 0.958] | 0.985 [0.959, 1.014] | **1.043 [1.015, 1.074]** | 0.997 [0.975, 1.019] | 1.008 [0.986, 1.030] |

- **Sizing the cache recovers all of IC's small-cell time lead.**
  - The sized rho is 0.875–0.913 of the lean rho's wall time at the four small cells.
  - `incumbent/rho` reproduces the earlier isolated re-timing (0.87–0.94 there, 0.86–0.92 here).
  - Against the sized rho, IC is between 0.98 and 1.03 on wall time.
- **`n19a0` is the one small cell where IC stays ahead: 0.980 [0.963, 0.995].**
  - That is about 2%.
  - It is inside the A/A residue. `sized_aa` runs 2.7–2.9% slower than `sized` at three of the
    five cells, with intervals that exclude 1. So 2% of wall time is not distinguishable here from
    an arm-position or binary-copy effect.
  - **Prediction P2 held anyway, on the registered "at least three of four" wording.**
- **At `n23a0` the sized rho is strictly faster than IC**: `incumbent/sized` = 1.043 on wall
  [1.015, 1.074], 1.066 on `cpu`, 1.126 on the `worker` timer. Round 0025 had that cell at parity
  (`incumbent/rho` 0.985 [0.959, 1.014] here).
  - The registered "strictly faster" list covers only the four small cells, so it reads empty.
    `n23a0` had no prediction.
  - It is an ordinary 4–7% wall-time margin at a cell whose rho cache shrank from 160 KB to 80 KB.

## The mechanism is visible in page faults

Median minor faults per process:

| cell | `incumbent` | `rho` | `sized` | `rho − sized` |
|:--|--:|--:|--:|--:|
| `n13a0` | 56 | 89 | 50 | 39 |
| `n17a1` | 61 | 91 | 53 | 38 |
| `n19a0` | 64 | 94 | 59 | 35 |
| `n19a1` | 65 | 94 | 59 | 35 |
| `n23a0` | 71 | 96 | 76 | 20 |

The saving is 35–39 pages at the four small cells, where the cache falls from 40 pages to 1–5.
At `n23a0` it is 20 pages, half of the old 40, which is what an 80 KB cache predicts. This is
the fixed start-up write the [throughput diagnostic](../throughput_gap_20260929/RESULTS.md) found
by simulation. It is now measured.

## Regression: R1–R3 hold

1,536 rho runs over round 0025's development, confirmation, replay, aa, smoke and selection
directories, against the recorded lean-rho answers.

- **R1.** All 1,536 verified, and all 1,536 recovered the recorded logarithm.
- **R2.** 1,509 runs had the identical walk. The 27 that differ are 15 runs at `n = 23` and 12 at
  `n = 31`. **All 60 timed cases had the identical walk** (`iterations` equal to the lean rho's,
  and `sized_aa` equal to `sized`), so the timing compares equal walks.
- **R3.** Where the walk differs, total iterations rise by 0.12% at both degrees, and no restart
  count changed.

## Scope, and what this does not say

- **A constant, in a measured cell range.** Five cells, `n` 13–23, one host class, one round's
  minimal-startup worker (AGENTS.md §10 conditions). Nothing here concerns how either cost grows
  with `n`.
- **Instruction counts are unchanged in ranking.** Rho was already lower in instructions at every
  cell ([RESULTS round 0025](../RESULTS.md#round-0025-a-lean-rho-and-the-pairtriple-switch)). This
  round removes the one place IC looked better, its time. It does not create an IC result, and
  the exponent question is untouched.
- **The patch changes the walk at some cells** (`n = 23`, `n = 31`), by 0.12% in total iterations.
  It is a research candidate for the tournament's rho, not a change to the library's rho.
- **`n19a0` stays 2% in IC's favour** on wall time, at a noise level the A/A cannot resolve.
- **Stock builds.** A `std` rho binary has its own fixed start-up cost, so this says nothing about
  the ratio for one.

## What this changes

- **The IC arm's small-cell time lead was rho's fixed start-up write.** It is gone once that
  write is sized to the walk, and the sized rho is ahead at `n23a0`. Classified under AGENTS.md §3
  as **engineering**, and it is a rho improvement, not an IC one.
- **The stop decision's loose end is closed.** Its statement that the sized rho "has not been
  built or measured" is superseded by this record, and its reopening item for a stronger rho at
  small `n` is met for the time lead.
- **For any reopened tournament round,** the reference rho at `n ≤ 23` is the sized-cache rho, not
  the lean rho. The IC arm must be timed against it on the isolated runner.
