# The IC/rho throughput gap: results

Written after the grid ran. [PROTOCOL.md](PROTOCOL.md) is unchanged. The readout is
[runs/grid/readout.txt](runs/grid/readout.txt), and the counts are in
[runs/grid/counts.json](runs/grid/counts.json).

## Answer

**The gap is mostly memory, and one buffer causes it.** Rho's walk (`signed_rho_fast_with`
in `koblitz_index_calculus.rs`) zero-fills a fixed cycle-detection cache on every solve:
`vec![None; RECENT_SLOTS]` with `RECENT_SLOTS = 1 << 12`, 40-byte slots, **163,840 bytes**.
That is 2,560 cache lines.

The simulation charges that function **2,619 last-level write misses** at `n19a0-000`, 72% of
rho's 3,631. It does so with 7% of rho's instructions. The size is fixed whatever the group
size. At `n ≤ 23` a walk is a few hundred steps, so this start-up write is a large share of
rho's time. It also fits the roughly 30 extra minor page faults per rho process that the
isolated re-timing measured: 163,840 bytes is 40 pages.

## Per 1,000 instructions, summed over four cases per cell

| cell | arm | mispredicts | D1 misses | LL misses | penalty cycles / instruction |
|:--|:--|--:|--:|--:|--:|
| `n13a0` | IC | 11.57 | 6.89 | 6.33 | 1.129 |
| | rho | 13.99 | 15.54 | 14.40 | 2.384 |
| `n17a1` | IC | 9.19 | 5.29 | 4.69 | 0.849 |
| | rho | 10.88 | 11.62 | 10.65 | 1.772 |
| `n19a0` | IC | 9.71 | 4.51 | 3.98 | 0.749 |
| | rho | 10.95 | 8.40 | 7.64 | 1.320 |
| `n19a1` | IC | 8.58 | 3.79 | 3.29 | 0.629 |
| | rho | 10.26 | 7.37 | 6.59 | 1.152 |
| `n23a0` | IC | 6.30 | 2.95 | 2.14 | 0.425 |
| | rho | 7.52 | 6.02 | 5.21 | 0.904 |

- **H-branch: minor.** Rho mispredicts only 1.2–2.4 more branches per 1,000 instructions,
  which is worth 0.02–0.04 cycles per instruction at 15 cycles each. Both arms' top
  mispredicting function is the same field reduction, `f2m::reduce`.
- **H-memory: dominant.** Rho has 2.0–2.5× IC's last-level misses per instruction. At 150
  cycles each, that is 0.5–1.2 extra cycles per instruction.

## Against the measured gap

The measured throughput gap `T` is the instruction ratio divided by the in-process time
ratio (IC over lean rho). `c0*` is the common base cost, in cycles per instruction, at
which the penalty model alone reproduces `T`:

| cell | `T` | `c0*`, half penalties | `c0*`, registered penalties | `c0*`, double penalties |
|:--|--:|--:|--:|--:|
| `n13a0` | 1.95 | 0.094 | 0.187 | 0.374 |
| `n17a1` | 2.07 | 0.007 | 0.015 | 0.030 |
| `n19a0` | 1.67 | 0.054 | 0.107 | 0.214 |
| `n19a1` | 1.65 | 0.089 | 0.177 | 0.354 |
| `n23a0` | 1.73 | 0.114 | 0.229 | 0.457 |

**A disclosed interpretation of the protocol.** §Reading asks whether "a common base rate"
reproduces the gap within 15%, but it did not fix that rate. I read it at the physical
bound, `c0 = 0.25` cycles per instruction, that is four instructions per cycle.

- **At `c0 = 0.25`** the model predicts `T` = 1.91, 1.84, 1.57, 1.59 and 1.71 against the
  measured 1.95, 2.07, 1.67, 1.65 and 1.73. That is within 11% at every cell, so under this
  reading **H-memory explains the gap**.
- **At a more typical `c0 = 0.5`** (two instructions per cycle), it explains 81–91% of `T`.
  That would be read as partial.
- **`n17a1` is the weakest fit.** Its `c0*` is implausibly small at any penalty scale, so
  some of its gap is not captured by the model.

## What this changes

- **The IC arm's small-cell time lead is largely rho's fixed start-up write, not IC
  efficiency.** Earlier results called IC's higher throughput unexplained; it is now
  attributed. The lead is still real for the rho as built, but it is a property of
  `RECENT_SLOTS`, not of either algorithm.
- **Sizing the cache to the walk would likely remove the lead.** For example, scale
  `RECENT_SLOTS` with `√r`, or allocate the cache lazily. This is **not measured here.**
  Changing the cache size can change when fruitless cycles are detected, and so change the
  walk. A fix needs its own regression, one that checks the answers and the walk, and an
  isolated re-timing.
- **Nothing changes for the stop decision.** Rho already wins in instructions everywhere;
  this removes the last advantage that looked like IC's own.

## Scope

- **Simulated, not measured.** Valgrind's model of this host's caches and a generic branch
  predictor. The cycle penalties are the protocol's assumed values, not measured ones.
- **One binary and five cells.**
