# What the step cap discards: the σ walk on a scale model

**Result.** The client's guard cuts exactly the trails longer than the cap
and changes nothing else, and trail lengths follow the geometric law the
fleet's report interval implies. At the campaign's old cap-to-trail ratio
(`maxIters = 2^30` at `dpWeight = 32`, 3.01 mean trail lengths) the scale
model discards the share of trail steps the model predicts, on the
development run id and on both holdouts. So the model's production figure
stands: **2^30 discarded 15.6% of the σ walk's steps, ×1.185 expected work;
2^32 discards 0.007%**. `aws/campaign.json` now sets `2^32`.

[PROTOCOL.md](PROTOCOL.md) was committed before any run, and its one
amendment, for the rare arm, before the first measurement. Nothing below was
tuned after a run.

## Boundary, unit and reference

- **Unit:** walk steps (one step is one `R ← R + σ^j(R)`); the table's
  column is the share of trail steps a guard discards,
  `D = cut steps / (reported + cut steps)`.
- **Floor:** `D = 0`, loss ×1.000. A guard can at best cut no honest trail.
- **Reference:** the same walk with no guard (`ref`), on identical seeds.

## The table

Scale model: the challenge curve and the production host client, at DP
weight 44 (mean trail `μ = 6907.3`), guard checked every 64 steps, 16,384
walks, 16 mean trail lengths per walk (`S = 110,592` steps). Three run ids:
30001 (development), 30002 and 30003 (holdouts). Figures from
[results.json](results.json).

| arm | production cap | scale cap `c` | `c / μ` | trails ended | cut | `D`, three run ids | predicted `D` (same horizon) | steady-state `D` | loss factor | checks |
|---|---|---:|---:|---:|---:|---|---:|---:|---:|---|
| `ref` | none | ∞ | ∞ | 730,879 | 0 | 0, 0, 0 | 0 | 0 | ×1.000 | pass |
| `c30` (old) | `2^30` | 20,800 | 3.011 | 760,507 | 32,631 | 14.23%, 14.22%, 14.29% | 14.21 ± 0.12% | 15.59% | ×1.166 | pass |
| `c31` | `2^31` | 41,600 | 6.023 | 732,021 | 1,195 | 1.071%, 1.036%, 1.039% | 1.076 ± 0.053% | 1.46% | ×1.011 | pass |
| `c32` (new) | `2^32` | 83,200 | 12.045 | 730,882 | 1 | 0, 0, 0.0053% | 0.0023 ± 0.0034% | 0.0071% | ×1.000 | pass |

Against the prediction the measured `D` sits at `z = 0.19, 0.08, 0.65` for
`c30` and `−0.09, −0.76, −0.69` for `c31`; the cut shares at `0.33, −0.04,
0.73` and `−0.05, −0.79, −0.65`. `c32` expects 0.43 cuts per run and saw 0,
0 and 1 (Poisson tails 1.0, 1.0, 0.35). The measured `D` is below the steady
state because a run ends after 16 trail lengths and the trail still running
at the end is not counted; the prediction simulates the same horizon, and the
steady state is what a campaign that runs for years sees.

The trail-length law is the one the fleet's interval implies. The even-weight
tail `θ = Σ_{even k ≤ w} C(131, k) / 2^130` gives `2^-28.409` at the
campaign's weight 32, where the fleet measured `2^-28.41`
(`../dp-interval/`), and `1.447744e-4` at weight 44, where the reference arms
measured `1.446494e-4 ± 0.0017e-4`, 0.09% low and `−0.7` standard errors.

## Checks

Every one holds on all twelve runs (`results.json`, `failures: []`):

1. **Exact accounting.** For every one of the 16,384 walks of every run, the
   last trail's start equals the sum of the reported trails' lengths, `c`
   per cut trail, and the idle steps each waited for its launch to end. The
   2.17 × 10^10 steps of the twelve runs are each assigned to exactly one
   reported trail, cut trail, idle stretch or trail in flight.
2. No capped arm reports a trail longer than its cap.
3. **Per seed**, against the reference on the same run id: every seed
   reported by both arms (701,581 under `c30`, 729,784 under `c31`, 730,878
   under `c32`) has identical length and orbit; of the reference
   trails longer than the cap that the capped arm finished (9,750, 9,692 and
   9,856 under `c30`; 366, 360 and 369 under `c31`; one under `c32`) none was
   reported, and no trail within the cap was cut.
4. Every record's eight witness counts sum to its length; the 48 reference
   rewalks (`--verify 4` per run) matched.
5. Every run exited 0, with no report dropped.

## Production projection

A model row, admissible because the scale model succeeded: `D(c)` at
`θ = 2^-28.41`, guard overshoot (at most 3,072 steps) and launch idle (at
most 1,024 per trail) left out as below `2^-16` of a trail. Work per solve is
the σ walk's `2^60.809 × [1.070, 1.082]` on completed trails
(WALK-CONSTANT.md §6) times the loss factor.

| `maxIters` | trail lengths | trails cut | steps discarded | loss factor | σ work per solve |
|---|---:|---:|---:|---:|---:|
| `2^30` (before) | 3.01 | 4.93% | 15.60% | ×1.1848 | `2^61.15–61.17` |
| `2^31` | 6.02 | 0.243% | 1.465% | ×1.0149 | `2^60.93–60.94` |
| `2^32` (now) | 12.04 | 0.0006% | 0.0071% | ×1.00007 | `2^60.91–60.92` |

`speedup = 1.1848 / 1.00007 = 1.185` in expected σ-walk iterations per solve,
from the model; the scale model measures the mechanism and the law it rests
on, not the campaign. No per-step cost changes: the guard is checked the same
way at any cap.

Class under AGENTS.md §3: **engineering.** The walk, its constant and the
floor are unchanged; a self-imposed loss is removed.

## What it does not cover

- The campaign runs the packed CUDA kernel; this ran the bitsliced host
  engine, which shares the kernel header's guard rule (`include/kernel.h`,
  mirrored in `include/packedkernels.cuh`). No GPU was used.
- Collectors that walk with no guard at all are unaffected by the cap: the
  Modal launcher passes no `--max-iters`, and the FPGA engine has no per-walk
  guard. Their trails can exceed any resolver's cap; see WALK-CONSTANT.md
  §11.3.
- A cairn claim is still capped by its job's `max_steps_per_walker`, which
  the production job in the cairn repository sets to `2^30`.

## Reproduce

```
./run.sh                                # about 8 minutes on four cores
WORK=/path ./run.sh analyze             # re-analyze kept raw files
```

The 24 raw files (12 corpora and 12 final checkpoints, 226 MB) are not in
the repository. A single-thread run is a function of its run id, so `run.sh`
regenerates them byte for byte; [manifest.json](manifest.json) holds their
SHA-256, the binary's and the exact compile command, and rerunning holdout
30002 from a fresh build reproduced all eight of its files' hashes. Client
logs are in [raw/](raw/). Wall time, as a practicality note: 150–154 s per
run at 11.8–12.1 M steps/s, four runs at a time.
