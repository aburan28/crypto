# Iterations per distinguished point at the running cutoff

The campaign distinguishes at `HW(x) <= 32` (`dpWeight` 32 in the live
`campaign.json`). The interval between reports at that cutoff had been quoted
as `2^27.9` iterations per point, from the 2026-09-15 ratio of two fleet-wide
counters: `7.93e15` checkpointed iterations against `31.6 M` uploaded points.
Neither counter was clean. The checkpoint sum multiplied every slot's per-walk
base by one campaign constant (`385,024 x 16`), which is not the grid an Ada
slot runs (`aws/README.md`, monitoring); the point total mixed in the first
two slots of the campaign, which walked at weight 34 under the CUDA 13.0
build (their retired state still reads `2^25.2` iterations per point, the
weight-34 interval). A ratio of a mis-weighted sum to a mixed count is not a
measurement of the running cutoff.

## Method

Every worker logs, once a minute, the client's own two counters for the
current process: iterations this run, computed by the client as
`(iterBase - timedIterBase) x walksPerLaunch` on the grid it built itself, and
distinguished points reported over the same interval (`aws/worker.py`,
`src/main.cu`). No campaign constant and no cross-slot sum enters, so the
ratio is exact for that worker's cutoff whatever its grid. The last such line
of every live worker was read from the log copies in the bucket at
2026-09-17 18:47 UTC and is frozen in
[`raw/worker-counters-2026-09-17.tsv`](raw/worker-counters-2026-09-17.tsv).

## Result

| family | workers | iterations this run | points | iterations per point | per-worker log2, min / median / max |
|---|---:|---:|---:|---:|---|
| RTX PRO 6000 (g7e) | 38 | 5.606e15 | 15,724,672 | 2^28.41 | 28.41 / 28.41 / 28.42 |
| L40S (g6e) | 34 | 2.723e15 | 7,643,600 | 2^28.41 | 28.40 / 28.41 / 28.41 |
| RTX PRO 4500 (g7) | 48 | 3.081e15 | 8,640,973 | 2^28.41 | 28.40 / 28.41 / 28.43 |
| CPU | 12 | 2.804e14 | 785,410 | 2^28.41 | 28.40 / 28.41 / 28.44 |
| all | 132 | 1.169e16 | 32,794,655 | **2^28.41** | |

One report per `2^28.41 = 3.56e8` iterations, the same on every GPU family
and on the CPU workers, with a per-worker spread of `0.04` in the exponent
over 132 workers. The L40S slots are the ones that auto-size their grid; they
agree with the Blackwell slots to two decimals, which is what a grid-free
measurement must do and what the 2026-09-15 ratio could not.

Cross-check against the weight-34 figure the plan is priced at: the plain
binomial tail `P(HW <= w)` over 131 bits gives `2^-29.01` at `w = 32` and
`2^-25.84` at `w = 34`; the measured `2^28.41` sits `2^0.60` below the first
and Bailey et al.'s `2^25.27` sits `2^0.57` below the second. The two
figures carry the same offset from the naive count, so the running cutoff
is `2^3.14 = 8.8x` sparser than the weight-34 rows in `aws/README.md`, not
the `2^2.6 = 6x` the retired figure implied.

## What this changes

Class: **accounting**, per AGENTS.md §3. The walk did not move; a figure
did. The superseded value stays visible as the before mark:

| figure | before (2026-09-15, retired) | now (2026-09-17, measured) |
|---|---:|---:|
| iterations per point at weight 32 | 2^27.9 | 2^28.41 |
| points per GPU-day at 14.1 B it/s | ~5 M | ~3.4 M |
| points per GPU | ~56/s | ~40/s, ~0.11 GB/day |
| corpus at the expected work 2^60.9 | 2^33.0 records | 2^32.5 records, ~0.19 TB |
| walk in flight per GPU at the collision (6.16 M walks x interval) | 30 GPU-h | 43 GPU-h |
| weight-34 -> weight-32 rarity factor | 2^2.6 | 2^3.14 |

The consumers of the constant are `aws/README.md`, `scripts/rho_status`
(`work_feed.py`, its README), the public dashboard
`docs/ecc2k130-status/index.html`, and `scripts/site/test_build.py`, which
pins every copy to the campaign document. All of them carry `2^28.4` now.

Not established here, and left for the merge step: the corpus in `dp/` holds
the weight-34 points of the first two slots next to the weight-32 points of
every later one. Weight 34 is the looser cutoff, so a weight-32 walk that
merges into one of those early trails stops at a later point than the early
walk did and the collision is seen only if the early walk's stopping point
also has `HW(x) <= 32`. `aws/README.md` already says a `dpWeight` change
breaks the collision guarantee; this is the concrete form of it in this
bucket.
