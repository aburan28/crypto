# GPU-wide resolver geometry: frozen matched screen

Source parent: `cf9ce34d711b145b4b4031a3a8bbaf7e6df15adb` on
`codex/ecc2k130-resolver-threads-20261002`.

Status: authorized bounded benchmark.  Dispatch is forbidden until the native
producer is committed, independent source review passes, and the three device
queue controls are present in the job.  Exactly one RTX PRO 6000 allocation is
allowed.  Owner decoding and every non-geometry input remain fixed.

The parent GPU-wide experiment is terminal `DO_NOT_PROMOTE`: block control
7.021655 B/s versus resolver-128 5.610974 B/s, geometric-mean ratio
0.79901828425.  This follow-up tests whether cold-resolver occupancy can recover
some of that loss.  It does not reinterpret the parent result or lower its
1.10 promotion threshold.

## Arms and invariant work

One NVIDIA RTX PRO 6000 Blackwell Server Edition, CUDA 13.3.73 and native
`sm_120` code.  All arms use the exact parent B16/T512/min1 table-v3 preset,
385,024 timing workers, 16 slots, 1,024 steps and 32 launches.  Every timed row
therefore completes exactly **201,863,462,912 scalar updates**.

| arm | hint schedule | resolver blocks | resolver threads | changed input |
|---|---|---:|---:|---|
| `control` | original block queue | n/a | n/a | reference |
| `r128` | GPU-wide queue | 188 | 128 | parent candidate/reference geometry |
| `r256` | GPU-wide queue | 188 | 256 | resolver width only |
| `r512` | GPU-wide queue | 188 | 512 | resolver width only |

The queue has 6,160,384 entries (24,641,536 bytes) and a four-byte counter at
the timing population.  Each resolver block copies the same 57,052-byte table.
The grid, exact-v3 tags, owner keys, `owner / workers`, `owner % workers`, cold
anchor, queue reset, stream ordering, select kernel, hot kernel, arithmetic,
state layout, reports and checkpoints do not change.  `64` threads and owner-
decode changes are out of scope.

## Pre-timing gates

Timing is suppressed unless all gates pass:

1. native ownership/phase model, table-v3 replay and packed arithmetic controls;
2. native `sm_120` builds for the block control and all three resolver widths;
3. one production device queue control per width, each reporting exactly 49
   empty/full/sparse/endpoint/repeated-reset cases, intact canaries, exact
   scalar-reference histories, its requested launch maximum, exactly one
   active block/SM and 57,052 dynamic shared bytes, preserving the registered
   4/8/16-warps-per-SM interpretation;
4. runtime resource receipts for hot, select and resolver stages, including
   registers, local bytes, static/dynamic shared bytes, launch bounds, active
   blocks/SM, 188 resolver blocks and exact queue/counter bytes;
5. seven odd 95-step launches at DP weight 48 with 300/300 host replay and zero
   drops for every arm at 96,256 workers;
6. identical complete sorted v3 corpora across all arms at full and 511/513
   partial populations;
7. deterministic 300-record spread replay across all full corpora with at least
   299 nonzero trails; and
8. four-launch prefixes and three-launch continuation for every prefix/arm
   combination at 513 workers, with byte-identical final checkpoints and
   sorted continuation corpora.

Kernel-resource and resolver-geometry lines captured by the first full
verification become the immutable per-arm resource receipts; every subsequent
log must reproduce them.  Queue entries and bytes are checked exactly against
each log's worker population.  Any mismatch, overflow, missing owner,
duplicate owner, stale queue entry, report mismatch, dropped record, partial
record or incorrect work count invalidates the run.

## Timing order

Four full-work warmups are excluded, one per arm.  Warmups never enter a
selection or ratio.

Five A/A pairs then run the block control twice, with odd pairs `A,B` and even
pairs `B,A`.  Maximum symmetric drift must stay below 1%.

Five ordered screen rounds follow:

```
round 1: control, r128, r256, r512
round 2: r512, r256, r128, control
round 3: r128, control, r512, r256
round 4: r256, r512, control, r128
round 5: control, r128, r256, r512
```

Each arm occupies every order position at least once and one position twice.
Ratios pair arms within the same round.  The native summarizer rejects missing,
duplicated or reordered rows and independently reopens all logs and exact work.

## Decisions

- The geometry statistic is the median of the five within-round `r256/r128` or
  `r512/r128` ratios, never a ratio of medians or geometric mean.  A non-default
  width qualifies as a **geometry-only** follow-up only when that paired median
  is at least 1.05, all five ratios exceed one, A/A drift is below 1%, and every
  gate passes.  If both qualify, the width with the higher paired median is the
  geometry candidate; an exact tie selects 256 deterministically.  Qualification
  does not promote either width in this screen.
- The unchanged parent GPU-wide-map gate is evaluated only on the preregistered
  `r128/control` pair: geometric mean of its five within-round ratios at least
  1.10, every ratio above one, the same noise/correctness gates, and no omitted
  cold cost.  No post-hoc best-of-three width may satisfy the parent gate.
- The 26 B/s objective is met only by a verified complete-update median above
  26,000 M/s on this one GPU.
- Otherwise retain resolver width 128 and the block-queue production path, and
  preserve the negative result.

The experiment stops after this panel.  No search, solver, collision recovery,
key recovery, CUDA graph, owner-decode optimization or second geometry run is
part of it.
