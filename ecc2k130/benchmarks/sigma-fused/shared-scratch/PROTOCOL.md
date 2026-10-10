# Fused-sigma shared scratch: frozen follow-up protocol

Status: **frozen before implementation; no GPU run is authorized by this
document.**  Parent source is current `main`
`18ee85c1075010ce83be60a4a02636e9aa4ac83a`.  The static source and traffic
audit is [`PROPOSAL.md`](PROPOSAL.md), and its native receipt is
[`model.json`](model.json).

## Scope and hypothesis

This is compatible kernel engineering for the ECC2K-130 sigma walk.  It does
not change the iteration function, scalar walks, automorphism quotient,
collision model, report format, checkpoint format or generic work.  No search,
collision collection, solver or key recovery is in scope.

The hypothesis is that routing the first `C` per-thread `pchain` and tagged
denominator slots from compact global storage to a full-word shared SoA can
reduce the selected fused kernel's hot global field traffic enough to overcome
the added shared-memory instructions and reduced L1 capacity.  Static logical
savings are 8.5, 12.75 and 17 bytes/update for `C = 2,3,4`.  They are not a
bandwidth measurement or a predicted speedup.

## Frozen arms

All arms use one NVIDIA RTX PRO 6000 Blackwell Server Edition, CUDA 13.3.73,
native `sm_120`, 385,024 workers, batch 16, 256-thread blocks and a minimum of
two resident blocks.  Every non-cache input is the selected counter-free fused
sigma preset:

```
BATCH=16 THREADS=256 MINBLOCKS=2 WITNESS=0 WALK_TABLE=0
PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1
PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1
PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1
PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=1 PACKED_STATE_TILE=256
PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1
PACKED_SHARED_SIGMA=1 PACKED_INLINE_POLY=3 PACKED_CHAINS=1
SIGMA_FUSED=1 SIGMA_FUSED_LATE_Y=0
```

| arm | `SIGMA_FUSED_SHARED_SLOTS` | registered function shared bytes | role |
|---|---:|---:|---|
| `control` | 0 | 1,792 | selected reference |
| `cache2` | 2 | 22,272 | low-footprint candidate |
| `cache3` | 3 | 32,512 | middle candidate |
| `cache4` | 4 | 42,752 | maximum two-block full-word candidate |

Values 1 and 5--16 are outside the experiment.  Five slots are the registered
capacity stop because two compiled extents would require 108,032 bytes against
the 102,400-byte SM ledger.  Compact shared packing, caching `x/y`, changing
the cached slot positions, resizing global allocations and combining another
knob are separate hypotheses.

## Implementation and compile gate

Timing is forbidden until a focused implementation PR:

1. exposes only the closed `0/2/3/4` knob and rejects every incompatible
   configuration at compile time;
2. uses the `[field][slot][word][thread]` SoA from the static model, routes
   both producer and consumer through one helper contract, and adds no
   inter-thread barrier;
3. preserves all arithmetic, alternating slot order, first/final fused
   boundaries, global allocations and checkpoint version;
4. passes the native model and a native host ownership/boundary model;
5. builds all arms from one clean exact source revision with CUDA 13.3.73 and
   retains source, command, compiler and binary hashes; and
6. records ptxas plus runtime attributes for every arm.

Every arm must report its registered function shared bytes, zero local bytes,
zero stack bytes, zero spill loads and stores, no more than 128 registers per
thread, launch bounds `256 x 2`, and exactly two active blocks/SM from the CUDA
occupancy API.  A missing or different value is a producer failure.  No subset
of arms proceeds adaptively.

## Device correctness gate

Before performance timing, one bounded RTX PRO 6000 allocation must pass all
of the following for all four arms:

1. the existing packed arithmetic, compact-storage and shared-sigma device
   suites;
2. a scratch-specific device test over basis, dense and tagged-top-word values,
   including partial blocks and canaries, comparing shared helpers with the
   compact global helpers;
3. seven odd 95-step launches at DP weight 48, with 300/300 reports replayed by
   the host reference, zero drops and no mismatch or overflow;
4. complete sorted corpus identity against the control;
5. one-, two-, three-, seven-, sixteen- and 95-step launch-boundary state
   comparisons; and
6. every control/candidate prefix and continuation combination at a partial
   worker count, with byte-identical final checkpoints and sorted continuation
   corpora.

Failure suppresses every timing row and is preserved as a producer failure.

## Equal-work timing

Each timing sample uses 385,024 workers, 16 slots, 1,024 steps and 32 launches:

```
385024 * 16 * 1024 * 32 = 201863462912 complete scalar updates
```

Four full-work warmups, one per arm, are excluded.  Five alternating A/A pairs
then run the control twice to measure session noise.  Five four-arm rounds use
this fixed order:

```
round 1: control, cache2, cache3, cache4
round 2: cache4, cache3, cache2, control
round 3: cache2, control, cache4, cache3
round 4: cache3, cache4, control, cache2
round 5: control, cache2, cache3, cache4
```

A native summarizer must reject missing, duplicate, reordered or unequal-work
rows and reopen every log.  Candidate/control ratios are paired within a
round.  Warmups never enter a statistic.

## Decisions and stop

- Noise is acceptable only when maximum A/A symmetric drift is below 1%.
- A candidate qualifies for a separate held-out confirmation only when its
  paired geometric mean over control is at least 1.01 and all five paired
  ratios exceed one.
- If multiple candidates qualify, the highest paired geometric mean wins; an
  exact tie selects the smaller cache to retain more unified L1 capacity.
- The 26 B/s objective is met only if a correctness-qualified candidate's
  five-sample median exceeds 26,000 M complete scalar updates/s.  Static bytes
  and an isolated stage rate cannot satisfy this gate.
- Otherwise retain the current fused preset and preserve the negative result.

The experiment stops after this panel.  There is no automatic second cache
layout, compact-shared arm, cache-position sweep, launch-geometry change or
second GPU allocation.  A qualifier is not promoted by this screen; its next
gate is a separately frozen held-out confirmation.
