# Fused sigma schedule: frozen RTX PRO 6000 comparison

Baseline commit: `efac7fc1f1a0eac8271313b3e9bea8194fe78f01`.

Class: compatible kernel engineering. This experiment keeps the ECC2K-130
sigma walk, branch rule, polynomial state, batch 16, 256-thread blocks,
minimum two resident blocks, worker population, witness format and update
count fixed. It neither runs a challenge search nor changes the collision
model.

## Hypothesis

The current kernel loads each `(x,y)` state in the forward pass, then again in
the reverse pass. The candidate begins a launch with the same forward pass,
but after every reverse update except the last it immediately computes the next
step's selection, denominator, numerator and prefix contribution from the new
point still in registers. Prefix order alternates between ascending and
descending slots so Montgomery's inversion is consumed in reverse order.

This removes one coordinate load pass per completed update. The arithmetic and
stored point after every launch remain identical. The candidate is selected by
`SIGMA_FUSED=1`; the default is zero.

## Frozen builds and hardware

Both arms use native `sm_120` code from CUDA 13.3.1 and this knob set:

```
BATCH=16 THREADS=256 MINBLOCKS=2
PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1
PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1
PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1
PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=1 PACKED_STATE_TILE=256
PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1
PACKED_INLINE_POLY=3 WITNESS=1 WALK_TABLE=0
```

The comparison runs on one NVIDIA RTX PRO 6000 Blackwell Server Edition. Rate
rows use 385,024 workers, 1,024 steps and 64 launches, so both arms complete
the same scalar-update count.

## Correctness and timing protocol

Before timing:

1. Build both binaries and retain compiler/resource logs and binary hashes.
2. Run each with an odd 95-step launch length, seven launches, DP weight 48,
   and host-reference verification of 300 reports.
3. Require 300/300, zero dropped, no mismatch or overflow, the requested fused
   marker, and byte-identical sorted complete v2 corpus records from a native
   comparator. Any failure suppresses timing.

Timing has four excluded warmups, two per binary. It then runs five alternating
64-launch control/control pairs to measure session noise and five alternating
control/candidate pairs. Pair order reverses each round. A native summarizer
reopens the TSV and computes all ratios.

## Decision

- Goal met: candidate median exceeds 26,000 M complete scalar updates/s, every
  candidate/control pair exceeds one, and all correctness gates pass.
- Promote as engineering: A/B median ratio is at least 1.02, every A/B pair
  exceeds one, and maximum absolute A/A pair drift is below 1%.
- Otherwise retain the experiment without selecting the knob. A result below
  26 B/s does not meet the active goal even if it is promoted.

The experiment stops after this matched panel. No search, solver, collision
recovery or multi-GPU aggregation is part of it.

