# Fused sigma headline composition: frozen RTX PRO 6000 comparison

Source parent: `efac7fc1f1a0eac8271313b3e9bea8194fe78f01`.

The first fused-schedule panel prices the campaign witness counters. The active
26 B/s objective and the current 14.8547595 B/s sigma benchmark use the
counter-free path. This follow-up composes only the same fused schedule with
`WITNESS=0`; every other walk, arithmetic, geometry and hardware input stays
fixed.

Class: compatible kernel engineering. No search, solver, collision recovery or
multi-GPU aggregation is run.

## Frozen builds

Both arms use CUDA 13.3.73, native `sm_120`, one NVIDIA RTX PRO 6000 Blackwell
Server Edition, batch 16, 256-thread blocks, minimum two resident blocks and
385,024 workers. The common arithmetic is:

```
PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1
PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1
PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1
PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=1 PACKED_STATE_TILE=256
PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1
PACKED_INLINE_POLY=3 WITNESS=0 WALK_TABLE=0
```

Control has `SIGMA_FUSED=0`; candidate has `SIGMA_FUSED=1`.

## Correctness and measurement

The job fails closed on the exact GPU, `sm_120`, compiler, worker counts,
DP weight and step count. Each arm first runs seven odd 95-step launches at DP
weight 48, replays 300 reports against the host reference with zero drops, and
must produce an identical sorted headerless-v1 corpus. Timing is suppressed on
any failure.

Four warmups are excluded. Five alternating control/control pairs establish
the session noise, then five alternating control/candidate pairs each complete
the same `385024 * 16 * 1024 * 64` scalar updates. The native summarizer reopens
the TSV.

## Decision

- Goal met only if the candidate median exceeds 26,000 M complete scalar
  updates/s, every A/B pair exceeds one, correctness passes and A/A maximum
  absolute drift stays below 1%.
- Promote as compatible engineering if the A/B median is at least 1.02 with
  every pair above one under the same correctness/noise gates.
- Otherwise retain the result and leave the knob experimental.

