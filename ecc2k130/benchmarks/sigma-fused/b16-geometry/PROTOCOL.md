# Fused sigma B16 block geometry: frozen five-pair protocol

Source parent: `60fb2ae47fcca9352d8bae2b280fce77371e57f0` (`origin/main`).

Class: compatible kernel geometry.  The earlier fused geometry screen retained
B16/T256/min2 while testing T512 only at batches 32 and 64.  A later,
separately allocated diagnostic produced a faster B16/T512 control, but a
cross-allocation rate is not paired evidence.  This experiment closes the
missing B16 comparison on one allocation.

This is benchmark-only defensive research.  It does not search for challenge
points, collect a campaign corpus, solve collisions, run a solver, recover a
key, or aggregate multiple GPUs.

## Frozen arms

Both arms use one NVIDIA RTX PRO 6000 Blackwell Server Edition, CUDA 13.3.73,
native `sm_120`, batch 16, 385,024 worker threads and 6,160,384 live walks.

| arm | block threads | launch bound | expected occupancy | resident warps/SM |
|---|---:|---:|---:|---:|
| `t256` | 256 | minimum 2 blocks/SM | 2 blocks/SM | 16 |
| `t512` | 512 | minimum 1 block/SM | 1 block/SM | 16 |

All arithmetic, walk and reporting knobs are identical:

```
PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1
PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1
PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1
PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=1 PACKED_STATE_TILE=256
PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1
PACKED_INLINE_POLY=3 SIGMA_FUSED=1 SIGMA_FUSED_LATE_Y=0
WITNESS=0 WALK_TABLE=0
```

The rejected direct-table implementation is absent from current main and is
therefore fixed off in both arms, equivalent to `DIRECT_SIGMA=0`.

## Pre-timing gates

The job fails closed on the exact GPU, compiler, architecture, source commit,
build flags, binary hashes, worker counts and launch bounds.  Before timing:

1. both arms pass the packed CUDA arithmetic, compact-storage and shared-sigma
   tests;
2. ptxas reports zero stack, zero spill stores and zero spill loads for the
   exact `eccPacked131::walk` symbol;
3. runtime reports 126 registers/thread, zero local bytes/thread and 1,792
   static shared bytes/block for each arm, with the expected occupancy above;
4. each arm runs seven odd 95-step launches at DP weight 48, replays 300/300
   reports against the host reference with zero drops, and produces the same
   complete sorted v1 corpus; and
5. uninterrupted checkpoints and both cross-arm resume directions finish at
   byte-identical state.

Any failure suppresses all timing.  The native post-run checker reopens every
sample log and validates its SHA-256, exact build identity, backend population,
completed update count, zero drops and terminal rate.

## Timing and decision

Every timing row runs 1,024 steps and 64 launches, exactly
`385024 * 16 * 1024 * 64 = 403,726,925,824` complete scalar updates.  Four
warmups are excluded, two per arm.  Five alternating `t256/t256` pairs measure
session noise, followed by five alternating `t256/t512` pairs.  Pair order is
AB, BA, AB, BA, AB in each panel.

- Select T512 only if median `t512/t256 >= 1.01`, every A/B pair is above one,
  A/A maximum absolute drift is below 1%, and all correctness/equal-work gates
  pass.
- Record goal met only if the selected median exceeds 26,000 M complete scalar
  updates/s under the same gates.
- Otherwise retain T256 and preserve the result.

No separate-allocation rate enters the decision.
