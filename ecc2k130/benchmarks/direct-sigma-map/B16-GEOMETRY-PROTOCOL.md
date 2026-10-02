# Fused sigma B16 block geometry: frozen RTX PRO 6000 comparison

Source parent: `85a6c793` on
`codex/ecc2k130-direct-sigma-gpu-20261001`.

Class: compatible kernel geometry.  The direct-map experiment's matched
control measured 15.893066 B/s at B16/T512/min1, while the separately
allocated published fused headline is 15.436677 B/s at B16/T256/min2.  The
prior fused geometry panel retained B16/T256 as its baseline but tested T512
only at batches 32 and 64.  Cross-allocation rates cannot select a geometry,
so this protocol compares the missing B16 pair on one allocation.

No search, collision recovery, solver, key recovery or multi-GPU aggregation
is run.  `DIRECT_SIGMA=0` in every arm.

## Builds and equal work

Both arms use CUDA 13.3.73, native `sm_120`, one NVIDIA RTX PRO 6000 Blackwell
Server Edition, batch 16, 385,024 worker threads and 6,160,384 live walks.

| arm | block threads | minimum blocks/SM | expected resident blocks/SM | resident warps/SM |
|---|---:|---:|---:|---:|
| t256 | 256 | 2 | 2 | 16 |
| t512 | 512 | 1 | 1 | 16 |

Every other build input is fixed:

```
PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1
PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1
PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1
PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=1 PACKED_STATE_TILE=256
PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1
PACKED_INLINE_POLY=3 WITNESS=0 WALK_TABLE=0 SIGMA_FUSED=1
DIRECT_SIGMA=0
```

The job records exact build commands, binary hashes, ptxas walk resources,
runtime occupancy, compiler/GPU identity and completed scalar counts.

## Correctness gate

Before timing, both arms must:

1. pass the packed CUDA arithmetic, compact-storage and shared-sigma tests;
2. compile the exact packed walk with zero stack and zero spills;
3. run seven odd 95-step launches at DP weight 48, replay 300/300 reports
   with zero drops and produce byte-identical sorted complete v1 corpora; and
4. exchange a checkpoint in both directions and reach the same byte-identical
   completed state as uninterrupted execution.

The t256 verification run must report two resident 256-thread blocks per SM;
t512 must report one resident 512-thread block.  Any gate failure suppresses
timing.

## Timing and decision

Each timing sample runs 1,024 steps and 64 launches: exactly
403,726,925,824 complete scalar updates.  Four warmups are excluded, two per
binary.  Three alternating t256/t256 pairs measure session noise, followed by
five alternating t256/t512 pairs.

- Select t512 geometry if median t512/t256 is at least 1.01, every pair is
  above one, A/A maximum absolute drift is below 1%, and every correctness and
  equal-work gate passes.
- Record goal met only if the selected median exceeds 26,000 M complete scalar
  updates/s under the same gates.
- Otherwise retain t256 and preserve the negative result.

This panel selects block geometry only.  It does not promote the rejected
direct table or reinterpret cross-allocation rates as paired evidence.
