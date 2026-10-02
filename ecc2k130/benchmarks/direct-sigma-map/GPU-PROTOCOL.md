# Direct sigma shared-table screen: frozen RTX PRO 6000 protocol

Source parent: `8ce1b207` on
`codex/ecc2k130-direct-sigma-gpu-20261001`.

Class: compatible kernel engineering.  The direct table computes the same
`toPolynomial((I + sigma^j) fromPolynomial(p))` map selected by the existing
walk.  This experiment runs validation and benchmark modes only.  It does not
search for challenge points, solve collisions, or recover a key.

## Frozen candidate

The native static screen admitted the three-input-bit table: 56,320 bytes for
all eight maps, all 1,048 basis cases exact, and a 30.2% median host static
instruction reduction against the reduced-input composition.

For the GPU screen, the same bits are laid out in dynamic shared memory as
`[chunk][output-group][(j-3)*8 + value]`.  Each lookup reads the low 128 bits
as one aligned vector and the top three bits as one word.  A block copies the
table once before processing.  Constant memory is excluded because unrelated
walks make the `(j,value)` address divergent.  Global lookup memory is
excluded by the static protocol.

The candidate replaces only the two `sigma^j + 1` maps used to form X and Y
differences.  X is still converted to normal basis for the unchanged Hamming
weight and branch rule.  All field products, inversion, point formulas,
distinguished-point test, restart behavior and scalar-update accounting stay
fixed.

## Builds and equal work

Both arms use CUDA 13.3.73, native `sm_120`, one NVIDIA RTX PRO 6000 Blackwell
Server Edition, batch 16, 512-thread blocks, minimum one resident block and
385,024 worker threads (6,160,384 live walks).  This retains 16 resident warps
per SM in both arms.  Common knobs are:

```
PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1
PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1
PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1
PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=1 PACKED_STATE_TILE=256
PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 PACKED_INLINE_POLY=3
WITNESS=0 WALK_TABLE=0 SIGMA_FUSED=1
```

Control has `DIRECT_SIGMA=0 SHARED_SIGMA=1`; candidate has
`DIRECT_SIGMA=1 SHARED_SIGMA=0`.  Candidate dynamic shared memory is exactly
56,320 bytes.  Control uses its existing small shared sigma masks.  The job
must record compiler versions, binary hashes, register/local/shared resources,
the requested knobs, worker counts and completed scalar counts.

## Correctness gate

Before timing:

1. regenerate the table with the native C++ generator and require the frozen
   table SHA-256;
2. run the existing packed CUDA arithmetic and storage tests for both arms;
3. run a dedicated GPU kernel over every 131 basis input for all eight powers
   plus deterministic dense inputs, comparing direct and composed outputs;
4. run seven odd 95-step launches at DP weight 48 for each arm, replay 300
   reports against the host reference, require zero drops, and compare complete
   sorted headerless-v1 corpora byte for byte; and
5. exchange a checkpoint between arms in both directions and require the same
   completed state after a common suffix.

Any mismatch, requested-knob mismatch, launch failure, spill, incomplete
sample or corpus difference suppresses performance conclusions.

## Timing and decision

Every sample uses 1,024 steps and 64 launches, exactly
403,726,925,824 complete scalar updates.  Four warmups are excluded, two per
binary.  Three alternating control/control pairs establish session noise,
then five alternating control/candidate pairs reverse order each round.

- Goal met: candidate median exceeds 26,000 M complete scalar updates/s, every
  A/B pair exceeds one, all correctness gates pass and A/A maximum absolute
  drift is below 1%.
- Qualify for confirmation: A/B median is at least 1.02, every A/B pair exceeds
  one, correctness passes and A/A drift is below 1%.
- Otherwise retain the static circuit and measured rejection without selecting
  the knob.

The 30% host instruction cut is a screening fact, not a predicted GPU gain.
The extra shared-memory traffic and one-block geometry can erase it.
