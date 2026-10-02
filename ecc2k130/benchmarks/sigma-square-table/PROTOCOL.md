# Sigma-fused lambda square table: frozen one-knob protocol

Protocol source parent: `18ee85c1075010ce83be60a4a02636e9aa4ac83a`.

This bounded engineering experiment asks whether the existing exact
polynomial-square table helps the selected counter-free sigma-fused walk when
it replaces only the hot `squarePolynomial131(lambdaPoly)` call.  The new
`SIGMA_SQUARE_TABLE` build flag is default-off and independent of the table
walk's `PACKED_SQUARE_TABLE` flag.  This protocol does not authorize a GPU run;
native field checks, source replay, and CUDA compile/resource evidence must
land first.

Class: compatible Pollard-walk kernel engineering.  The iteration map, field
products per update, distinguished-point rule, checkpoint representation, and
generic-group boundary do not change.  A positive result would be a same-walk
throughput gain on one RTX PRO 6000, not progress toward a cryptanalytic
exponent change or a claim that 26 B/s is reached.

## History and novelty audit

Commit `c81c310deef144fd75354dca53be61f99df9ebd2` introduced
`fillSquareTable131` and `squarePolynomialTable131`.  The table has only been
wired through `ECC_WALK_TABLE`: `twFillConsts` appends it to the table-walk
constant buffer and `twSquare` consumes it.  The later block-v3 experiment
measured `PACKED_SQUARE_TABLE=1` only together with `PACKED_INV_POLY=2` on a
different B16/T512 table-walk kernel; the five-pair median ratio was 1.006625.
That combined result does not isolate the square table and does not transfer to
the sigma-fused resource balance.

The sigma-fused star explicitly fixed `PACKED_SQUARE_TABLE=0`, described that
flag as table-walk-only, and every retained sigma-fused build/runtime record
confirmed zero.  Its `PACKED_ALU_SQUARE=1` arm was a different implementation
and lost with a 0.988597 paired geometric mean.  No retained experiment uses a
shared square table for the sigma-fused lambda square.  This one-knob route is
therefore new rather than a replay of a rejected arm.

## Fixed implementation contract

The implementation must add `SIGMA_SQUARE_TABLE`, mapped to
`ECC_SIGMA_SQUARE_TABLE`, with values zero or one and default zero.  Value one
is valid only with:

```
SIGMA_FUSED=1 WALK_TABLE=0 PACKED_POLY_STATE=1
PACKED_SQUARE_TABLE=0 PACKED_ALU_SQUARE=0
```

The first comparison does not combine this flag with ALU polynomial squaring,
the table walk, polynomial-inversion changes, or any other tuning arm.

The candidate reuses the existing 13-window, 5-bit layout exactly:

```
5 output words * 13 windows * 32 entries * 4 bytes = 8,320 bytes/block
```

For fixed output word and window, entry `e` is at a 32-word-aligned base plus
`e`.  Lanes choosing different entries hit different banks; lanes choosing the
same entry read the same address and use shared-memory broadcast.  No new table
layout is admitted.

The host fills those 2,080 words with `fillSquareTable131`, stores them through
the already ABI-compatible `WalkParams::twConsts` pointer, and frees them
through the existing engine member.  The kernel copies all 2,080 words to
dynamic shared memory and synchronizes before the partial-final-block
`tid >= p.threads` return.  The selected shared sigma network separately uses
1,792 static bytes, so the source allocation is:

| resource | bytes/block | two blocks/SM |
|---|---:|---:|
| shared sigma masks | 1,792 | 3,584 |
| lambda square table | 8,320 | 16,640 |
| combined | **10,112** | **20,224** |

The dynamic component is below CUDA's default 48 KiB limit.  The compile gate
must nevertheless show the actual ptxas registers, local bytes, static shared
bytes, dynamic launch bytes, and occupancy API result.  The candidate is
inadmissible if B16/T256 cannot retain two resident blocks per SM.

## Exact operation ledger

There is exactly one lambda square per complete scalar update in the fused
reverse pass.  With `PACKED_ALU_SQUARE=0`, the two arms perform:

| per update | control | candidate | exact delta |
|---|---:|---:|---:|
| carry-less 32-bit spreads | 5 | 2 | **-3 `CLMAD`** |
| dense 9-word-to-5-word beta reduction | 1 | 0 | **-1 reduction** |
| 5-bit window extractions | 0 | 13 | **+13** |
| 32-bit shared table reads | 0 | 65 | **+65 `LDS`** |
| table-result XOR accumulations | 0 | 65 | **+65 XOR** |

This is an exact source and memory-operation inventory for the frozen
algorithm, not a SASS instruction count or throughput prediction.  Compiler
fusion can change the ALU instruction count.  The one-time cooperative table
copy costs exactly 2,080 global reads and 2,080 shared stores per block per
kernel launch, plus one block barrier; it is launch setup rather than
per-update arithmetic and remains charged in every timed candidate sample.

## Native and compile admission gates

Before any GPU execution:

1. Re-run the existing 20,134-case square control: 131 basis vectors, three
   edge values, and 20,000 deterministic dense values.  Each table square must
   equal both `squarePolynomial131` and the independent field reference.
2. Compile the same sigma-lambda wrapper twice, once with
   `ECC_SIGMA_SQUARE_TABLE=0` and once with `=1`.  Their complete deterministic
   replay streams and SHA-256 digests must match byte for byte.
3. A separate native audit must reopen the stream, table contents, exact
   operation ledger, flag guards, and source manifest without calling the
   producer's decision code.
4. Compile both exact B16/T256/minBlocks2 CUDA clients for native `sm_120` with
   CUDA 13.3.73.  Both must have zero spills and zero local bytes.  The
   candidate must bind `squarePolynomialTable131` in the walk binary, report
   8,320 dynamic shared bytes, and retain two resident 256-thread blocks/SM.
5. Compile the existing packed arithmetic, compact-storage, and shared-sigma
   CUDA correctness probes under both arm configurations.

Failure of any gate stops the route before a GPU benchmark.  Passing these
gates authorizes only the separately dispatched matched run below.

## Frozen future GPU comparison

Both arms use one NVIDIA RTX PRO 6000 Blackwell Server Edition, native
`sm_120`, CUDA 13.3.73, batch 16, 256-thread blocks, minimum two blocks,
96,256 correctness threads, 385,024 timing threads, `WITNESS=0`, and the
selected fused-sigma settings:

```
PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1
PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1
PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1
PACKED_GENERATED_PRODUCT=1 PACKED_INLINE_POLY=3 PACKED_CLMAD=1
PACKED_CLMAD_SQUARE=0 PACKED_KARAT3=0 PACKED_WEIGHTED_PREFIX=2
PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1 PACKED_TOP_CLMAD=0
PACKED_STATE_TILE=256 PACKED_ADD_COMBINE=0 PACKED_PAIR_ILP=0
UNROLL_SLOTS=1 PACKED_SLOT_PREFETCH=0 PACKED_SLOT_PIPELINE=0
SIGMA_FUSED=1 SIGMA_FUSED_LATE_Y=0 PACKED_L2_PERSIST=0
PACKED_TOP_HOIST=0 PACKED_ONB_INV=0 PACKED_FROM_REDUCED=0
PACKED_CLMUL_FLAT=0 PACKED_PAIR_CLMUL=0 PACKED_ALU_SQR=0
PACKED_ALU_SQUARE=0 PACKED_SQUARE_TABLE=0 PACKED_INV_POLY=0
PACKED_CHAINS=1 WALK_TABLE=0 WITNESS=0 PHASE_PROFILE=0
```

Control sets `SIGMA_SQUARE_TABLE=0`; candidate sets only
`SIGMA_SQUARE_TABLE=1`.

Before timing, both binaries run seven odd 95-step launches at DP weight 48
and run id 43.  Each must replay 300/300 reports with zero drops and produce
the same sorted headerless-v1 corpus.  Bidirectional whole/prefix checkpoint
replay must also be byte-identical.  A build, marker, resource, replay, corpus,
or checkpoint mismatch suppresses timing.

Each timing row completes exactly:

```
385024 threads * 16 slots * 1024 steps * 32 launches
= 201863462912 complete scalar updates
```

Two control/candidate warmups are excluded.  Five A/B pairs alternate order,
interleaved with five A/A pairs of byte-identical control binaries.

## Decision and stop rule

Promote the candidate as compatible engineering only when all correctness and
resource gates pass, every A/B ratio is above 1.0, the paired median is at
least 1.005, every A/A ratio lies in `[0.995, 1.005]`, and the A/A median lies
in `[0.998, 1.002]`.

Reject and leave the flag off if correctness/resource identity fails or the
A/B median is at most 1.0.  A positive median below 1.005 is retained as
inconclusive/optional.  The 26 B/s objective is reported separately and is met
only by a measured complete-walk candidate rate at or above 26 B/s; no source
ledger or isolated primitive rate can satisfy it.
