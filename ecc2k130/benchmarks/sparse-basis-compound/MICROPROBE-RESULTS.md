# Exact-layout sparse-table multicast microprobe: measured result

Decision: **the three-bit sparse-basis lookup route fails its prerequisite; do
not build the whole-walk candidate.**  The exact table composition supports
only **11.448092 billion complete-update equivalents/s** before charging any
field multiplication, reduction, state traffic, loop/control work, report
logic, or initialization.  That is 0.741616 times the confirmed 15.436677 B/s
fused sigma reference and below the frozen 16.208511 B/s admission threshold.

This is primitive timing, not sparse-walk throughput.  One RTX PRO 6000 was
used only for the bounded table microprobe.  No sparse walk, search,
distinguished-point collection, collision recovery, or key recovery ran.

## Single table

Every timing sample applies one complete 44-chunk map to 6,308,233,216 lane
inputs.  Each application executes exactly 44 `LDS.128` and 44 scalar `LDS`
instructions.  Medians are over five ranked samples after one excluded warmup
per arm.

| arm | key pattern | median B map applications/s | role |
|:--|:--|--:|:--|
| `joint_varying` | lane- and round-varying `(j,value)` over 64 entries | **28.638983** | two applications/update; binding route term |
| `fixed_normal` | varying value over one eight-entry map | **64.225526** | one application/update |
| `fixed_to_beta` | varying value over one eight-entry map | **64.237372** | one per batch |
| `fixed_from_beta` | varying value over one eight-entry map | **64.245682** | one per batch |
| `joint_fixed_j` | joint layout, fixed `j`, varying value | 64.250690 | matched bank-conflict control |
| `joint_uniform_key` | one joint key broadcast within the warp | 82.976896 | multicast ceiling control |
| `fixed_normal_clone` | exact A/A clone | 64.239253 | noise control |
| composed lookup-only route | `2/joint + normal + (to+from)/16` | **11.448092 B updates/s** | prerequisite fails |
| confirmed fused reference | complete scalar update | **15.436677 B updates/s** | boundary |
| frozen admission | `1.05 * reference` | **16.208511 B updates/s** | boundary |

The joint varying key is only 0.445738 times the same joint table with fixed
`j`, and the fixed-j control is 2.24347 times faster.  The measured lookup time
per update is 0.0873508 ns: 79.95% comes from the two joint maps, 17.82% from
the normal conversion, and 2.23% from the two batch-boundary conversions.
This directly identifies joint-key bank conflicts as the dominant lookup
wall.  The earlier transferred random-`LDS.U8` estimate of 15.283375 B/s was
optimistic for the complete route.

The five A/A ratios are
`1.000376/0.999700/1.000008/1.000554/1.000166`; maximum absolute drift is
0.0554%, well inside the frozen 1% gate.  Every ranked arm reproduced its
complete-output SHA-256 in all five samples.

## Exact layout and resources

The runtime constructed the actual selected pentanomial tables natively:

- eight joint `L_j` maps;
- sparse-to-normal;
- sparse-to-beta; and
- beta-to-sparse.

The table payload is 77,440 bytes, SHA-256
`0d8662afe5b454c028c6fb38125fadaf1ab836a22414e27b6babe3fc7ec9d8f2`.
With the 1,024-byte driver reservation it occupies 78,464 bytes/block.  Every
kernel residents exactly one 512-thread block/SM and reports zero local bytes.
The varying-j arm uses 121 registers/thread, below the 128-register limit for
that geometry; fixed maps use 110 and the uniform control 98.

The retained `sm_120` cubin has, in each of all seven benchmark kernels:

- 44 `LDS.128` instructions;
- 44 scalar `LDS` instructions;
- zero local load/store instructions; and
- zero calls.

All 64 joint keys were observed by the GPU coverage kernel.  The GPU checked
16,908 one-map outputs: four map families over all 131 basis inputs and 4,096
dense inputs.  One spread-selected thread in every timed arm was then replayed
natively through all 65,536 rounds.  All gates passed before ranked timing.

## Sparse ALU-square accounting

The frozen sparse route has 38.125 `CLMAD` lane-instructions/update.  Moving
the five-word lambda square to the ALU gives:

```
38.125 - 5 = 33.125 CLMAD/update
```

At the ideal two-lane rate and maximum recorded clock, that isolated ceiling
is 27.582792 B/s, 1.0609 times the 26 B/s objective.  The source circuit adds
75 ALU operations/update, moving the charged sparse reducer/map/spread
subledger to 1,542.25 operations/update versus 2,276.25 for the reference
subledger.  These are source counts, not SASS or measured walk speed.

The exact lookup wall is already 11.448092 B/s.  Therefore the favourable
carryless-only ceiling cannot rescue this table route, and a sparse
`PACKED_ALU_SQUARE=1` whole-walk build is not justified.  This conclusion does
not generalize to a sparse basis with a materially different linear circuit.

## Preserved failed producers

Three pre-result failures are retained additively:

1. `microprobe-attempt-1-failure.json`: runtime `blockDim.x` in table
   initialization emitted an out-of-line division call; the SASS gate stopped
   before any GPU kernel ran.
2. `microprobe-attempt-2-failure.json`: the call was removed and all kernels
   had the exact 44+44 loads, but the native SASS parser split repeated
   nvdisasm metadata into false zero-length functions; no GPU kernel ran.
3. `microprobe-attempt-3-failure.json`: SASS and 16,908 GPU correctness cases
   passed, but the compiler deleted repeated uniform-control maps whose prior
   outputs were dead.  Five excluded warmups were recorded; no ranked row or
   decision was produced.  The final source carries every mapped value through
   a stored dependency accumulator in all arms and in the native replay.

None of these attempts is timing evidence for the decision.

## Evidence and reproduction

The successful source revision is
`9603d7653b216fd09d0e3bdc914e494cf40bb19c`.  The native executable, cubin,
SASS, tables, complete output files, sample ledger, manifests, and logs are in
the durable Modal archive:

- token: `7619dbf0d3534fad9f841e308e4b1f6a`
- archive: 13,971,051 bytes
- SHA-256: `b300054b3f4ba63a59997a70cfa91e4029034ba4682ce89e14df165dbc06637d`

Fetch it with:

```sh
modal run modal_job.py \
  --fetch-token 7619dbf0d3534fad9f841e308e4b1f6a \
  --out /tmp/ecc2k-sparse-table-microprobe-retrieved
```

The independent native audit reopened all 35 ranked rows, the seven resource
rows, seven SASS rows, 64-key histogram, output digests, correctness receipt,
launch identity, and manifest.  It recomputed 11.448092 B/s and the 0.741616
ratio.  Its SHA-256 is
`6366ef4628b2723ce15a929b9116dd388a36128ba04a1fd6081ebf363bbc8145`.
