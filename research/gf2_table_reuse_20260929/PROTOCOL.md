# Reuse GF(2) elimination tables: frozen F5 complete-call experiment

## Hypothesis and reference

The direct packed-row F5 experiment reached 1.955× on the frozen primary
complete call, with 63.98 ms still in echelon elimination. Each pivot block
currently clears and zero-fills its entire Gray-code table before rebuilding
it. Every nonzero entry that a row can address is overwritten in that block;
only entry zero of each table must be cleared because table generation reads
it. Reusing the allocation and clearing only those zero entries should
remove table writes without changing a single XOR or row. The remaining
2× gap is small enough that a measurable complete-call saving here may
close it. The hypothesis must be tested directly; memory-traffic estimates
are not a measured gain.

The unmodified reference is commit
`dbc48b990f5cb745f19f0144747a54b3e3e24913` (dependent on merged
PR #919 and pending PR #922). Source SHA-256 values are
`7598ae38a45d07a076995ad2eda10ba684730ea9ea33797c7631ffbd6e6cf388`
for `src/cryptanalysis/gf2_elim.rs`,
`e8c66b4d9163160ed32b848820beb5013eb384ddde09a8bc7b4919565ae93959`
for `src/cryptanalysis/matrix_f5_f2.rs`, and
`35cafca3dc64e8cf25a1ec63b1338eddb8da4725d0ca1c0f9dedc98313cb16aa`
for `examples/f4_f2_bench.rs`.

One release binary runs the three arms in separate processes:

| Arm | F5 form | Fused build | Direct pack | AVX2 XOR | Table reuse |
| --- | --- | ---: | ---: | ---: | ---: |
| Current default | reduced | 0 | 0 | 0 | 0 |
| Prior combined | selective echelon | 1 | 1 | 1 | 0 |
| New combined | selective echelon | 1 | 1 | 1 | 1 |

Set `KIC_F5_AVX512_UNPACK=0`, `KIC_GF2_SIMD=1`,
`KIC_GF2_DEFER_ABOVE=0` and `KIC_GF2_WORD_BATCH=0` in all arms. Only the
new arm sets `KIC_GF2_REUSE_TABLE=1`. The default remains unchanged.

## Frozen workload and measurement

Use all seven `f5` cases from `examples/f4_f2_bench.rs` with `1 24 f5`;
primary `f5_n24_m24_d4`. The original seed XOR zero and holdouts
`badc0de1`, `5eed2026`, `f5c02a28` are fixed. Require a Linux x86-64
host with AVX2 and BMI2, otherwise record `unsupported_host` with zero
timed calls. Pin one allowed CPU and set `RAYON_NUM_THREADS=1`. On each
seed, warm all arms once, run five default/default A/A pairs, then five
rotating three-arm triads. Preserve every process output, status, load,
affinity, CPU feature set, source and binary digest, phase timing, row
signature, rank, pruning, word operation and output term count. No
eligible-host selection or timing-driven repeat is allowed.

`wall_ms` is the full F5 call from criterion through row unpacking.
Fixture construction and post-call fingerprints are outside it. Report
medians of paired default/prior, default/new and prior/new complete-call
ratios with exact five-pair bootstrap 95% intervals and A/A ranges.
Prior/new must match raw and canonical row fingerprints, row and column
counts, rank, pruning, criterion work, reduction word operations and term
counts on every case and seed.

## Gate and stop rule

The 2× goal passes for this one-thread opt-in combination only if frozen
primary default/new median and lower 95% bound both exceed 2.00, every
holdout primary median exceeds its A/A maximum, and no smaller-case
default/new complete-call median falls below its own A/A minimum. If it
passes, freeze and run a four-thread control before considering a default
change. An incremental success short of 2× needs prior/new frozen median
and lower 95% bound above 1.03 with the same holdout and smaller-case
guards; then retain it as an opt-in for further combinations. If that
incremental gate fails, archive the negative result, remove the one-off
workflow and leave the default untouched without a timing-driven rerun.
This is a matrix-F5 solver-stage diagnostic; it makes no IC online-time
or DLP speedup claim.
