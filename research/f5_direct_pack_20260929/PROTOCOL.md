# F5 direct packed-row construction: frozen experiment

## Hypothesis and reference

The frozen AVX-512-unpack study measured 18.9 ms of build time for the
`f5_n24_m24_d4` prior combined arm, after fused row counting. That path
still sorts each polynomial product, materializes all monomial rows, sorts
their flattened terms to discover columns, builds a hash index, then packs
rows. On dense quadratic systems, every monomial of degree at most four
appears. A direct builder can enumerate the full graded column universe
once, compute each monomial's combinatorial column rank, and XOR product
terms straight into the packed row. XOR handles collisions exactly. If a
column is missing, it falls back to the existing builder so the output
contract remains exact on sparse systems. The hypothesis is an incremental
complete-call gain above measurement noise; this is one step toward the
user's 2× goal, not a projected 2× result.

The unmodified reference is commit
`984a089efb015b672e157a57f0ff504e15214e56`, with SHA-256
`056fa31dc32f203e9efaa8dafc7672a632e13efab21340467ef98049ac81d89d`
for `src/cryptanalysis/koblitz_groebner.rs` and
`4615c6c6d04342209002b571cb8ff4ff8ee0574d5547bd1e750a995304ef4cb1`
for `src/cryptanalysis/matrix_f5_f2.rs`. This branch depends on PR #919;
measure after it is merged or against its exact source head. One release
binary supplies all arms in separate processes:

| Arm | F5 form | Fused build | AVX2 XOR | Direct pack | AVX-512 unpack |
| --- | --- | ---: | ---: | ---: | ---: |
| Current default | reduced | 0 | 0 | 0 | 0 |
| Prior combined | selective echelon | 1 | 1 | 0 | 0 |
| New combined | selective echelon | 1 | 1 | 1 | 0 |

Set `KIC_GF2_SIMD=1`, `KIC_GF2_DEFER_ABOVE=0` and
`KIC_GF2_WORD_BATCH=0` in all arms. The direct option is gated to at most
24 variables and degree at most four; it falls back if the full column
universe is not present. It must return the identical packed matrix,
`F5Report` and raw row fingerprint as the prior combined arm.

## Frozen workload and accounting

Use the seven F5 cases emitted by `examples/f4_f2_bench.rs` with
`1 24 f5`, primary `f5_n24_m24_d4`. Use seed XOR zero plus holdouts
`badc0de1`, `5eed2026` and `f5c02a28`. Require a Linux x86-64 host with
AVX2 and BMI2; record `unsupported_host` and zero calls otherwise. Pin one
CPU, set `RAYON_NUM_THREADS=1`, warm all arms once per seed, take five
default/default A/A pairs, then five rotating three-arm triads. Keep all
process outputs, failures, timeouts, host data, source and binary digests,
phase timings, row signatures, ranks, pruning, word operations and term
counts, plus whether the direct path was actually used. The primary new
arm must report direct-path use; a fallback there invalidates the timing
hypothesis. Make no timing-driven change or eligible-host selection.

`wall_ms` is the entire matrix-F5 call: criterion, build, elimination and
unpacking. Input generation and post-call fingerprinting are outside that
interval. For each seed/case, report medians of paired default/prior,
default/new and prior/new complete-call ratios, exact five-pair bootstrap
95% intervals and A/A ranges. Preserve all raw processes in the receipt.

## Correctness, gate and stop rule

Every case on every seed must match canonical row-space fingerprint, rank,
pruning and criterion work across arms. Prior/new must also match raw
fingerprint, output term count and reduction word operations. The
incremental gate requires a frozen primary prior/new median and 95% lower
bound above 1.10, all holdout primary medians above their A/A maxima, and
no smaller complete-call default/new median below its own A/A minimum.
The user's direct 2× target is met only if the frozen primary default/new
median and 95% lower bound both exceed 2.00, with the same holdout and
smaller-case conditions. If that happens, freeze and run a four-thread
control before considering a default change. If the incremental gate
fails, archive the result, remove the one-off timing workflow and stop this
candidate without a timing-driven repeat. If it passes only the
incremental gate, retain the direct-pack opt-in for subsequent combinations
but do not change the default. This measures a solver stage, not IC online
time or DLP speedup.
