# AVX-512 F5 unpack experiment: frozen complete-call comparison

## Hypothesis and reference

The primary `f5_n24_m24_d4` call in the prior direct combination study took
about 56 ms to expand packed pivot rows into monomials, even when echelon
output and the faster row builder and XOR kernel were enabled. The output
representation remains unchanged in this experiment. An AVX-512
compress-store can expand the selected monomials of eight consecutive columns
in one operation. The hypothesis is that the resulting lower unpack cost can
bring the *complete* F5 call to 2× against the current default on the same
host and frozen inputs. The previously measured 1.530× combination is a
motivation, not a speedup factor to multiply by a new phase gain.

The unmodified reference is main commit
`5654c42f6c8d4c8b5f7f25c762097e6fa257ecc2`; its
`src/cryptanalysis/matrix_f5_f2.rs` SHA-256 is
`fbc966fa00f3ef5bab9e49756ef6079db4dfd548f40cb48283843ba8cff3e167`
and `examples/f4_f2_bench.rs` SHA-256 is
`812fdc75bd89c7ac7944a2631158dc02bf1861fed25aa72fd629774daeb6f37a`.
One candidate release binary supplies all three arms in separate processes:

| Arm | `KIC_F5_ECHELON` | `KIC_F5_FUSED_BUILD` | `KIC_GF2_FORCE_AVX2` | `KIC_F5_AVX512_UNPACK` |
| --- | ---: | ---: | ---: | ---: |
| Default reference | 0 | 0 | 0 | 0 |
| Prior combined option | 2 | 1 | 1 | 0 |
| New combined option | 2 | 1 | 1 | 1 |

`KIC_GF2_SIMD=1`, `KIC_GF2_DEFER_ABOVE=0` and
`KIC_GF2_WORD_BATCH=0` in every arm. The prior and new arms use the same
echelon row form; their raw fingerprints and term counts must match. Their
canonical row-space fingerprints must also match the reduced default.

## Frozen workload and measurement

Use all seven F5 cases from `examples/f4_f2_bench.rs` with `1 24 f5`.
`f5_n24_m24_d4` is the sole primary case. The original seed XOR zero and
holdouts `badc0de1`, `5eed2026` and `f5c02a28` are fixed before any timing.
Require a Linux x86-64 host reporting AVX-512F, AVX2 and BMI2; otherwise
record `unsupported_host` and zero timed calls. Pin one allowed CPU and set
`RAYON_NUM_THREADS=1`. For each seed, warm all arms once, take five default
versus default A/A pairs, then five paired triads containing all arms, rotating
the triad order. Do not discard an eligible host or timing sample. Run every
case, retain failures and timeouts, and make no timing-driven changes to this
protocol or candidate.

The charged `wall_ms` is the full matrix-F5 solver call inside the example:
criterion, row build, elimination and unpacking. Fixture construction and
post-call fingerprints are outside that interval. Retain every process's
stdout/stderr, exit status, affinity, load, CPU features, source and binary
digests, internal phase times, rank, pruning, word operations, raw and
canonical fingerprints, and output term count. Compare the five paired
default/new complete-call ratios by their median and exact five-pair
bootstrap 95% interval; also report default/prior and prior/new ratios on
the same host. Preserve the A/A ratio range for each case and seed.

## Gate and stop rule

The 2× one-thread opt-in solver-call target passes only if the primary frozen
median **and** its 95% lower bound exceed 2.00, each holdout primary median
exceeds its own A/A maximum, and no smaller case's default/new median falls
below its own A/A minimum. All seven cases on all four seeds must have equal
canonical row-space fingerprints, rank, pruning and criterion word operations
across arms. Raw fingerprints and output term counts must match between the
prior and new arms; otherwise the experiment fails correctness. If the 2×
gate passes, freeze and run a four-thread control before considering any
default change. If it fails, retain the opt-in code and negative receipt,
remove the one-off timing workflow, and stop this candidate without a
timing-driven rerun. A complete F5 call is a solver-stage diagnostic; this
experiment makes no IC online-time or DLP speedup claim.
