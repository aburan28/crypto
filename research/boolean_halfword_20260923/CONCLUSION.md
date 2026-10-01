# 16-bit syndrome implementation verified; qualified comparison pending

The complete 16-bit syndrome implementation preserves the original model and
assignment order and checks every projected zero on the original equations. The
current source additionally provides runtime-gated three-input-XOR kernels for
both 32-bit and 16-bit lanes, with portable fallbacks. No global CPU feature flag
changes the retained control implementations.

The current code passes 79 Rust correctness checks. Resource/evidence tests cover
the exact historical replay, source identity, full-word agreement, corrupted work,
missing/contended isolation records, false resource flags, and A/A work drift.
These are producer checks, not an independent external review or a performance gain.

## Historical discovery, preserved without promotion

Each of the two September discovery probes completed 24 systems and 10,944
observations across 57 methods. Probe 1's 64-point native variant passed six
incremental groups; probe 2 passed five. Neither produced a dramatic group. Probe
2's outlined mutable-result helper regressed. Its source-based explanation and
the subsequent read-only-helper change are in `SOURCE_EVOLUTION.md`.

The repository's current isolation policy requires pinned/reserved CPUs and A/A
calibration. Those older probes lack that evidence. Their raw timings and original
decisions remain immutable historical diagnostics; they do not discharge the
current performance gate and are not pooled with new runs. Their contention is
unknown. `ACCOUNTING_NOTE.md` corrects one prose exclusion: repeat-signature checks
were actually charged inside validation and total time. No numbers were rewritten.

| Probe 2 arm | Historical groups above 2.0 | Historical groups above 1.0 |
|---|---:|---:|
| half16_native | 0/9 | 0/9 |
| half16_scalar | 0/9 | 0/9 |
| half64_native | 0/9 | 5/9 |
| half64_scalar | 0/9 | 0/9 |
| half_dispatch | 0/9 | 3/9 |

## Preserved n24 timing table

Units are milliseconds per complete generic solve plus validation, pooling two
discovery seeds and eight repetitions. These are **unisolated historical values**,
not eligible new performance evidence. The first planted column retains the initial
figure beside the second probe. New hardware variants have no historical value.

| Method | Probe 1 planted ms | Probe 2 planted ms | Probe 2 cross-planted ms | Probe 2 unplanted ms | Qualification |
|---|---:|---:|---:|---:|---|
| search | 101.416396 | 102.724750 | 108.779604 | 133.462229 | historical; unqualified |
| flat | 408.947938 | 402.651979 | 131.306959 | 595.335729 | historical; unqualified |
| bucket | 1274.376062 | 1380.611875 | 455.870958 | 1732.489979 | historical; unqualified |
| hybrid | 612.422395 | 606.334396 | 202.741584 | 931.873896 | historical; unqualified |
| small_flat | 115.250208 | 118.017167 | 113.367458 | 165.375667 | historical; unqualified |
| word_tail | 104.836500 | 103.836938 | 115.828666 | 151.515250 | historical; unqualified |
| merge_search | 88.800854 | 89.838354 | 85.660083 | 116.898208 | historical; unqualified |
| quadratic_state | 35.610583 | 35.604729 | 39.018521 | 52.941770 | historical; unqualified |
| packed_state | 21.929624 | 22.794167 | 21.761084 | 28.722104 | historical; unqualified |
| basis_list | 284.050230 | 290.962083 | 128.860792 | 479.320791 | historical; unqualified |
| basis_wide | 60.371625 | 62.479000 | 31.137062 | 108.958229 | historical; unqualified |
| tail_list | 200.494249 | 207.070500 | 58.788354 | 244.863625 | historical; unqualified |
| tail_wide | 57.122729 | 60.393000 | 23.265167 | 81.784937 | historical; unqualified |
| packed_untraced | 15.863479 | 15.992812 | 17.645729 | 23.250729 | historical; unqualified |
| affine_sl_basis_list | 216.260250 | 209.653521 | 152.298875 | 530.979958 | historical; unqualified |
| affine_sl_basis_fast | 22.711042 | 23.484375 | 14.810250 | 63.287708 | historical; unqualified |
| gray_scalar | 18.706479 | 19.366521 | 17.761250 | 21.973208 | historical; unqualified |
| gray_simd | 2.714458 | 2.830563 | 2.619209 | 3.183708 | historical; unqualified |
| packed_gray12_scalar | 19.034792 | 19.753979 | 20.877833 | 27.696396 | historical; unqualified |
| packed_gray12_simd | 8.260542 | 8.607146 | 9.273271 | 12.215708 | historical; unqualified |
| packed_gray16_scalar | 15.060521 | 15.662729 | 17.933000 | 22.474791 | historical; unqualified |
| packed_gray16_simd | 2.612291 | 2.668395 | 3.164583 | 3.998813 | historical; unqualified |
| gray_delta_scalar | 18.367687 | 19.442646 | 18.824729 | 22.033021 | historical; unqualified |
| gray_delta_simd | 2.553021 | 2.442208 | 2.441146 | 3.035792 | historical; unqualified |
| initial_list | 18.195354 | 18.004708 | 8.899126 | 21.777896 | historical; unqualified |
| initial_simd | 2.550458 | 2.463146 | 1.330520 | 2.973354 | historical; unqualified |
| packed_gray16_delta_scalar | 15.061917 | 15.157396 | 17.861020 | 22.224354 | historical; unqualified |
| packed_gray16_delta_simd | 2.601084 | 2.672896 | 2.955583 | 3.890708 | historical; unqualified |
| fiber_rows | 4.374042 | 4.612646 | 4.081354 | 9.362208 | historical; unqualified |
| fiber_columns | 3.657437 | 4.676501 | 3.775584 | 8.625813 | historical; unqualified |
| fiber_simd | 1.694438 | 1.935167 | 1.483437 | 3.509500 | historical; unqualified |
| fiber_zero_simd | 1.841355 | 2.001459 | 2.011062 | 4.836042 | historical; unqualified |
| gray_quiet | 2.462063 | 2.432729 | 2.527979 | 2.916270 | historical; unqualified |
| wide64_quiet | 1.603417 | 1.555979 | 1.649167 | 1.903687 | historical; unqualified |
| byte_scalar | 6.590541 | 6.640729 | 6.593458 | 7.650167 | historical; unqualified |
| byte_simd | 1.519313 | 1.584146 | 1.506625 | 2.067687 | historical; unqualified |
| leaf16_quiet | 2.351479 | 2.435563 | 2.613563 | 3.566104 | historical; unqualified |
| leaf16_byte_scalar | 5.719208 | 5.814979 | 6.102021 | 8.599521 | historical; unqualified |
| leaf16_byte_simd | 1.779521 | 1.838062 | 1.792521 | 2.763771 | historical; unqualified |
| byte_single_quiet | 1.984958 | 2.061500 | 1.808833 | 2.379812 | historical; unqualified |
| leaf16_single_quiet | 2.024875 | 2.114187 | 2.036521 | 3.252750 | historical; unqualified |
| byte_planes | 1.745604 | 1.701833 | 1.774854 | 2.129271 | historical; unqualified |
| leaf16_byte_planes | 1.867501 | 1.944958 | 2.274167 | 2.965416 | historical; unqualified |
| byte_unrolled | 1.332876 | 1.268042 | 1.188458 | 1.602438 | historical; unqualified |
| wide64_unrolled | 0.750896 | 0.714376 | 0.719541 | 0.869750 | historical; unqualified |
| leaf16_byte_unrolled | 1.518125 | 1.603145 | 1.474646 | 2.379230 | historical; unqualified |
| word16_unrolled | 1.376980 | 1.316104 | 1.311020 | 1.588083 | historical; unqualified |
| word_dispatch | 0.744750 | 0.733479 | 0.719896 | 0.869792 | historical; unqualified |
| leaf16_word_unrolled | 1.579271 | 1.622770 | 1.782458 | 2.397167 | historical; unqualified |
| projected4 | 2.470396 | 2.574687 | 2.192979 | 4.043209 | historical; unqualified |
| projected5 | 2.781146 | 2.630584 | 2.190895 | 3.333063 | historical; unqualified |
| projected6 | 2.762167 | 2.875146 | 2.200250 | 3.295084 | historical; unqualified |
| half16_scalar | 7.970291 | 9.046854 | 9.009229 | 10.381625 | historical; unqualified |
| half64_scalar | 1.238541 | 1.291416 | 1.294250 | 1.560812 | historical; unqualified |
| half16_native | 1.170334 | 1.880917 | 1.887833 | 2.298729 | historical; unqualified |
| half64_native | 0.597188 | 0.648229 | 0.606313 | 0.739854 | historical; unqualified |
| half_dispatch | 0.596876 | 0.614834 | 0.612375 | 0.737208 | historical; unqualified |
| half16_eor3 | pending | pending | pending | pending | not measured |
| half64_eor3 | pending | pending | pending | pending | not measured |
| eor3_word16 | pending | pending | pending | pending | not measured |
| eor3_word64 | pending | pending | pending | pending | not measured |

## Current measurement contract

`isolated_run.py` refuses unsupported platforms. It freezes the sources and protocol,
keeps the rebuilt unmodified baseline and actual comparison/test binaries, builds
and tests under the repository's busy lock, and then runs A/A before A/B with pinned
CPU reservations. It records host identity/features, memory, compiler, exact commands,
pressure, context switches, faults and other-process CPU use. Any failed, contended
or missing-isolation stage is retained and stops that campaign; it is not pooled
or silently retried. The A/A symmetric noise spread is an additional gate.

The comparison retains the 52 September 23 solver implementations as a frozen
reference roster. It does not claim to include every later repository optimization
or every known solver. Nine new treatments include matched 32-bit EOR3 controls so
the hardware operation is not attributed solely to narrowing. Runtime capability
records distinguish an actual feature path from a portable fallback.

The intended discovery grid has 61 A/B arms and 24 fixtures: 11,712 comparison
observations plus 384 A/A calibration observations. The full grid has 240 fixtures,
117,120 comparison observations and 3,840 calibration observations. It retains all
216 earlier inputs and adds 24 unused holdouts. The full run must bind unchanged
timed source from qualified discovery. A positive primary comparison still requires
confirmation on further unused holdouts. No source tuning follows holdout timing.

Native timing currently requires Linux affinity and pressure interfaces, so the
local macOS correctness checks do not supply new qualified timings. The Linux ARM64
CI route is described in `RESOURCE_PLAN.md`. VM neighbours and frequency remain
outside the tool's control; the measured A/A spread and hardware class must accompany
any result.

`QUALIFIED_RUNS.json` binds any accepted new run; null means no such measurement
is recorded. `RUN_LEDGER.json` binds this report and the preserved probes. All
full-IC, production, calibrated-operation and rho costs remain **null**. The active
dramatic-gain objective is not achieved by this implementation or its tests.
