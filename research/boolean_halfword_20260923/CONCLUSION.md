# 16-bit specialization discovery complete; full comparison pending

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
Before startup and every fixture, a bounded readiness wait retains every rejected
and accepted observation at unchanged CPU/PSI thresholds. The locked worker launch
checks again. Waiting contributes to campaign duration; no timed solve is retried.

The comparison retains the 52 September 23 solver implementations as a frozen
reference roster. It does not claim to include every later repository optimization
or every known solver. Nine new treatments include matched 32-bit EOR3 controls so
the hardware operation is not attributed solely to narrowing. Runtime capability
records distinguish an actual feature path from a portable fallback.

The current paired-worker discovery grid has 61 A/B arms and 24 fixtures: 14,640
comparison observations plus 480 A/A calibration observations. It uses sixteen
repetitions at n12 and eight at larger sizes. The full grid has 240 fixtures,
146,400 comparison observations and 4,800 calibration observations. It retains all
216 earlier inputs and adds 24 unused holdouts. The full run must bind unchanged
timed source from qualified discovery. A positive primary comparison still requires
confirmation on further unused holdouts. No source tuning follows holdout timing.

Schema 3 runs A/A followed by A/B for one fixture inside one reserved worker.
The resource receipt covers that real paired computation, while each cold solve
retains its own timer. The combined stdout and exact phase slices are preserved.
No padding or resource-threshold relaxation is used. This packaging change follows
the retained short-worker resource failures and needs a fresh qualified run.

Native timing currently requires Linux affinity and pressure interfaces, so the
local macOS correctness checks do not supply new qualified timings. The Linux ARM64
CI route is described in `RESOURCE_PLAN.md`. VM neighbours and frequency remain
outside the tool's control; the measured A/A spread and hardware class must accompany
any result.

`QUALIFIED_RUNS.json` binds any accepted new run; null means no such measurement
is recorded. `RUN_LEDGER.json` binds this report and the preserved probes. All
full-IC, production, calibrated-operation and rho costs remain **null**. The active
dramatic-gain objective is not achieved by this implementation or its tests.

## Qualified predecessor: fa01d9b8d63ef2f489a6497ffaea33342a3dd27c

The preserved Linux ARM64 discovery completed 24 fixtures, 11,712 comparison observations and 384 A/A observations. Every resource record passed. The 64-point native policy is about 1.55x the retained dispatcher on pooled n24 medians, but no dramatic group passed. These results belong to the predecessor source, before the complete-budget specialization; they do not measure the current source or discharge the full/holdout gate.

All values below are n24 discovery milliseconds per cold solve plus validation. The ratio is the pooled planted dispatcher median divided by the arm median; the paired gates and A/A noise floors are retained in the result.

| Method | Planted ms | Cross-planted ms | Unplanted ms | Dispatcher / arm, planted | Correctness |
|---|---:|---:|---:|---:|---|
| search | 105.719880 | 107.963231 | 137.869802 | 0.0199 | PASS |
| flat | 476.711574 | 163.982665 | 626.913101 | 0.0044 | PASS |
| bucket | 1340.397402 | 470.012779 | 1788.952407 | 0.0016 | PASS |
| hybrid | 703.921217 | 236.472423 | 999.120916 | 0.0030 | PASS |
| small_flat | 129.712210 | 131.661622 | 168.501168 | 0.0162 | PASS |
| word_tail | 112.532685 | 114.386076 | 145.444884 | 0.0187 | PASS |
| merge_search | 84.860065 | 86.962409 | 109.129018 | 0.0248 | PASS |
| quadratic_state | 46.918948 | 47.823172 | 60.066727 | 0.0449 | PASS |
| packed_state | 26.498610 | 26.852876 | 34.119647 | 0.0795 | PASS |
| basis_list | 330.944323 | 155.390369 | 519.508548 | 0.0064 | PASS |
| basis_wide | 91.233736 | 43.040775 | 141.659374 | 0.0231 | PASS |
| tail_list | 234.058086 | 72.263843 | 297.199323 | 0.0090 | PASS |
| tail_wide | 80.957523 | 26.014319 | 100.748394 | 0.0260 | PASS |
| packed_untraced | 19.078602 | 19.308741 | 24.453173 | 0.1105 | PASS |
| affine_sl_basis_list | 223.632966 | 163.458516 | 557.585184 | 0.0094 | PASS |
| affine_sl_basis_fast | 31.086813 | 19.759242 | 77.204048 | 0.0678 | PASS |
| gray_scalar | 32.100125 | 32.093979 | 39.265472 | 0.0657 | PASS |
| gray_simd | 6.748621 | 6.719289 | 8.213806 | 0.3123 | PASS |
| packed_gray12_scalar | 32.352067 | 34.004120 | 43.774976 | 0.0651 | PASS |
| packed_gray12_simd | 13.714227 | 14.593191 | 19.072804 | 0.1537 | PASS |
| packed_gray16_scalar | 27.407441 | 30.200672 | 40.184866 | 0.0769 | PASS |
| packed_gray16_simd | 6.273331 | 6.906911 | 9.209210 | 0.3359 | PASS |
| gray_delta_scalar | 24.521331 | 24.515345 | 30.102401 | 0.0859 | PASS |
| gray_delta_simd | 3.913597 | 3.908260 | 4.775121 | 0.5385 | PASS |
| initial_list | 24.561029 | 12.303752 | 29.927939 | 0.0858 | PASS |
| initial_simd | 3.937794 | 1.977740 | 4.798617 | 0.5352 | PASS |
| packed_gray16_delta_scalar | 21.149373 | 23.319103 | 31.187461 | 0.0996 | PASS |
| packed_gray16_delta_simd | 4.100207 | 4.520645 | 6.026293 | 0.5140 | PASS |
| fiber_rows | 6.743840 | 6.112463 | 13.364536 | 0.3125 | PASS |
| fiber_columns | 5.994147 | 5.611097 | 12.768343 | 0.3516 | PASS |
| fiber_simd | 1.769000 | 1.599959 | 3.772489 | 1.1913 | PASS |
| fiber_zero_simd | 2.413365 | 2.393546 | 5.887454 | 0.8732 | PASS |
| gray_quiet | 6.022172 | 6.041904 | 7.335682 | 0.3499 | PASS |
| wide64_quiet | 2.936964 | 2.966862 | 3.619991 | 0.7175 | PASS |
| byte_scalar | 7.181916 | 7.130006 | 8.784897 | 0.2934 | PASS |
| byte_simd | 2.092758 | 2.005104 | 2.614246 | 1.0070 | PASS |
| leaf16_quiet | 5.745887 | 6.225156 | 8.307172 | 0.3668 | PASS |
| leaf16_byte_scalar | 6.456023 | 6.950722 | 9.460914 | 0.3264 | PASS |
| leaf16_byte_simd | 2.476602 | 2.449918 | 3.698079 | 0.8509 | PASS |
| byte_single_quiet | 2.536446 | 2.424253 | 3.176274 | 0.8308 | PASS |
| leaf16_single_quiet | 2.855276 | 2.815314 | 4.287693 | 0.7381 | PASS |
| byte_planes | 1.997032 | 1.930870 | 2.472877 | 1.0553 | PASS |
| leaf16_byte_planes | 2.459707 | 2.507347 | 3.669624 | 0.8568 | PASS |
| byte_unrolled | 1.841151 | 1.720751 | 2.328106 | 1.1446 | PASS |
| wide64_unrolled | 2.109215 | 2.106968 | 2.576904 | 0.9991 | PASS |
| leaf16_byte_unrolled | 2.275157 | 2.214836 | 3.366536 | 0.9263 | PASS |
| word16_unrolled | 3.269214 | 3.263733 | 3.995912 | 0.6446 | PASS |
| word_dispatch | 2.107387 | 2.108615 | 2.573934 | 1.0000 | PASS |
| leaf16_word_unrolled | 3.455238 | 3.810814 | 5.065165 | 0.6099 | PASS |
| projected4 | 3.041375 | 2.803546 | 4.244109 | 0.6929 | PASS |
| projected5 | 4.061677 | 3.271200 | 5.014383 | 0.5188 | PASS |
| projected6 | 5.470330 | 4.220270 | 6.572756 | 0.3852 | PASS |
| half16_scalar | 11.201064 | 11.197945 | 13.700282 | 0.1881 | PASS |
| half64_scalar | 2.796897 | 2.181459 | 2.658514 | 0.7535 | PASS |
| half16_native | 3.520116 | 3.519623 | 4.313834 | 0.5987 | PASS |
| half64_native | 1.359292 | 1.363018 | 1.665304 | 1.5504 | PASS |
| half_dispatch | 1.359587 | 1.361387 | 1.659069 | 1.5500 | PASS |
| half16_eor3 | 3.715054 | 3.726886 | 4.549054 | 0.5673 | PASS |
| half64_eor3 | 1.385261 | 1.397932 | 1.694177 | 1.5213 | PASS |
| eor3_word16 | 3.146378 | 3.144396 | 3.832539 | 0.6698 | PASS |
| eor3_word64 | 1.755348 | 1.750443 | 2.137465 | 1.2006 | PASS |

## Current qualified discovery

Source `127f9830252f238f3fd5fcfb99a046bb9daf5ff2` completed 24 fixtures, 14,640 A/B observations and 480 A/A observations. All resource receipts passed and every arm completed verification. Discovery only; no promotion. The group A/A symmetric noise floors range from 1.0037 to 1.4500. Campaign duration, including readiness waits, was 351.248 seconds; peak worker RSS was 14,295,040 bytes. Host: Linux aarch64, 4 logical CPUs; reserved CPUs [3]; EOR3 available: True. Host/compiler details and exact per-group confidence intervals remain in the sealed artifact.

The gate uses the pointwise fastest of all 52 frozen reference methods and requires the lower 95% bound to exceed the A/A floor as well as the numeric threshold.

Predecessor and current campaigns used separate hosted VMs. Their absolute before/after times do not isolate the effect of budget specialization; only within-campaign comparisons to the frozen reference roster are admitted here.

| Candidate | Groups above 2x | Groups above 1x | Primary dramatic verdict |
|---|---:|---:|---|
| eor3_word16 | 0/9 | 0/9 | discovery only |
| eor3_word64 | 0/9 | 2/9 | discovery only |
| half16_eor3 | 0/9 | 0/9 | discovery only |
| half16_native | 0/9 | 0/9 | discovery only |
| half16_scalar | 0/9 | 0/9 | discovery only |
| half64_eor3 | 0/9 | 4/9 | discovery only |
| half64_native | 0/9 | 4/9 | discovery only |
| half64_scalar | 0/9 | 0/9 | discovery only |
| half_dispatch | 0/9 | 1/9 | discovery only |

All methods below use n24 discovery milliseconds per complete solve plus validation. The displayed ratio divides pooled planted dispatcher medians; it is descriptive, not the primary fastest-reference gate.

| Method | Planted ms | Cross-planted ms | Unplanted ms | Dispatcher / arm, planted | Correctness |
|---|---:|---:|---:|---:|---|
| search | 104.553748 | 106.881303 | 136.359259 | 0.0200 | PASS |
| flat | 487.309635 | 167.425981 | 641.582068 | 0.0043 | PASS |
| bucket | 1335.511633 | 468.448710 | 1783.931259 | 0.0016 | PASS |
| hybrid | 711.492731 | 239.108532 | 1007.641847 | 0.0029 | PASS |
| small_flat | 129.447555 | 131.286347 | 167.912848 | 0.0162 | PASS |
| word_tail | 111.626884 | 113.428723 | 144.069997 | 0.0187 | PASS |
| merge_search | 85.081701 | 87.188851 | 109.238916 | 0.0246 | PASS |
| quadratic_state | 47.169623 | 48.128025 | 60.369006 | 0.0444 | PASS |
| packed_state | 26.777940 | 27.083433 | 34.339220 | 0.0781 | PASS |
| basis_list | 330.361128 | 155.108903 | 518.784663 | 0.0063 | PASS |
| basis_wide | 91.609356 | 43.146001 | 142.323885 | 0.0228 | PASS |
| tail_list | 233.379305 | 72.119328 | 296.231722 | 0.0090 | PASS |
| tail_wide | 80.487677 | 25.895419 | 100.309107 | 0.0260 | PASS |
| packed_untraced | 18.984797 | 19.275845 | 24.459836 | 0.1102 | PASS |
| affine_sl_basis_list | 224.297336 | 163.738504 | 558.450240 | 0.0093 | PASS |
| affine_sl_basis_fast | 31.090414 | 19.872093 | 76.315011 | 0.0673 | PASS |
| gray_scalar | 32.051380 | 32.052146 | 39.271082 | 0.0653 | PASS |
| gray_simd | 6.742654 | 6.760577 | 8.253766 | 0.3103 | PASS |
| packed_gray12_scalar | 32.346322 | 33.999190 | 43.771149 | 0.0647 | PASS |
| packed_gray12_simd | 13.729873 | 14.623010 | 19.179692 | 0.1524 | PASS |
| packed_gray16_scalar | 27.409869 | 30.170752 | 40.191046 | 0.0763 | PASS |
| packed_gray16_simd | 6.294831 | 6.922709 | 9.248537 | 0.3323 | PASS |
| gray_delta_scalar | 24.439484 | 24.384396 | 29.814465 | 0.0856 | PASS |
| gray_delta_simd | 3.860983 | 3.848304 | 4.732903 | 0.5419 | PASS |
| initial_list | 24.411808 | 12.183588 | 29.825623 | 0.0857 | PASS |
| initial_simd | 3.851884 | 1.942227 | 4.758058 | 0.5431 | PASS |
| packed_gray16_delta_scalar | 21.007015 | 23.223705 | 30.798582 | 0.0996 | PASS |
| packed_gray16_delta_simd | 4.062911 | 4.451808 | 5.952367 | 0.5149 | PASS |
| fiber_rows | 6.649902 | 6.008488 | 13.163386 | 0.3146 | PASS |
| fiber_columns | 5.878600 | 5.517244 | 12.580246 | 0.3559 | PASS |
| fiber_simd | 1.762293 | 1.595330 | 3.837722 | 1.1871 | PASS |
| fiber_zero_simd | 2.432231 | 2.441140 | 5.981624 | 0.8601 | PASS |
| gray_quiet | 5.956030 | 5.930328 | 7.210914 | 0.3513 | PASS |
| wide64_quiet | 2.960438 | 2.900965 | 3.635572 | 0.7067 | PASS |
| byte_scalar | 6.914671 | 6.823484 | 8.361628 | 0.3026 | PASS |
| byte_simd | 2.086683 | 2.001860 | 2.600838 | 1.0026 | PASS |
| leaf16_quiet | 5.595007 | 6.149719 | 8.191797 | 0.3739 | PASS |
| leaf16_byte_scalar | 6.603875 | 6.918960 | 9.700167 | 0.3168 | PASS |
| leaf16_byte_simd | 2.470855 | 2.455643 | 3.697033 | 0.8467 | PASS |
| byte_single_quiet | 2.556373 | 2.442363 | 3.218628 | 0.8184 | PASS |
| leaf16_single_quiet | 2.886058 | 2.833537 | 4.321749 | 0.7249 | PASS |
| byte_planes | 2.014807 | 1.953040 | 2.493056 | 1.0384 | PASS |
| leaf16_byte_planes | 2.480596 | 2.531086 | 3.689154 | 0.8434 | PASS |
| byte_unrolled | 1.850756 | 1.721880 | 2.333818 | 1.1304 | PASS |
| wide64_unrolled | 2.105225 | 2.113662 | 2.569311 | 0.9938 | PASS |
| leaf16_byte_unrolled | 2.279297 | 2.222987 | 3.375850 | 0.9179 | PASS |
| word16_unrolled | 3.280844 | 3.281335 | 3.994311 | 0.6377 | PASS |
| word_dispatch | 2.092079 | 2.103739 | 2.573336 | 1.0000 | PASS |
| leaf16_word_unrolled | 3.477716 | 3.802679 | 5.093454 | 0.6016 | PASS |
| projected4 | 3.068396 | 2.826070 | 4.277967 | 0.6818 | PASS |
| projected5 | 4.093949 | 3.305119 | 5.053212 | 0.5110 | PASS |
| projected6 | 5.345794 | 4.072516 | 6.501048 | 0.3914 | PASS |
| half16_scalar | 11.216853 | 11.212701 | 13.743119 | 0.1865 | PASS |
| half64_scalar | 2.202784 | 2.207496 | 2.685017 | 0.9497 | PASS |
| half16_native | 3.552208 | 3.535434 | 4.358923 | 0.5890 | PASS |
| half64_native | 1.362587 | 1.402410 | 1.625218 | 1.5354 | PASS |
| half_dispatch | 1.359452 | 1.401814 | 1.622493 | 1.5389 | PASS |
| half16_eor3 | 3.716743 | 3.723616 | 4.549493 | 0.5629 | PASS |
| half64_eor3 | 1.373667 | 1.388459 | 1.672346 | 1.5230 | PASS |
| eor3_word16 | 3.071503 | 3.066976 | 3.755245 | 0.6811 | PASS |
| eor3_word64 | 1.767829 | 1.767456 | 2.162423 | 1.1834 | PASS |

## Retained execution and validation failures

- GitHub run 36922805760, attempt 1: Other-process CPU exceeded the unchanged isolation threshold in n12/seed17/unplanted A/B; no performance result admitted. The complete artifact is retained in `failed_isolation_01` and contributes no accepted timing samples.
- GitHub run 36930298574, attempt 1: Kernel RCU CPU tick exceeded the unchanged isolation threshold during n12/seed17/cross-planted A/A; no performance result admitted. The complete artifact is retained in `failed_isolation_02` and contributes no accepted timing samples.
- GitHub run 36930298574, attempt 2: Same-source retry stopped at n20/seed17/unplanted A/A because the unchanged resource threshold was exceeded; no samples admitted. The complete artifact is retained in `failed_isolation_03` and contributes no accepted timing samples.
- GitHub run 36935271707, attempt 1: After 19 completed fixture pairs, n24/seed17/cross-planted was refused before worker launch because CPU PSI avg10 was 18.36, above the unchanged limit of 5.0. No performance samples admitted. The complete artifact is retained in `failed_isolation_04` and contributes no accepted timing samples.
- GitHub run 36936730092, attempt 1: All 24 fixture pairs passed resource admission; the final verifier used the separate-worker suffix for a paired A/A receipt and raised FileNotFoundError. Raw evidence is preserved; repaired analysis remains diagnostic and is not an admitted discovery binding. The complete artifact is retained in `failed_analysis_01` and contributes no accepted timing samples.
- GitHub run 36938560536, attempt 1: After eight fixture pairs, the n16/seed17/unplanted paired worker completed normally but provisioning-service CPU use raised total other-process CPU to 0.17 seconds in 1.129265813 seconds, above the unchanged 10% threshold. No performance samples admitted. The complete artifact is retained in `failed_isolation_05` and contributes no accepted timing samples.
- GitHub run 36943520046, attempt 1: Full run stopped after 44 of 240 fixture pairs before holdouts: n12/regression/seed20261022/unplanted worker exited normally in 0.075414677 seconds, but 0.01 seconds of other-process CPU exceeded the unchanged 10% limit. No full performance sample is admitted or pooled. The complete artifact is retained in `failed_full_01` and contributes no accepted timing samples.

`ISOLATION_ATTEMPTS.md` records the exact failure and any subsequent complete same-source retry.
