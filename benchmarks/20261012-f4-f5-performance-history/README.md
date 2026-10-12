# Separate F4/F5 historical instruction/runtime/allocation plots

[Source study](../../research/f4_f5_history_20261012) and [run/result record](../../research/f4_f5_history_20261012/reports/run-20261012T0009Z/report.pdf) remain the source of truth. The PDF and SVG files below are byte-identical browsing copies added during integration on 2026-10-12. PNGs render the corresponding exact SVG at 150 dpi; no experiment or native result generation was run. The [copy manifest](copies.json) records source paths and hashes.

Seven frozen perf-history snapshots; F4/F5 instruction and diagnostic-runtime series, plus distinct allocation groups; no new timing measurement

Host and isolation: Historical hosts/compilers/profilers differ; connections break at changed fingerprints; worker-thread coverage incomplete. Original host/isolation wording and measurement qualifications are preserved in the copied PDFs and [verbatim source report](SOURCE_REPORT.md).

History rendering quotes existing records; runtime has no stored spread or isolation eligibility.

## Visual copies

- [boolean-matrix-instructions.svg](visuals/research/f4_f5_history_20261012/reports/run-20261012T0009Z/boolean-matrix-instructions.svg) — `research/f4_f5_history_20261012/reports/run-20261012T0009Z/boolean-matrix-instructions.svg`
- [boolean-matrix-runtime.svg](visuals/research/f4_f5_history_20261012/reports/run-20261012T0009Z/boolean-matrix-runtime.svg) — `research/f4_f5_history_20261012/reports/run-20261012T0009Z/boolean-matrix-runtime.svg`
- [boolean-matrix.pdf](visuals/research/f4_f5_history_20261012/reports/run-20261012T0009Z/boolean-matrix.pdf) — `research/f4_f5_history_20261012/reports/run-20261012T0009Z/boolean-matrix.pdf`
- [boolean-solver-instructions.svg](visuals/research/f4_f5_history_20261012/reports/run-20261012T0009Z/boolean-solver-instructions.svg) — `research/f4_f5_history_20261012/reports/run-20261012T0009Z/boolean-solver-instructions.svg`
- [boolean-solver-runtime.svg](visuals/research/f4_f5_history_20261012/reports/run-20261012T0009Z/boolean-solver-runtime.svg) — `research/f4_f5_history_20261012/reports/run-20261012T0009Z/boolean-solver-runtime.svg`
- [boolean-solver.pdf](visuals/research/f4_f5_history_20261012/reports/run-20261012T0009Z/boolean-solver.pdf) — `research/f4_f5_history_20261012/reports/run-20261012T0009Z/boolean-solver.pdf`
- [complete-allocations-bytes.svg](visuals/research/f4_f5_history_20261012/reports/run-20261012T0009Z/complete-allocations-bytes.svg) — `research/f4_f5_history_20261012/reports/run-20261012T0009Z/complete-allocations-bytes.svg`
- [complete-allocations-requests.svg](visuals/research/f4_f5_history_20261012/reports/run-20261012T0009Z/complete-allocations-requests.svg) — `research/f4_f5_history_20261012/reports/run-20261012T0009Z/complete-allocations-requests.svg`
- [complete-allocations.pdf](visuals/research/f4_f5_history_20261012/reports/run-20261012T0009Z/complete-allocations.pdf) — `research/f4_f5_history_20261012/reports/run-20261012T0009Z/complete-allocations.pdf`
- [criterion-allocations-bytes.svg](visuals/research/f4_f5_history_20261012/reports/run-20261012T0009Z/criterion-allocations-bytes.svg) — `research/f4_f5_history_20261012/reports/run-20261012T0009Z/criterion-allocations-bytes.svg`
- [criterion-allocations-requests.svg](visuals/research/f4_f5_history_20261012/reports/run-20261012T0009Z/criterion-allocations-requests.svg) — `research/f4_f5_history_20261012/reports/run-20261012T0009Z/criterion-allocations-requests.svg`
- [criterion-allocations.pdf](visuals/research/f4_f5_history_20261012/reports/run-20261012T0009Z/criterion-allocations.pdf) — `research/f4_f5_history_20261012/reports/run-20261012T0009Z/criterion-allocations.pdf`
- [report.pdf](visuals/research/f4_f5_history_20261012/reports/run-20261012T0009Z/report.pdf) — `research/f4_f5_history_20261012/reports/run-20261012T0009Z/report.pdf`

## Captured planning source

The roadmap references in the verbatim source report describe the original planning snapshot captured at PR #1677 head [baa624218809b821400cfe33155c6bacb542e9a8](https://github.com/aburan28/crypto/tree/baa624218809b821400cfe33155c6bacb542e9a8/research/f4_f5_history_20261012). Its five original planning files are preserved in the [captured planning archive](../../research/f4_f5_history_20261012/archives/pr-1677-baa6242188/), with source hashes and provenance. The original planning PDF has SHA-256 `c09a1753fdf9d9af9caebcef48c87e3ba15434007a7f3eb94966b4fdf1a2c3e1`.

The [current roadmap](../../research/f4_f5_history_20261012/NEXT_STEPS.md) records the later native hashing execution separately; its refreshed planning PDF has SHA-256 `2ce7eb49b71ef823a25481f47af4503d7903407e051b1eae3a952c6e33093eb3`. The historical measured plots and browsing PDF/SVG copies in this entry retain their original bytes and measurement scope. The [native hashing entry](../20261012-f4-f5-native-hashing/) carries that follow-up's instruction and live-allocation reports.
