#set page(paper: "a4", margin: 17mm)
#set text(size: 9pt)
= Binary field multiply: two Modal x86-64 CPU panels

*Status:* two completed field-level panels; descriptive timing only.
The full report and raw receipts are in research/ecc_arithmetic_20261013/reports/modal-20261013/REPORT.md.

== Execution and correctness

The same frozen six Rust files and a corrected launcher ran twice on separate
Modal workers. Baseline source SHA-256: 2bc2a1f55887ef7db1fdb1ce5a6cd05c86b7c75e9eb8ade36497ff5686b37fd3.
Candidate source SHA-256: 8b32c5467ee8c87752b9cb86594fb742aa9053fd5918172686c4270c9e936364.
Each worker requested two vCPUs, 2,048 MiB and a 1,200-second timeout. It
reported Linux/gVisor x86-64, PCLMULQDQ and VPCLMULQDQ, but no CPU model or
reserved-core receipt, remote timestamps, or per-phase resource data. Each
passed 34 tests, then completed 24 alternating
AB/BA pairs for four fields at 500,000 public multiplies per arm per pair.
All 384 arm rows are retained; all paired checksums match. The first launcher
attempt failed at a remote source path and has no timed rows.

== Descriptive timing

Median ms for 500,000 public field multiplies, baseline / candidate:

#table(
  columns: (1.2fr, 1fr, 1fr, 1.2fr),
  [Field], [Run 1], [Run 2], [Product count],
  [GF(2^113)], [28.767 / 27.471], [23.996 / 23.200], [4 → 3],
  [GF(2^131)], [36.079 / 33.580], [28.825 / 27.043], [9 → 6],
  [GF(2^163)], [38.869 / 35.389], [32.390 / 29.583], [9 → 6],
  [GF(2^233)], [39.378 / 40.410], [33.142 / 33.702], [unchanged control],
)

#image("../../modal_paired.svg", width: 90%)

Absolute times shifted markedly between workers, including the unchanged
control. The model, placement, frequency and thermal state are unreported;
there is no A/A or production-library comparison. This does not pass the
controlled wall-time gate and establishes no qualified speedup or worldwide
rank. Timing covers public field multiply allocation, dispatch and reduction,
not point/scalar arithmetic or an ECC2K walk.
