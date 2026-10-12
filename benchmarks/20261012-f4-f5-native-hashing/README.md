# Native Boolean F5 hashing: instruction and live-allocation comparison

[Source study](../../research/f4_f5_hash_20261012/REPORT.md) and [combined run report](report.pdf) retain the measured claims. These browsing copies pin commit [31088bfce](https://github.com/aburan28/crypto/tree/31088bfce9d9cc86fb55a655b064b5472db263ae/research/f4_f5_hash_20261012), study tree `e00b12589743ceeb250255f81aa204027e738980` (991 source blobs). [SOURCE.json](SOURCE.json), [SOURCE_TREE.json](SOURCE_TREE.json) and [copies.json](copies.json) record the captured source and exact file hashes. [SOURCE_REPORT.md](SOURCE_REPORT.md) preserves the original report verbatim.

Five paired A/B/A′ rounds across eleven frozen Boolean systems record 3.715% fewer complete-F5 instruction references (95% interval 3.714–3.726%) and 45.136% fewer standalone criterion instruction references (45.121–45.155%). F4 is neutral under the 2% gate. Every result digest and word-XOR count matched; allocation requests, cumulative requested bytes and tracked peak live bytes did not regress.

The units are Callgrind instruction references and allocator-requested bytes. The Linux x86-64 host is an AMD EPYC 9V74 with five logical CPUs; the driver fixes one Rayon worker and Callgrind 3.24 merges worker coverage. Controlled runtime and process RSS are unmeasured; missing CPU-pressure telemetry leaves controlled timing eligibility unearned. These measurements retain their own eleven-system/profiler contract, separate from the older CI history.

The four PDFs and three SVGs below are byte-exact copies from the pinned commit. Each PNG renders its corresponding copied SVG at 150 dpi. Integration added no experiment, native result generation or new measurement.

- [complete-instructions.pdf](complete-instructions.pdf) — `research/f4_f5_hash_20261012/reports/run-20261012T0215Z/complete-instructions.pdf`
- [complete-instructions.svg](complete-instructions.svg) — `research/f4_f5_hash_20261012/reports/run-20261012T0215Z/complete-instructions.svg`
- [criterion-instructions.pdf](criterion-instructions.pdf) — `research/f4_f5_hash_20261012/reports/run-20261012T0215Z/criterion-instructions.pdf`
- [criterion-instructions.svg](criterion-instructions.svg) — `research/f4_f5_hash_20261012/reports/run-20261012T0215Z/criterion-instructions.svg`
- [peak-live.pdf](peak-live.pdf) — `research/f4_f5_hash_20261012/reports/run-20261012T0215Z/peak-live.pdf`
- [peak-live.svg](peak-live.svg) — `research/f4_f5_hash_20261012/reports/run-20261012T0215Z/peak-live.svg`
- [report.pdf](report.pdf) — `research/f4_f5_hash_20261012/reports/run-20261012T0215Z/report.pdf`
