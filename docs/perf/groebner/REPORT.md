# Separate F4/F5 performance history

Frozen rendering of **16 recorded main revisions**, plus the native allocation baseline/candidate audit. No new performance run or speedup claim.

Plots retain each algorithm as its own series. CI runtime is diagnostic; incomplete historical profiler coverage is retained as recorded. Connecting lines stop when fingerprints or measurement conditions change. X spacing is revision order, not elapsed time.

[Combined PDF](report.pdf) · [Editable source](report.typ) · [Native plot data and provenance](plot-data.json)

## Boolean matrix step: F4 and F5

One shared K_1 / 2^23 m=2 system, 22 Boolean variables, maximum degree 4. Matrix construction, reduction and readback; F5 includes its criterion.

[Group PDF](boolean-matrix.pdf) · [Editable source](boolean-matrix.typ)

![Recorded profiler instructions](boolean-matrix-instructions.svg)

![Historical CI runtime diagnostics](boolean-matrix-runtime.svg)

| UTC date / revision | F4 Ir | F5 Ir | F4 ms | F5 ms |
| --- | --- | --- | --- | --- |
| 2026-10-08 / 93ff0c72 | 860,633,889 | 766,940,469 | 249.241 | 208.250 |
| 2026-10-09 / 6f81260f | 860,555,567 | 766,776,087 | 333.279 | 264.687 |
| 2026-10-11 / 8a41eca7 | 808,263,705 | 718,576,235 | 97.380 | 94.658 |
| 2026-10-11 / 0cfba411 | 808,263,705 | 718,578,523 | 86.390 | 82.214 |
| 2026-10-11 / f88e5909 | 807,989,935 | 718,391,055 | 82.597 | 81.948 |
| 2026-10-11 / fb926c02 | 808,263,705 | 718,576,189 | 101.379 | 97.323 |
| 2026-10-11 / e478b952 | 803,017,657 | 719,516,063 | 107.585 | 102.927 |
| 2026-10-11 / aa0cee7b | 808,281,869 | 718,492,761 | 58.723 | 60.673 |
| 2026-10-12 / 81ca1c78 | 803,015,877 | 719,509,003 | 100.762 | 98.775 |
| 2026-10-12 / a8194d61 | 808,072,539 | 718,282,030 | 101.827 | 96.904 |
| 2026-10-12 / 5cbe16d5 | 802,824,065 | 719,316,911 | 81.784 | 83.144 |
| 2026-10-12 / eb70862e | 802,824,065 | 719,316,876 | 63.944 | 62.505 |
| 2026-10-12 / 31088bfc | 776,242,363 | 686,315,894 | 102.929 | 99.489 |
| 2026-10-12 / 0d34f782 | 774,379,185 | 686,521,131 | 96.521 | 93.790 |
| 2026-10-12 / 89e3233c | 776,490,367 | 686,543,371 | 81.628 | 80.581 |
| 2026-10-12 / 80f5391d | 774,379,185 | 686,519,896 | 98.501 | 93.239 |

- Historical Python-orchestrated perf.history/v1 snapshots; no new performance run or speedup claim.
- Runtime: one recorded median per kernel and revision; sample counts retained per point, no stored spread, no isolated timing eligibility. Runtime is diagnostic.
- Profiler Ir is the recorded historical scope. Worker-thread coverage is incomplete for these kernels; it is not a complete solver operation count.
- Lines stop at missing values, changed output fingerprints, or changed host/compiler/harness/profiler. Different algorithm fingerprints are expected and are not equated.
- X positions are measured revision order, labelled with UTC commit date and revision; spacing is not elapsed time.

## Boolean solver calls: F4 and F5

Shared K_1 / 2^17 m=2 workload, 16 Boolean variables, four planted and four random targets, maximum degree 3. The solve_boolean_system entry point, including root readback.

[Group PDF](boolean-solver.pdf) · [Editable source](boolean-solver.typ)

![Recorded profiler instructions](boolean-solver-instructions.svg)

![Historical CI runtime diagnostics](boolean-solver-runtime.svg)

| UTC date / revision | F4 Ir | F5 Ir | F4 ms | F5 ms |
| --- | --- | --- | --- | --- |
| 2026-10-08 / 93ff0c72 | 1,930,989,200 | 1,953,108,992 | 316.570 | 297.297 |
| 2026-10-09 / 6f81260f | 1,929,374,570 | 1,952,598,745 | 324.403 | 435.004 |
| 2026-10-11 / 8a41eca7 | 1,953,021,023 | 1,974,789,730 | 192.203 | 194.561 |
| 2026-10-11 / 0cfba411 | 1,953,020,766 | 1,975,267,516 | 162.266 | 164.161 |
| 2026-10-11 / f88e5909 | 1,949,981,977 | 1,973,219,339 | 161.982 | 164.584 |
| 2026-10-11 / fb926c02 | 1,953,020,766 | 1,974,789,781 | 191.988 | 195.374 |
| 2026-10-11 / e478b952 | 1,962,376,386 | 1,985,422,672 | 200.058 | 202.996 |
| 2026-10-11 / aa0cee7b | 1,952,402,641 | 1,975,976,609 | 131.376 | 127.933 |
| 2026-10-12 / 81ca1c78 | 1,962,376,386 | 1,985,933,163 | 191.630 | 194.567 |
| 2026-10-12 / a8194d61 | 1,949,694,281 | 1,972,772,104 | 190.899 | 193.455 |
| 2026-10-12 / 5cbe16d5 | 1,961,124,313 | 1,981,687,200 | 163.135 | 164.619 |
| 2026-10-12 / eb70862e | 1,961,124,313 | 1,981,686,883 | 129.169 | 141.847 |
| 2026-10-12 / 31088bfc | 1,783,643,232 | 1,803,568,510 | 184.204 | 186.912 |
| 2026-10-12 / 0d34f782 | 1,787,314,371 | 1,807,797,314 | 176.437 | 178.824 |
| 2026-10-12 / 89e3233c | 1,786,351,448 | 1,806,842,750 | 146.640 | 148.159 |
| 2026-10-12 / 80f5391d | 1,787,314,371 | 1,807,244,283 | 176.353 | 178.878 |

- Historical Python-orchestrated perf.history/v1 snapshots; no new performance run or speedup claim.
- Runtime: one recorded median per kernel and revision; sample counts retained per point, no stored spread, no isolated timing eligibility. Runtime is diagnostic.
- Profiler Ir is the recorded historical scope. Worker-thread coverage is incomplete for these kernels; it is not a complete solver operation count.
- Lines stop at missing values, changed output fingerprints, or changed host/compiler/harness/profiler. Different algorithm fingerprints are expected and are not equated.
- X positions are measured revision order, labelled with UTC commit date and revision; spacing is not elapsed time.

## Complete Boolean calls: allocation history

Matching frozen Boolean allocation suite: 8 seeded systems, 3 repetitions, 24 calls per arm and revision. Polynomial dimensions, degrees and seeds are retained in the source audit protocol.

[Group PDF](complete-allocations.pdf) · [Editable source](complete-allocations.typ)

![Cumulative requested allocation bytes](complete-allocations-bytes.svg)

![Allocation and reallocation requests](complete-allocations-requests.svg)

| UTC date / revision | F4 bytes | F5 bytes | F4 requests | F5 requests |
| --- | --- | --- | --- | --- |
| 2026-10-11 / e957e48c | 11,667,696 | 12,400,920 | 48,414 | 42,801 |
| 2026-10-11 / 2e3dffda | 6,712,608 | 6,762,696 | 48,294 | 39,507 |

- Native allocation audit. Original baseline/candidate rows, accepted comparison, exact outputs and accounting, shared inputs and all totals revalidated by the Rust plot-data tool.
- One baseline and one candidate point per accepted audit, not a daily sequence. A/A and replication pairs are controls, not additional algorithm revisions.
- Bytes are cumulative requested traffic, not peak memory. No timing is measured by this allocation series.
- Criterion construction is a separate part of F5 and is never compared as a complete solver.
- Audit revisions come from their frozen records. They are not inserted as measurements of historical CI workloads; the workloads and metrics remain separate.

## F5 criterion: allocation history

Matching frozen Boolean allocation suite: 8 seeded systems, 3 repetitions, 24 calls per arm and revision. Polynomial dimensions, degrees and seeds are retained in the source audit protocol.

[Group PDF](criterion-allocations.pdf) · [Editable source](criterion-allocations.typ)

![Cumulative requested allocation bytes](criterion-allocations-bytes.svg)

![Allocation and reallocation requests](criterion-allocations-requests.svg)

| UTC date / revision | F5 criterion bytes | F5 criterion requests |
| --- | --- | --- |
| 2026-10-11 / e957e48c | 1,329,024 | 8,748 |
| 2026-10-11 / 2e3dffda | 660,480 | 5,676 |

- Native allocation audit. Original baseline/candidate rows, accepted comparison, exact outputs and accounting, shared inputs and all totals revalidated by the Rust plot-data tool.
- One baseline and one candidate point per accepted audit, not a daily sequence. A/A and replication pairs are controls, not additional algorithm revisions.
- Bytes are cumulative requested traffic, not peak memory. No timing is measured by this allocation series.
- Criterion construction is a separate part of F5 and is never compared as a complete solver.
- Audit revisions come from their frozen records. They are not inserted as measurements of historical CI workloads; the workloads and metrics remain separate.

## Historical provenance

Raw snapshots: `docs/perf/history/` on the `perf-history` branch; frozen copies: `research/f4_f5_history_20261012/records/history/` on the implementation branch. Every SHA-256 below was calculated natively from the original file bytes.

| Snapshot | SHA-256 | Recorded samples |
| --- | --- | --- |
| `20261008T150820Z-93ff0c72333c.json` | `161c9799f49c7c11059a869e0e09ede2d494980a23a32a7b59e4b28e25b3a65e` | 5 |
| `20261009T163719Z-6f81260f6ab5.json` | `325262d8694be901bf434abdf120933f9bcd071e26b63cbdfb9d602c72a4de03` | 5 |
| `20261011T173144Z-8a41eca79320.json` | `48b949e3d6db2930f266b8095297297e870b02abfedb2b4315b411953048062b` | 5 |
| `20261011T181908Z-0cfba41136ff.json` | `e2128df08370e31a34c66a783b16ed10ae7537f22ac396b259a7384441cbf340` | 5 |
| `20261011T204041Z-f88e5909f6cd.json` | `5a5e9ca300950b03163b9dbf71539e05a258e7120d7141c276ed52472ac1e8e4` | 5 |
| `20261011T220542Z-fb926c0206ab.json` | `f3ffddab49ee073b0ffda6ecd2e6513e5292a82035532020c0fb19e027ca9df5` | 5 |
| `20261011T222036Z-e478b952a890.json` | `d3195a2144c1e534f53d41a791c31a48462076a405bf2e4ae6c84a6a5f6967a9` | 5 |
| `20261011T235155Z-aa0cee7bd176.json` | `b9d970e259c97fb5868fd17228cc85bf4c535919e5996d041e5ff0b09baf9859` | 5 |
| `20261012T001856Z-81ca1c78df0a.json` | `3664215805eab84b3bb31f1088c8e404c0bf0d217bafaac661d9ae32e7c1d44a` | 5 |
| `20261012T003523Z-a8194d616ed8.json` | `c0418532025d47b30fffdcac3b9891f5ca8fdd147d4c7f26a047756ebf1eba10` | 5 |
| `20261012T025918Z-5cbe16d515f8.json` | `be4270af1b7d7d5a2c869367e9e438395d3c3ddd2e85e5644d3be14108d93872` | 5 |
| `20261012T040712Z-eb70862e4c69.json` | `fc4d98246f15b49e00dd99e995a019d93ae8de9eec50a95784faad01cc611325` | 5 |
| `20261012T042724Z-31088bfce9d9.json` | `ab40418dc93ef669658769c6b0a9b68615d2b92c344a7ff05c82e62cc15196b2` | 5 |
| `20261012T043519Z-0d34f782ec21.json` | `7807facf6ec02fb9768b0500adab861a905d54aef28ae934ccfb2b12787218ec` | 5 |
| `20261012T044326Z-89e3233c0647.json` | `89b8ab077974d3645755c273a2007ac97dc54a5f8482ea7f9da25cf38a4c4df3` | 5 |
| `20261012T052003Z-80f5391d4f70.json` | `213694afc2098d82e95a9afd842c65dd83155cd9db71d0bc0795435d435dc585` | 5 |

Allocation provenance and hashes are in `plot-data.json`. Millisecond labels are the native conversion `wall_ns / 1,000,000`; no ratios, uncertainty estimates or fitted trends are calculated by the renderer.

Future history regeneration: native `tools/groebner-history`, followed by render-only `scripts/perf/render_groebner_history.py`; the workflow saves SVGs, PDFs, sources and JSON together.
