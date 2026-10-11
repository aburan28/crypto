#set page(paper: "a4", margin: (x: 1.8cm, y: 2cm), numbering: "1")
#set text(size: 9.5pt)
#set table(stroke: 0.4pt, inset: 4pt)
#show heading: set block(above: 1.2em, below: 0.6em)

#align(center)[
  #text(size: 15pt, weight: "bold")[Run report: run-20261011T1827Z-sidechannel] \
  Eigen-Séta v2: constant-time verification of the decapsulation secret phase (descriptive record) \
  #text(size: 8.5pt)[Study: research/eigen_seta_rust_20261011 · sources: commit 19c641d121097e80f727bcc650a08632e204ba1d]
]

= What ran

#table(columns: (auto, 1fr),
  [Start (UTC)], [2026-10-11T18:27Z, after the benchmark run finished],
  [Sources], [Commit 19c641d12. The only uncommitted file in `eigen_seta/` was README.md, which is not compiled.],
  [Build], [rustc 1.97.0, release with debuginfo. valgrind-3.22.0 (callgrind), glibc malloc.],
  [Host], [Intel Xeon \@ 2.80 GHz, 4 logical cores, Linux 6.18.44 x86_64, shared cloud VM],
  [Step 1], [`valgrind --tool=callgrind --toggle-collect='*probe_core*' ./ct_probe <key> <ct> <valid|foreign>`. Only `decaps_core` runs inside `probe_core`; the public phase runs outside.],
  [Grid], [Key seeds 1–4, each with valid ct 1, valid ct 2 and a foreign ct. A foreign ct is well formed but made for another key, so FO rejects it.],
  [Summary], [`ctsummary.sh`: callgrind exclusive counts. `alloc_mem` = glibc malloc, arena and mem\* routines. `rest` = total − `alloc_mem`.],
  [Step 2], [`./dudect_decaps 1000`: one key. Fixed ct vs. a pool of 16 cts (8 FO-rejected), interleaved at random. Welch t.],
  [Records], [`records/ctprobe_summary.txt`, `records/callgrind/*.xz`, `records/probe_stdout/`, `records/dudect_output.txt`, `records/host_and_command.txt`, `records/build.log`],
)
Requested: 12 probe runs and 1 timing test. Completed: 13/13. Failures: none.

= Step 1: instructions executed in `probe_core`

#table(columns: 6,
  [*key*], [*ciphertext*], [*accepted*], [*total*], [*alloc_mem*], [*rest*],
  [1], [valid 1], [yes], [1,147,512,199], [31,326,196], [1,116,186,003],
  [1], [valid 2], [yes], [1,147,512,199], [31,326,196], [1,116,186,003],
  [1], [foreign], [no], [1,147,807,708], [31,621,705], [1,116,186,003],
  [2], [valid 1], [yes], [1,147,639,699], [31,453,696], [1,116,186,003],
  [2], [valid 2], [yes], [1,147,639,699], [31,453,696], [1,116,186,003],
  [2], [foreign], [no], [1,147,306,795], [31,120,792], [1,116,186,003],
  [3], [valid 1], [yes], [1,147,928,378], [31,742,375], [1,116,186,003],
  [3], [valid 2], [yes], [1,147,928,378], [31,742,375], [1,116,186,003],
  [3], [foreign], [no], [1,147,958,199], [31,772,196], [1,116,186,003],
  [4], [valid 1], [yes], [1,147,848,114], [31,662,111], [1,116,186,003],
  [4], [valid 2], [yes], [1,147,848,114], [31,662,111], [1,116,186,003],
  [4], [foreign], [no], [1,147,424,282], [31,238,279], [1,116,186,003],
)
Observation: `rest` takes one value in all 12 runs. `alloc_mem` takes 8 distinct values, and the two valid-ct runs of each key agree exactly.

#figure(image("ct_probe.svg", width: 90%), caption: [Instruction counts per run: rest (top) and glibc allocator plus mem\* (bottom). Measured; records/ctprobe_summary.txt.])

= Step 2: timing test (ns per call, 1,000 samples per class)

```
samples/class 1000; median ns: fixed 122309563, random 122318589
t(raw) = +0.75
t(crop <= p50) = +0.84  (n = 501, 500)
t(crop <= p75) = -0.19  (n = 743, 758)
t(crop <= p90) = +0.68  (n = 898, 903)
```

= Scope

Toy16 parameters only, one compiler and one target. The check counts instructions, not memory addresses. The timing test uses one key, so it compares ciphertext classes, not keys.
