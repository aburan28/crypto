# Run report: run-20261011T1827Z-sidechannel

A descriptive record of the constant-time verification run for the secret
phase of Eigen-Séta v2 decapsulation (`decaps_core` in `eigen_seta/src/kem.rs`).
It states observations only. Interpretation is in `../../README.md` and in
`eigen_seta/SIDE_CHANNELS.md`. The PDF is `report.pdf`.

## What ran

| Field | Value |
|---|---|
| Run id | `run-20261011T1827Z-sidechannel` |
| Start (UTC) | 2026-10-11T18:27Z. The script waited for the benchmark run to finish, then started. |
| Sources | Commit `19c641d121097e80f727bcc650a08632e204ba1d`. The only uncommitted file in `eigen_seta/` was `README.md`, which is not compiled (`records/worktree_status.txt`). |
| Build | rustc 1.97.0, release profile with debuginfo. The build log is `records/build.log`. |
| Tools | valgrind 3.22 (callgrind), glibc malloc |
| Host | Intel Xeon Processor @ 2.80 GHz, 4 logical cores, Linux 6.18.44 x86_64, shared cloud VM |
| Step 1: instruction counts | `valgrind --tool=callgrind --toggle-collect='*probe_core*' ./ct_probe <key_seed> <ct_seed> <valid\|foreign>`, run as 3 parallel streams of 4 probes. Only `decaps_core` is inside `probe_core`. Parsing, point lifting and the canonical bases run outside it. |
| Grid | Key seeds 1–4, each with valid ciphertext seed 1, valid ciphertext seed 2, and a foreign ciphertext. A foreign ciphertext is well formed but made for another key, so the FO check rejects it. |
| Summary | `records/ctsummary.sh` uses `callgrind_annotate` exclusive counts. `alloc_mem` is glibc malloc, arena, `mem*`, `string_fortified` and the unresolved realloc helper. `rest` = total − `alloc_mem`. |
| Step 2: timing test | `./dudect_decaps 1000`. One OsRng key. Class 0 is a fixed ciphertext. Class 1 is a pool of 16 ciphertexts, 8 of them FO-rejected. The classes are interleaved at random, and Welch's t is computed. It ran alone, after step 1. |

Requested: 12 probe runs and 1 timing test. Completed: 13 of 13. Failures:
none.

## Measured values

Step 1 is recorded in `records/ctprobe_summary.txt`, with raw outputs in
`records/callgrind/*.xz` and acceptance in `records/probe_stdout/`. Units are
instructions executed in `probe_core`.

| key seed | ciphertext | accepted | total | glibc allocator and `mem*` | rest |
|---|---|---|---|---|---|
| 1 | valid 1 | yes | 1,147,512,199 | 31,326,196 | 1,116,186,003 |
| 1 | valid 2 | yes | 1,147,512,199 | 31,326,196 | 1,116,186,003 |
| 1 | foreign | no | 1,147,807,708 | 31,621,705 | 1,116,186,003 |
| 2 | valid 1 | yes | 1,147,639,699 | 31,453,696 | 1,116,186,003 |
| 2 | valid 2 | yes | 1,147,639,699 | 31,453,696 | 1,116,186,003 |
| 2 | foreign | no | 1,147,306,795 | 31,120,792 | 1,116,186,003 |
| 3 | valid 1 | yes | 1,147,928,378 | 31,742,375 | 1,116,186,003 |
| 3 | valid 2 | yes | 1,147,928,378 | 31,742,375 | 1,116,186,003 |
| 3 | foreign | no | 1,147,958,199 | 31,772,196 | 1,116,186,003 |
| 4 | valid 1 | yes | 1,147,848,114 | 31,662,111 | 1,116,186,003 |
| 4 | valid 2 | yes | 1,147,848,114 | 31,662,111 | 1,116,186,003 |
| 4 | foreign | no | 1,147,424,282 | 31,238,279 | 1,116,186,003 |

Observation: `rest` takes a single value in all 12 runs. `alloc_mem` takes 8
distinct values. For each key, the two valid-ciphertext runs agree exactly.

Step 2 is recorded in `records/dudect_output.txt`. Units are ns per call,
with 1,000 samples per class.

```
samples/class 1000; median ns: fixed 122309563, random 122318589
t(raw) = +0.75
t(crop <= p50) = +0.84  (n = 501, 500)
t(crop <= p75) = -0.19  (n = 743, 758)
t(crop <= p90) = +0.68  (n = 898, 903)
```

## Graph

`ct_probe.svg` shows `rest` above and `alloc_mem` below for each of the 12
runs. It is drawn from `records/ctprobe_summary.txt` by
`../../figures/plot.py`.

## Scope

This run uses Toy16 parameters only, one compiler, and one target. It counts
instructions, not memory addresses. The timing test uses one key, so it
compares ciphertext classes, not keys.
