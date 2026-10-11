# Run report: run-20261011T1813Z-bench

A descriptive record of one benchmark run of the Eigen-Séta v2 Rust KEM
(`eigen_seta/`). It states observations only. Interpretation is in the study
report, `../../README.md`. The PDF is `report.pdf` in this directory.

## What ran

| Field | Value |
|---|---|
| Run id | `run-20261011T1813Z-bench` |
| Start (UTC) | 2026-10-11T18:13:51Z |
| Sources | `eigen_seta/` at commit `19c641d121097e80f727bcc650a08632e204ba1d`. The binary was built from that tree just before the commit, and a cargo rebuild after the commit was a no-op. |
| Toolchain | rustc 1.97.0 (2d8144b78 2026-07-07), release profile: opt-level 3, thin LTO, 1 codegen unit, default target-cpu |
| Host | Intel Xeon Processor @ 2.80 GHz, 4 logical cores, 15 GiB, Linux 6.18.44 x86_64. Flags include adx, avx2, avx512f, bmi2, popcnt. |
| Isolation | Shared cloud VM. No CPU pinning and no frequency control. One benchmark process at a time, with no other heavy process. |
| Command | `bench_field`, then `run Toy16 5`, `run Toy32 3`, `run L1 3`, `run L3 2`, `run L5 2` (the binaries in `eigen_seta/examples/`) |
| Workload | Per parameter set: one key generation, then the stated number of encapsulation and decapsulation trials. Keys and ciphertexts come from OsRng, so they are fresh for every run. |
| Output pin | Every trial checks that decapsulation returns the encapsulated 32-byte key (`ok true`). |

Runs requested: 5 parameter sets with 15 trials in total, plus the field
micro-benchmark. Runs completed: all of them. Failures: none.

## Measured values

The record is `records/bench_output.txt`, and the host and command are in
`records/host_and_command.txt`. Times are wall clock, one thread.

| Set | keygen (s), n = 1 | encaps (s), per trial | decaps (s), per trial | correct |
|---|---|---|---|---|
| Toy16 | 0.381 | 0.0055, 0.0070, 0.0070, 0.0084, 0.0052 | 0.131, 0.159, 0.186, 0.126, 0.118 | 5/5 |
| Toy32 | 2.341 | 0.0095, 0.0116, 0.0091 | 0.510, 0.446, 0.591 | 3/3 |
| L1 | 32.279 | 0.2448, 0.2588, 0.2394 | 14.223, 14.997, 13.424 | 3/3 |
| L3 | 118.589 | 1.0256, 1.0304 | 51.834, 53.940 | 2/2 |
| L5 | 271.791 | 2.2045, 2.1651 | 123.953, 120.180 | 2/2 |

Public key and ciphertext sizes, from the same record, in bytes:

| Set | Toy16 | Toy32 | L1 | L3 | L5 |
|---|---|---|---|---|---|
| bytes | 216 | 376 | 1344 | 1992 | 2632 |

The secret degree f had 51, 80, 278, 409 and 536 bits in this run's keys,
from Toy16 to L5.

The field micro-benchmark used the same binary. Multiplication is the
minimum over 15 repetitions of 20,000 dependent multiplications. Each
inversion figure is the mean over 200 dependent inversions.

| Set | limbs | mul (ns) | inversion, Fermat (ns) | inversion, safegcd (ns) |
|---|---|---|---|---|
| Toy16 | 4 | 54 | 18,403 | 4,741 |
| Toy32 | 6 | 72 | 67,657 | 10,316 |
| L1 | 21 | 611 | 1,426,022 | 35,608 |
| L3 | 32 | 1,516 | 6,708,390 | 109,506 |
| L5 | 42 | 2,449 | 13,453,052 | 96,327 |

No spread is recorded for the field figures beyond what the table shows.
The L3 safegcd figure is above the L5 one. That is what was recorded, and
this run does not explain it.

## Graphs

* `kem_times.svg` shows wall time per operation and parameter set, with
  every trial plotted on a log scale.
* `field_ops.svg` shows field multiplication and the two inversions
  against the limb count.

Both are drawn from `records/bench_output.txt` by
`../../figures/plot.py`.
