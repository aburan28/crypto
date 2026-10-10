# C2800 batched S3 target lookup prefetch

The exact N53 C2800 target extractor now computes the canonical partner keys
for an eight-state batch, prefetches their root-table home slots, and checks
the candidate witnesses in the original state/root order. The frozen T025
control recovered the expected scalar with the same first relation witness,
complete rank 2800, and 24.08 GiB peak RSS in a 32 GiB envelope.

## Implementation and identity

- Source: `examples/s3_partial_prepared_fullrank_c2800prefetch_targetspan.rs`,
  commit `0e91397c7`, SHA-256
  `43705e308be40c6340c67c96312bb77ebe4fa4bc07ec16ef2a769bdf8a276422`.
- Candidate:
  `IC1N53Ckb1fb296800PDP4rootRCsampleLAgaussTDdirectISO0h9d7b111b5b0f`.
  The exact manifest is under `identities/candidates/`.
- Workload: N53 public point T025, `Wa238237883e9`, from the cross-curve
  freeze. One target is supplied as a public point; its known fixture scalar
  is checked after extraction.
- The implementation handles at most 32 partner keys per eight-state batch.
  It issues an L1 read prefetch for each table home slot on ARM64/x86-64 and
  uses a portable no-op on other architectures. Group lifting still follows
  the prior lookup order. Look-ahead S3 calls in the final batch are counted.

## Frozen T025 gate

`CONTROL_FREEZE.json` and `run_control.py` produced the first row. The solver
completed with matching stdout/report and exact witness, but macOS
`/usr/bin/time -l` exited 1 because its `sysctl kern.clockrate` probe was
denied in the sandbox. That row is preserved under
`local/first_wrapper_failure/`; its RSS is unavailable and its run status is
failed at the wrapper layer.

`CONTROL_FREEZE_V2.json` and `run_control_v2.py` launched the solver directly
and sampled child `getrusage`. The verified row and compressed raw report are
under `local/verified_t025/`. The new report differs from the prior mixed
C2800 report only in timings, RSS, and the three counters affected by
final-batch look-ahead (`target_state_probes`, `target_s3_calls`, and
`total_s3_calls`). Both reports use relation indices
`[55035, 67241, 57152, 189185]` and 995 partner probes. The prefetch row
records 256 state probes and 512 S3 calls; the prior row records 249 and 498.

| T025 execution | Solver exit | Online ms | Charged target PDP ms | Peak RSS GiB | Status |
| --- | ---: | ---: | ---: | ---: | --- |
| Prior mixed C2800 | 0 | 5.110 | 4.412 | 24.0 | Verified reference |
| Prefetch, time wrapper | 1 | 1.818 | 1.633 | unavailable | Wrapper failure; solver output matched |
| Prefetch, direct child | 0 | 0.310 | 0.161 | 24.08 | Verified control |

These local target intervals vary sharply even with identical source and
input. They are exploratory diagnostics. Host-wide exclusive CPU allocation,
frequency, IRQ, throttling, steal-time, and paired order were not established;
the controlled speedup field is therefore `null`.

## Complete eight-point paired screen

`SCREEN_FREEZE.json` fixed T173–T180, one previously unseen public point per
run, one repetition per arm, the exact source and binary hashes, 32 GiB memory
envelope, and alternating arm order. All 16 solver runs verified the fixture
scalar, reached rank 2800 before target input, stayed within the memory
envelope, and gave the same first relation witness within each pair. The full
raw ledger and 16 reports are under `local/paired_screen/`; the archive
receipt binds each compressed artifact by SHA-256.

| T | Mixed online ms | Prefetch online ms | Mixed/prefetch online | Mixed/prefetch charged PDP | Setup negative-control ratio |
| --- | ---: | ---: | ---: | ---: | ---: |
| 173 | 0.482 | 0.270 | 1.78x | 3.01x | 1.55x |
| 174 | 4.782 | 0.311 | 15.38x | 12.38x | 0.61x |
| 175 | 0.775 | 0.501 | 1.55x | 6.55x | 1.44x |
| 176 | 1.302 | 0.940 | 1.39x | 0.57x | 0.95x |
| 177 | 0.308 | 0.953 | 0.32x | 1.76x | 0.97x |
| 178 | 0.537 | 0.202 | 2.65x | 1.76x | 1.00x |
| 179 | 0.371 | 0.566 | 0.66x | 0.60x | 1.11x |
| 180 | 0.929 | 2.888 | 0.32x | 0.15x | 0.96x |

The target-resampled online geometric mean is **1.34x** with bootstrap 95%
interval **[0.62x, 3.17x]**. Charged PDP geometric mean is **1.58x** with
interval **[0.63x, 3.86x]**. Peak RSS was 25,859,375,104 bytes (24.08 GiB).
The screen's prespecified exploratory advancement threshold was 2x, so this
lookup placement does not advance. The unchanged reusable setup differed by
more than 10% in several pairs; its negative control fails, and no controlled
CPU speedup is inferred from these local timings. The full five-phase costs,
per-S3-call diagnostics, and individual failures/statuses are in
`analysis.json`.

The next candidate issues each root-slot prefetch immediately when its
canonical key is formed, overlapping memory fetch with the remaining S3
root arithmetic in the batch. Its source is
`examples/s3_partial_prepared_fullrank_c2800earlyprefetch_targetspan.rs`.
That separate variant must pass T025 scalar, witness, rank, and memory gates
before a paired screen. The owner threshold for a large speedup remains at
least 2x on the frozen N41/N53 workloads.
