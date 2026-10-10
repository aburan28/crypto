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
the controlled speedup field is therefore `null`. The next gate is an
alternated, same-point C2800 mixed/prefetch comparison with per-child resource
counters, followed by an isolated-host one-target panel against C1200 if the
lookup improvement persists. The owner threshold for a large speedup is at
least 2x on the frozen N41/N53 workloads.
