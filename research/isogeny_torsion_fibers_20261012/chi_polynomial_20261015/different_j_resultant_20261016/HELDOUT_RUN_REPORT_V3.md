# Run record: S3 exhaustive and seven-arm held-out panel

- **Date:** 2026-10-16. **Host:** Darwin arm64, shared L0; no physical-core
  or thermal isolation receipt. **Revision:** feature branch
  `codex/isogeny-portable-trend-20261015`; exact solver, checker, meter,
  script, public and private input SHA-256 values are in the two frozen V2/V3
  specs. **Prior run:** the corrected five-arm V2 and summary-only S3 V1
  remain separate and unchanged.
- **Hypothesis and stop rule:** [protocol](HELDOUT_SEVEN_ARM_PROTOCOL_V3.md).
  Independent exhaustive pair mismatch would stop the seven-arm panel.
- **Exhaustive command:** `target/release/examples/different_j_chi_resultant exhaustive-s3 FOUR_TARGET_PUBLIC_V1.json > S3_EXHAUSTIVE_V2.jsonl`;
  independent `different_j_chi_resultant_check check-s3` command consumed
  that raw JSONL and the public input. Exit statuses were zero. Receipt:
  [S3_EXHAUSTIVE_V2.check.json](S3_EXHAUSTIVE_V2.check.json).
- **Panel command:** `bash RUN_HELDOUT_SEVEN_ARM_V3.sh`; each child used
  `different_j_chi_heldout_meter` and the frozen solver binary, with
  individual JSON, stderr, meter, and exit records in
  [heldout_seven_arm_v3](heldout_seven_arm_v3/). The separate
  `different_j_chi_heldout_check PUBLIC CHECKER DIRECTORY` replayed all
  84 rows, and `different_j_chi_heldout_summary DIRECTORY` reduced them.
- **Limits and failures:** 300-second hard wall timeout per child;
  512 MiB observed peak RSS. No child failed or timed out in this V3
  panel. Earlier Darwin `RLIMIT_AS` and `RLIMIT_DATA` failures are
  preserved under the V2 run. Timing has shared-host noise and no L2
  admission.
- **Results:** [V3 report](HELDOUT_RESULTS_V3.md),
  [native summary](HELDOUT_SEVEN_ARM_V3.summary.json),
  [independent receipt](HELDOUT_SEVEN_ARM_V3.check.json), and
  [graph](HELDOUT_GRAPH_V3.svg). Source and quotient S3 returned the
  exact 229-case pair sets and the four relation arms had identical
  12-cell trial/row support. Portable primitive `S` is unresolved.
