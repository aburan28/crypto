# Run record: nine-arm held-out fiber comparison

- **Date:** 2026-10-16. **Host:** Darwin arm64, shared L0; no reserved
  physical core, affinity or thermal receipt. **Branch:**
  `codex/isogeny-portable-trend-20261015`. Exact source, binary, meter,
  orchestration, public and checker-only input hashes are in
  [HELDOUT_NINE_ARM_SPEC_V4.json](HELDOUT_NINE_ARM_SPEC_V4.json).
- **Hypothesis, frozen targets, controls and stop conditions:**
  [HELDOUT_NINE_ARM_PROTOCOL_V4.md](HELDOUT_NINE_ARM_PROTOCOL_V4.md).
  The corrected seven-arm V3 and its exhaustive S3 check remain frozen.
- **Commands:** `bash RUN_HELDOUT_NINE_ARM_V4.sh` launched 108 cold
  native children with `different_j_chi_heldout_meter`; outputs are
  under [heldout_nine_arm_v4](heldout_nine_arm_v4/). The separate
  `different_j_chi_heldout_check PUBLIC CHECKER DIRECTORY` replayed
  all children; `different_j_chi_heldout_summary DIRECTORY` reduced
  them and asserted identical six-arm relation support.
- **Preserved implementation correction:** first compile of the new
  inverse held-out output hit Rust's `json!` macro recursion limit.
  Splitting the nested relation-setup object into a separate value
  fixed the build before any V4 child ran. The source and binary hashes
  in the spec are the successful build; the V3 and four-order source
  snapshots remain byte identical to their earlier specs.
- **Limits and outcome:** 300-second hard wall timeout per child,
  512 MiB observed peak RSS on Darwin. No V4 child failed or timed out.
  The independent receipt reports 108 verified scalars and 996
  relation witnesses. Every inverse and full-image base equalled the
  frozen eight columns. [V4 results](HELDOUT_RESULTS_V4.md),
  [summary](HELDOUT_NINE_ARM_V4.summary.json),
  [checker receipt](HELDOUT_NINE_ARM_V4.check.json), and
  [graph](HELDOUT_GRAPH_V4.svg) contain the measured evidence.
