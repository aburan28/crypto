# Compact-orbit producer merge repair: correctness replay

The repaired `koblitz_orbit_dlp_fast` example parses, builds, and emits a rank trace that the existing independent Python verifier accepts. The repair removes duplicated statements introduced by union merge `fe07f6088f35dd2150c960ef2a0ca14ae0d11641` while retaining its rank-trace addition and the first parent's wide-field path.

This is a correctness smoke on the repository's n=13 public fixture, not a CPU performance comparison. The full n=53 fixed-base restart specified in the cryptanalysis repository remains a separate experiment.

## Frozen input and execution

| Item | Value |
| --- | --- |
| Parent revision | `845fd323015e28d32b73fa25e551f830d5826e6c` |
| Repaired example source SHA-256 | `627b18beff7352fbd092fb284412e67f5fb9ce8a5ad1f3837378795d742e0eb7` |
| Release executable SHA-256 | `74ff60142cc90fcfb92281d0247e5ee75f000978ddd80b8ef16d4714eeaab6a0` |
| Target input SHA-256 | `10159baf262b43a92d95db59dae1f72c645127301661e0a3ce4e38b295a97c58` |
| Host and toolchain | Darwin arm64; Homebrew `rustc` and `cargo` 1.93.1 |
| Construction | `construct:13:0:2`, one known-answer scalar `7`, rank seed `7` |
| Curve and base | GF(2^13), Koblitz a=0, subgroup order 2003; 52 points and 2 orbit columns |

The source hashes are of the worktree file and executable used for this run. From the repository root, build with `CARGO_INCREMENTAL=0 cargo build --offline --release --example koblitz_orbit_dlp_fast`. Then, from this evidence directory, run:

```sh
KIC_DUMP_BASE=base.jsonl KIC_DUMP_RANK=rank.jsonl \
  ../../../../target/release/examples/koblitz_orbit_dlp_fast \
  construct:13:0:2 targets_n13.txt 7 targets.jsonl > summary.jsonl
python3 ../compact_orbit_rank_evidence_20260929/verify_rank.py \
  --trace rank.jsonl --base base.jsonl --summary summary.jsonl --out replay.json
```

The command above shows the intended path layout; the recorded run used an isolated `CARGO_TARGET_DIR` and absolute output paths. `verify_rank.py` resolves its independent group arithmetic from the repository root. `summary.jsonl` is producer stdout; `targets.jsonl` is the path passed as the fourth producer argument.

## Result and validation

The producer inserted 2 verified relations in 2 attempts and reached rank 2/2. The independent replay checked both representative logarithms and 7 used base points (`replay.json`: `status=PASS`). The single target record reports recovered scalar 7, matching its published fixture scalar, with group replay true. The 1,029-byte summary retains raw stage timing, which is an uncontrolled correctness diagnostic on a contended host.

`rustfmt --edition 2021 --emit stdout` parsed the repaired source; `git diff --check` passed. `cargo test --offline --release --example koblitz_orbit_dlp_fast` passed all 3 wide JSON tests. The compact-orbit archived n=37 source42 library test passed after its tracked fixture was materialized in the sparse worktree. A full release library test was attempted; its unrelated failures and exact test environment are recorded in the PR validation report.

The first replay invocation incorrectly supplied the target-record file as `--summary`, yielding `replay_input_error.json`. Its input is preserved as `targets_attempt1.jsonl`; the corrected stdout capture produced the passing replay above. This was a command-routing error, not a producer or verifier arithmetic failure.
