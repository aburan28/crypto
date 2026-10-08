# Reproducing the n37 whole-solve instruction census

The [protocol](PROTOCOL.md) was pushed before measurement. The
[decision](RESULT.md), [machine-readable result](DECISION.json),
[frozen jobs](JOBS.json), [compact raw census](CENSUS.json), and
[full CI archive](RAW-CALLGRIND.zip) are tracked in this PR. The original
host/toolchain/binary receipt is [PROVENANCE.txt](PROVENANCE.txt).

Recompute the decision from the committed raw Callgrind parts using only
native Rust and an archive extractor:

```sh
mkdir -p /tmp/ecbench-callgrind-replay
unzip -q research/ecbench_callgrind_solve_20261004/RAW-CALLGRIND.zip -d /tmp/ecbench-callgrind-replay
cargo build --release --example ecbench_callgrind_analyze
target/release/examples/ecbench_callgrind_analyze \
  /tmp/ecbench-callgrind-replay \
  research/ecbench_callgrind_solve_20261004/JOBS.json \
  research/ecbench_n37_rank_columns_20261004/sessions/mac_arm64_l0_01/records.jsonl \
  /tmp/ecbench-callgrind-decision.json
cmp /tmp/ecbench-callgrind-decision.json research/ecbench_callgrind_solve_20261004/DECISION.json
```

The analyzer refuses a different archived session hash, missing or duplicate
jobs, any wrong scalar or method, a failed profile, a raw-file hash mismatch,
or a discrepancy between Callgrind parts and the compact census. The
Callgrind unit is **simulated user-space instructions** over one complete
`methods::solve` call per public point; no CPU wall-time claim follows from
this reproduction.
