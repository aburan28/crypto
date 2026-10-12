# Local integration smoke before the hosted cold panel

Status: **smoke only; not timing evidence**. The PR's frozen public Q and
`d8d67d3b683334db105195b0fc74d8e660dc4b5e` source were used. The
binary was built with the frozen Cargo lock using `cargo build --release
--locked --offline` for the two named examples. The host was macOS, with
Homebrew `rustc 1.93.1`; it did not have the Linux reserved-core monitor.
No local CPU number is a cold-panel result or a rho comparison.

The first exact-source worktree preflight did not materialize its sparse
checkout after `git worktree add --no-checkout`, so Cargo reported that
`Cargo.toml` was missing. `git read-tree -mu HEAD` populated the configured
sparse paths; the resulting worktree had zero tracked changes and built
successfully. No solver or timed arm ran before that repair.

| Raw directory in archive | Producer | Independent result | Scope |
|:--|:--|:--|:--|
| `n37_L1_v1` | exit 0 | `SMOKE_PASS` under the earlier runner | Preserved development smoke; runner hash predates final metadata fields |
| `n37_L1_v2` | exit 0 | `SMOKE_PASS` | One new Q, full rank, one recovered log |
| `n41_L1024` | exit 0 | `SMOKE_PASS` | Full rank and 1,024 independently checked logs |
| `n53_L1024` | exit 0 | `SMOKE_PASS` | Full rank and 1,024 independently checked logs |

`LOCAL_SMOKE.tar.gz` contains all 28 raw and receipt files, including
stdout, stderr, base, rank, target and run metadata for every smoke. Its
SHA-256 is
`2dac40e561aba42724a9aae085855b5ee32e5e8bbb2fd6fee316c1c4f686ff4e`
(1,670,556 bytes). `LOCAL_SMOKE_FILES.json` lists every member's byte count
and hash. To replay on a checkout containing the frozen inputs, extract with
`tar -xzf research/notes/ecc2k130/point_sum_cold_20260930/LOCAL_SMOKE.tar.gz
-C /tmp/point-sum-smoke` after making that directory, then run
`verify_cold.py --cell n37_L1 --run-dir /tmp/point-sum-smoke/n37_L1_v2
--out /tmp/point-sum-smoke/n37_L1_v2/replay.json --relocated` (or the
corresponding n41/n53 cell). Use a new output path each time. The v1 trace
is retained as raw history and is not a current-runner replay target.

The relocation check on an unchanged copy of `n37_L1_v2` returned
`SMOKE_PASS`. A separate copy with one byte of the raw target JSON changed
was rejected by the verifier's SHA-256 check with exit status 1;
`TAMPER_REJECTION.json` preserves that expected failure receipt. These
checks exercise the artifact boundary, not the timing eligibility gate.
The six Linux cold cells are still pending the hosted workflow, A/A and
isolation gates, independent second-host replay, and the outcome PR.
