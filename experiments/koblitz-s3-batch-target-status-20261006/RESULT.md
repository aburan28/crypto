# Failed compact S3-batch targets now fail the process

The [frozen protocol](PROTOCOL.md) was committed as `fa25c1ea2` and opened
in PR #1475 before either release control ran or the producer source changed.
PR #1475 later merged the protocol alone; this follow-on PR carries the
producer correction, independent verifier, and measured control artifacts.
The n41/K20 public point reproduced the contract defect on exact old main
`5874ce8b3`: the producer wrote `targets_failed:1`, no target scalar and
`group_verified:null`, yet wrote target `exit_code:0` and returned process
status **0**. The corrected producer returned process status **1** and target
`exit_code:1` on that same point. The old and corrected base and complete
rank traces are **byte-identical**, including 43 attempts, 23 misses and
full rank 20. The source edit changes only status reporting after each
target's timed interval and the final process return.

| Frozen control | K / actual base points | Rank attempts / misses | Target scalar | Target `group_verified` | Target `exit_code` | Process exit | Independent replay |
| --- | ---: | ---: | --- | --- | ---: | ---: | --- |
| Old source, failed pilot Q | 20 / 1,640 | 43 / 23 | null | null | **0** | **0** | Same base/rank and target decision as fixed arm |
| Corrected source, same failed Q | 20 / 1,640 | 43 / 23 | null | null | **1** | **1** | General-law rank and base-log replay passed |
| Corrected source, successful holdout Q | 85 / 6,970 | 85 / 0 | 449905522971 | true | **0** | **0** | General-law rank, four-point witness and scalar replay passed |

The Rust [independent verifier](../../examples/koblitz_s3_batch_status_replay.rs)
recomputed all 105 accepted rank relations and modular log systems, checked
all 8,610 factor-base point logs with the general binary-curve law, and
checked the successful four-point sum and `[449905522971]G=Q`. It also
checked that both failed arms used the same public point and had identical
base/rank bytes and target decision. Its compact [receipt](REPLAY.json) is
`verified:true`. The failed target has no scalar to replay and remains a
failed run. These controls do not estimate natural PDP yield or support, and
no CPU timing or IC/rho ratio is claimed.

The unit control covers `Some(true)`, `Some(false)` and `None` status, and
all seven S3-batch example tests passed. Local Clippy reported no
example-specific findings after the verifier cleanup; warnings in the
existing library remain. Both release producer builds used
the same offline [Cargo lock](producer-Cargo.lock) on Darwin 25.6.0 arm64,
Homebrew `rustc`/`cargo` 1.93.1. The old source SHA-256 was
`7b2393731cd7da047b9404585c85dc2f384ee857e43778661e806649c6cf552d`,
old binary `ffa6a5b3d9401b1f07cae5ccd5b4081ccd902ec0a8bd2d3fec163da773a1601a`,
corrected source `f81d4c3e4028d7c6a67141ff84d9797cfbfacf5433a232f853537b430832927e`,
and corrected binary
`d0035d0d634fe611be6bc4d3ff058a4b3f9e75df9520b9b842a0eb3749a73182`.
The [archived old producer source](old-producer.rs) is the exact mainline
input. The verifier source SHA-256 was
`e898eff13584be2849551ad94d5b93f3329bff13b6bd789e08bc20b7c852f66b`
and its release binary SHA-256 was
`fcdb59edd4302b111b412e1d85a3470d0ab2160ccd4e494bcf87ba77e0fa3982`.
The first local baseline build failed because the newly created sparse
worktree omitted four checked-in `docs/` files used by Rust `include_str!`;
those paths were materialized, and the unchanged old source then built.
This was a build-environment failure before any producer run.

The exact producer invocation for each row was the following command with
`K`, point file and output directory from the table above, built once from
each source state. `KIC_DUMP_BASE` and `KIC_DUMP_RANK` wrote the archived
base and rank files under that output directory:

```sh
env RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1 LC_ALL=C \
  KIC_S3_BATCH_WINDOW=64 KIC_S3_PREFILTER=blocked KIC_QUERY_BACKEND=point_sum \
  KIC_DUMP_BASE=OUT/base.jsonl KIC_DUMP_RANK=OUT/rank.jsonl \
  timeout 45 target/release/examples/koblitz_orbit_dlp_s3_batch \
  construct:41:0:K POINT.jsonl 410041 OUT/targets.jsonl \
  > OUT/summary.jsonl 2> OUT/stderr.txt
```

The archived `runs/*/process_exit.txt` files record the shell results 0,
1 and 0. The independent native replay is reproducible with:

```sh
cargo build --offline --locked --release --example koblitz_s3_batch_status_replay
target/release/examples/koblitz_s3_batch_status_replay \
  experiments/koblitz-s3-batch-target-status-20261006 \
  experiments/koblitz-s3-batch-target-status-20261006/REPLAY.json
shasum -a 256 -c experiments/koblitz-s3-batch-target-status-20261006/RAW_SHA256SUMS --quiet
```

The 45-second cap was not reached by any control. Raw point, target, base,
rank, summary, stderr, process-status, source, lock and replay hashes are
retained in [`RAW_SHA256SUMS`](RAW_SHA256SUMS). This correction is required
before a fresh public-Q full-point pair-sum evaluation can preserve failures
honestly. That evaluation still needs its own frozen workload, independent
replay and same-Q rho comparison; this PR does not supply those results.
