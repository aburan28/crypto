# Compact-orbit target failure exit contract

The frozen n41 factor-base-size pilot in PR #1353 observed three `K=20`
processes that reported `targets_failed: 1` but exited 0 and wrote target
rows with `exit_code: 0`. The independent replayer correctly rejected all
three. A caller that trusts only the process status, however, could accept
an unsolved logarithm. This focused fix makes each target row's `exit_code`
reflect scalar verification and makes the process return failure if **any**
target in its input fails, after writing the summary, base dump and rank
trace. Successful targets and successful processes retain exit code 0.

[`input/n41_public_q.jsonl`](input/n41_public_q.jsonl) is the exact public
pilot point `[1119268120096,860021860491]` from that frozen experiment.
The no-overwrite regression script exercises both an unsolved `K=20` and
solved `K=85` run on this same point, asserts the process and target-row
statuses, and checks that failure evidence remains available:

```sh
cargo build --offline --locked --release --example koblitz_orbit_dlp_fast_online
experiments/koblitz-failed-target-exit-20261004/verify.sh target/release/examples
```

The test has a 30-second cap per process. Its temporary outputs are erased
only after assertions; the original PR #1353 raw failure rows remain
immutable, with their old source/binary hashes. This change corrects the
interface contract and does not reinterpret or replace those measurements.
