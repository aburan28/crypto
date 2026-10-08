# Worker source validation

This validates the new ordinary-query worker and capsule contract as source.
There is no frozen worker/controller, accepted scientific registration, native
panel, new natural-yield result, one-target solve or performance claim yet.
The full goal and all remaining protocol gates stay open.

All 11 worker/transport controls passed locally, zero failures. Four controls
exercise the new strict configuration and compiled identity gate, chronological
create-only evidence, new claim/parent/durable-query binding, and fixed accepted
native argv. Seven existing shared transport controls cover immutable file
checks, accepted asset extraction without execution, bounded export reads,
registered stdin/environment delivery, nonreading-worker deadlines and process
group termination. These latter controls launch small shell fixtures only;
none dispatches F5, an exporter or CryptoMiniSat search.

The final source command exited zero:

```sh
/Volumes/SSD990/llm/tmp/ic-native-busy-20261002 busy -- sh -c '
  rustfmt --edition 2021 src/bin/prepared_ordinary_worker.rs &&
  CARGO_INCREMENTAL=0 CARGO_TARGET_DIR=/Volumes/SSD990/llm/tmp/ic-native-preparation-build-cache-20261004 cargo test --locked --offline --bin prepared_ordinary_worker -- --test-threads=1 &&
  CARGO_INCREMENTAL=0 CARGO_TARGET_DIR=/Volumes/SSD990/llm/tmp/ic-native-preparation-build-cache-20261004 cargo clippy --locked --offline --bin prepared_ordinary_worker
'
```

The earlier v1 and v2 logs are retained, not replaced; they covered three and
four new contract controls respectively before the final shared-transport run.
Decompress `local-tests-v3-20261004.txt.gz` with `gzip -cd` for all 11 controls.
Clippy emitted no changed-source warning but reported inherited warnings and
an unknown configured newer lint on local Rust 1.93.1. Its zero exit is not a
strict clean Clippy pass. The repository's pinned CI toolchain and exact-head
strict checks remain the acceptance gate. Linux/macOS CI runs the worker tests.
Formatting and `git diff --check` passed locally.

The actual ordinary debug binary was built through the busy launcher and its
`build-identity` command executed once. It returned null source/build pins, as
retained in `unfrozen-worker-identity.json`. The identity gate rejects such a
build. This identity-only invocation performs no preparation or solver work.
The raw build output and hashes of the source, lockfile, logs and unfrozen
binary are retained alongside it. This is not a scientific build freeze or a
source-runtime execution receipt. Incremental builds were disabled to limit
the disposable cache footprint on the nearly full shared volume.

The worker now writes its producer report before surfacing post-run drain or
immutable-binding failures. It retains a separate terminal diagnostic even if
the producer returns an error. An outer kill still requires the future frozen
controller's terminal receipt and audit of the durable prefix; these unit
controls do not prove that complete interrupted-panel integration.

No historical capsule was dispatched or restored. No IC1 identity, curve,
ecbench session, measured ratio or tournament result is introduced. Autolab
status and skill links change; no scoreboard winner or leaderboard result is
added. This PR is dependent, unmerged source until its exact-head validation
and stack gates are resolved.
