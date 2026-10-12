# `ecc2k130-dp-ingest`

Standalone Rust port of the offline half of [`aws/dp_ingest.py`](../aws/dp_ingest.py).

This crate is **outside** the root `crypto` package on purpose: the IC autolab
jobs freeze `research/ic_candidate_tournament_20260915/ci/Cargo.lock` onto the
root, and ingest needs Postgres + AWS SDK crates that are not on that lock.
Its own `Cargo.lock` lives here.

## Phase A (this tree)

- Record decode, key shapes, envelope / table-v3 / WITNESS=1 checks
- Coverage counting and the weight-32 cutoff verdict
- Status helpers (`state`, window sums, walk rate, banned-field guard)
- CLI stub: same flags as the Python, plus `--decode-file` for offline use
- Live daemon: still `aws/dp_ingest.py`

Protocol: [`research/ecc2k130_dp_ingest_rust_20261007/PROTOCOL.md`](../../research/ecc2k130_dp_ingest_rust_20261007/PROTOCOL.md).

```sh
cargo test --manifest-path ecc2k130/dp-ingest/Cargo.toml
cargo clippy --manifest-path ecc2k130/dp-ingest/Cargo.toml --all-targets -- -D warnings
cargo run --manifest-path ecc2k130/dp-ingest/Cargo.toml -- --decode-file some.bin
```

## Later

- **B:** Postgres/S3 adapters and a parity harness against Python on frozen fixtures
- **C:** cutover of systemd / Lambda / `ingest.sh`; keep Python for rollback
