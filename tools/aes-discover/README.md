# Native AES differential discovery laboratory

This integration mirrors the native sources, tests, example and lockfile
from the companion `aburan28/cryptanalysis` AES laboratory. Source parity should
be checked before applying future changes. It is an independent standalone Rust
crate; it does not replace existing AES research code or official campaigns.

Build and verify from the repository root:

```sh
cargo build --release --locked --manifest-path tools/aes-discover/Cargo.toml
cargo test --release --locked --manifest-path tools/aes-discover/Cargo.toml
cargo clippy --release --locked --all-targets --manifest-path tools/aes-discover/Cargo.toml -- -D warnings
cargo test --release --lib
```

Rust 1.87+ and OpenSSL for independent correctness fixtures are required.
Run `tools/aes-discover/target/release/aes-discover --help` for commands.
It supplies differential endpoint discovery with frozen, hash-bound validation,
exact DDT/LAT/BCT, bounded trail proposals, integral and local quartet controls,
key-schedule and fixed-delta boomerang diagnostics, and key-disjoint datasets.
Every new output path is write-once. Seeds describe public synthetic fixtures.
Use `bash tools/aes-discover/check-mirror.sh /path/to/cryptanalysis` to compare
the sources, tests, example and Cargo files with the companion checkout.

Validation records each fresh fixed key, disjoint pairs, a matched lazy PRP,
conservative Chernoff/finite-PRP coupling bounds and family correction. Candidate
acceptance requires the frozen direction on every tested key and passing matched
controls. Discovery scores and DDT-product trail scores are proposal diagnostics.
The reported sensitivity is an IID-model diagnostic, not measured power.

This is table-based research AES, not production encryption. No full-AES break,
SAT certificate backend, exhaustive hull solver, trained model, GPU kernel or
complete published quartet/key-recovery attack is included. The original frozen
suite, API documentation, observations and review record are in the companion:
https://github.com/aburan28/cryptanalysis/tree/codex/aes-distinguisher-research/tools/aes-discover
https://github.com/aburan28/cryptanalysis/blob/codex/aes-distinguisher-research/experiments/aes-distinguisher-lab-20261017/RESULTS.md

This task's standalone integration is explicitly authorized by the current
user. The Conductor service remains unavailable and no reservation is claimed.
See the local [integration record](../../research/aes_distinguisher_20261017/INTEGRATION.md)
and the companion access record:
https://github.com/aburan28/crypto-autoresearcher/blob/codex/aes-distinguisher-research/problems/aes-lab-conductor-and-campaign-admission.md
Software integration does not grant official
scientific campaign admission or establish novelty.
