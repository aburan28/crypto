# Validation commands and evidence

Commands run from the repository root, with cargo/rustfmt from Rust 1.99.1.
The local installation uses `CARGO_HOME=/workspace/.sea-cargo` and
`RUSTUP_HOME=/workspace/.sea-rustup`. These paths are environment settings,
not dependencies of the sources.

```sh
cargo test --release --manifest-path isogeny_algos/Cargo.toml --lib --test sea --test api --test kernel
cargo test --release --manifest-path isogeny_algos/Cargo.toml --example gls_loop_validate --example gls_loop_bench
cargo clippy --release --manifest-path isogeny_algos/Cargo.toml --lib -- -D warnings
cargo test --release --manifest-path research/sea_point_count_20261016/native/Cargo.toml --bin analyze
cargo clippy --release --manifest-path research/sea_point_count_20261016/native/Cargo.toml --bins -- -D warnings
cargo run --release --manifest-path research/sea_point_count_20261016/native/Cargo.toml --bin cm_regression
```

The release suite passes 23 library, seven API, three kernel and eight SEA
tests; one existing library test remains ignored. Both CM example test targets
compile (they contain zero unit tests). The analyzer's rejection unit test
passes. Library and native Clippy checks pass. Logs are under `evidence/`,
including `rust-validation-final.txt`, `examples-validation.txt`, `clippy.txt`,
`native-clippy.txt` and `analyzer-validation.txt`.

The CM regression was run before and after the implementation fix using the
same native source. The old constructor at 44f8d17f yielded 17 failures
(exit 1); the corrected constructor yielded zero (exit 0). Raw JSONL and exits
are retained separately. The valid cases comprise both trace signs for
(p,D,v)=(11,-7,2) and (59,-23,2), eight seeds per sign. The final unit suite
also checks invalid inputs and certificate-based twist selection at 40 bits.

Touched files pass rustfmt --check: gls_loop.rs, find/sea.rs, curve.rs and
the four native binaries. The repository-wide crate format gate has an
unrelated pre-existing lib.rs module-order failure, also present at cf274521
and origin/main. It is reported in the PR rather than patched around.

```sh
TLA_JAR=/workspace/sea-tools/tla2tools-1.8.0.jar \
  FORMAL_LOG_DIR="$PWD/research/sea_point_count_20261016/formal/logs/final" \
  bash research/sea_point_count_20261016/formal/verify.sh
```

For the recorded run, Java 21 and the Lean 4.19.0 toolchain were on PATH
(`ELAN_HOME=/workspace/.sea-elan`). Formal exit codes are 0, 0, 12, 0, 12
for Fixed, Wide, Legacy, Certification and SampleOnly respectively. The
script requires the expected invariant in each negative control. Final
Lean build and ten axiom inspections pass. See VERIFICATION.md for the
proof map, tools, bounds and implementation obligations. For a fresh replay,
choose a new FORMAL_LOG_DIR; do not overwrite this frozen evidence.

Report rendering uses only existing native summaries:
`MPLCONFIGDIR=/tmp/sea-matplotlib XDG_CACHE_HOME=/tmp/sea-font-cache python3 research/sea_point_count_20261016/render_reports.py`.
The results PDF (three pages) and run PDF (two pages) were rendered with
pdftoppm and visually inspected, including plots, table continuations, source
references and verification scope. SVG/150 dpi PNG and browsing copies are
generated from the same rendering source. See README.md for benchmark replay.
