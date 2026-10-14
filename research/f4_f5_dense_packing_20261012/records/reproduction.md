# Native reproduction

Use the recorded baseline and production candidate revisions with the exact
common driver, Cargo.lock, Rust compiler and baseline release flags. The
driver is byte-identical to the previous study. Scratch binaries were built
and retained at equally long `a/bench` and `b/bench` paths; the committed source
revisions and hashes permit rebuilding them. Binary files stay in scratch.

```sh
export PATH=/workspace/toolchains/cargo/bin:$PATH
export RUSTUP_HOME=/workspace/toolchains/rustup
export CARGO_HOME=/workspace/toolchains/cargo
cargo build --locked --release --example f4_f5_row_build_bench -j2
# Save this output at the appropriate a/bench or b/bench path before changing code.
bash research/f4_f5_dense_packing_20261012/run-profiles.sh baseline
bash research/f4_f5_dense_packing_20261012/run-profiles.sh
bash research/f4_f5_dense_packing_20261012/run-profiles.sh hash-control
bash research/f4_f5_dense_packing_20261012/verify.sh
```

The runner requires Valgrind 3.24 at its recorded extracted path and
`VALGRIND_LIB`, includes worker instructions and refuses existing output paths.
Native-time isolation was not earned, so no timings are compared.
The verification commands require permission to create loopback sockets;
the default sandbox's nine mock-network failures and network-enabled repeat
remain separate evidence in `VERIFICATION-ENVIRONMENT.md`.

Native analysis was compiled directly while Cargo built the library:

```sh
clippy-driver --edition=2021 -O -L dependency=target/release/deps \
  --extern serde_json=target/release/deps/libserde_json-d00c0fe8816b70e5.rlib \
  --extern blake3=target/release/deps/libblake3-e60f6d89c18cbd96.rlib \
  -D warnings examples/f4_f5_hash_compare.rs -o /workspace/scratch/f4-f5-dense-packing/compare
/workspace/scratch/f4-f5-dense-packing/compare research/f4_f5_dense_packing_20261012/records
/workspace/scratch/f4-f5-dense-packing/compare check research/f4_f5_dense_packing_20261012/records primary-a 3 hash-control-b 1
/workspace/scratch/f4-f5-dense-packing/compare inspect research/f4_f5_dense_packing_20261012/records primary-a 3 f5
/workspace/scratch/f4-f5-dense-packing/compare inspect research/f4_f5_dense_packing_20261012/records primary-b 3 f5
```

The statistics, profile parser and math-output validation remain native Rust.
The comparator adds one acceptance field for this preregistered F5-only change;
all existing fields, result arithmetic, bootstrap, workloads and dump-coverage
requirements stay intact. Cargo also runs both relevant examples' tests.

The test-only fixture-size adjustment changes no production logic. A final
normal benchmark rebuild is checked against the measured candidate by the
previous study's frozen native ELF helper. All allocated runtime sections and
loader layouts must match, ignoring build-ID contents only; raw whole-file
hashes are preserved separately. Verification logs and receipt are in records.

Library coverage is checked by matching the complete test catalog against
completed names, retaining ignored names and unmatched/unfinished lists.
Rust's ` - should panic` display suffix is stripped for name matching only;
the raw labels, first unnormalized comparison and original logs are retained.

Rendering only, using pinned `scripts/perf/groebner-render-requirements.txt`:

```sh
/tmp/f4-f5-report-venv/bin/python research/f4_f5_dense_packing_20261012/render.py \
  --out research/f4_f5_dense_packing_20261012/reports/run-20261012T1503Z
/tmp/f4-f5-report-venv/bin/python research/f4_f5_dense_packing_20261012/render-prior.py \
  --out /workspace/scratch/f4-f5-dense-packing/prior-visuals
```

Run the current renderer again in a fresh scratch directory, compare SVG/PNG,
Typst and PDF bytes, inspect every group and combined PDF page with `pdftoppm`,
and retain artifact hashes and the inspection receipt. Copy all graphs and
report PDFs to the corresponding `benchmarks/` entry; the study remains the
source of truth. The prior study's frozen SVG/PDF bytes must remain unchanged.
