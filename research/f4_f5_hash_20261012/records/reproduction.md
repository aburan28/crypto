# Native execution and validation commands

Production/library checks used the locked root Cargo dependencies and recorded Rust 1.99 toolchain. Baseline and candidate binaries use the identical `examples/f4_f5_row_build_bench.rs`; their source and binary hashes are retained. Callgrind 3.24 was extracted from the Debian package whose SHA-256 is recorded here. No package binary or credential is committed.

```sh
cargo build --locked --release --example f4_f5_row_build_bench --example f4_f5_hash_compare -j 2
cargo test --locked --release --lib cryptanalysis::matrix_f5_f2::tests -j 2 -- --test-threads=3
timeout -k 5 1200 cargo test --locked --release --lib -j 2 -- --test-threads=3
cargo test --locked --release --example f4_f5_row_build_bench --example f4_f5_hash_compare -j 2
```

`run-profiles.sh` records every measurement command and exit. The baseline primary was launched before that orchestration script existed, using the same Callgrind flags and the saved `a/bench profile 3` binary. Its exact start timestamp was not captured. All 18 process runs completed successfully; no timed host comparison was attempted. `comparison.json` hashes every used dump and native-record file.

During the library check, Cargo holds its build-directory lock. The comparison example was also built, tested and linted directly against the already compiled locked dependencies; this does not affect the measured driver binaries:

```sh
clippy-driver --edition=2021 -O -L dependency=target/release/deps \
  --extern serde_json=target/release/deps/libserde_json-d00c0fe8816b70e5.rlib \
  --extern blake3=target/release/deps/libblake3-e60f6d89c18cbd96.rlib \
  -D warnings examples/f4_f5_hash_compare.rs -o /workspace/scratch/f4-f5-hash/compare
```

The driver lint uses its recorded library, Rayon, JSON and BLAKE3 artifacts and also passed. Five comparison tests passed using both the native test build and the Cargo example target. Clippy's initial style findings and the failed inspector/compiler attempts remain beside the successful logs. After cleanup, the accepted comparison matched the previous output byte for byte.

PDFs were rasterized with `pdftoppm`; all six study pages and the refreshed planning-snapshot page were opened with `view_image`. The final memory caption was shortened before reinspection. `render.py` regenerated all study outputs into a separate directory; `diff -qr` returned zero. It runs no cryptographic computation or result arithmetic.

Named library coverage is checked against the actual test binary's `--list` catalog. Ignored lines with reasons use `ignored, ...`, so coverage parsing checks that status field rather than requiring `ignored` as the final token. The final named sets cover all 4,016 catalog entries, matching the successful runner summary; none is unfinished.
