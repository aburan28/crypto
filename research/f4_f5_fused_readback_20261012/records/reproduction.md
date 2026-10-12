# Native execution and report reproduction

All crypto execution, correctness checks, profile parsing and result arithmetic
are Rust. Shell only builds/runs. The Python script renders recorded native
fields into SVG/Typst/PDF; it is not a benchmark harness or result calculator.

Production baseline `c2963dcf`, first candidate `07151ff0`, second candidate
`07e568b1`, corrected production `d0a1d2d8`, and final analyzer `c8e1bcb6` are
recorded separately. The later `05c92eb4` repair changes only a curve test under
`#[cfg(test)]`. Every measured driver uses the identical frozen source
`examples/f4_f5_row_build_bench.rs`. Reproduction requires the recorded compiler,
locked Cargo dependencies and profiler; a different host/compiler is a new study.

The executed build/test commands, with the workspace toolchain on PATH, include:

```sh
cargo build --locked --release --example f4_f5_row_build_bench -j2
timeout -k 5 1200 cargo test --locked --release --lib cryptanalysis::koblitz_groebner::tests -j2 -- --test-threads=3
cargo test --locked --release --lib cryptanalysis::matrix_f5_f2::tests -j2 -- --test-threads=3
cargo test --locked --release --lib cryptanalysis::pq_groebner_f2::tests -j2 -- --test-threads=3
timeout -k 5 1200 cargo test --locked --release --lib -j2 -- --test-threads=3
timeout -k 5 1200 cargo test --locked --release --lib cryptanalysis::curve_traits::tests -j2 -- --test-threads=3
cargo test --locked --release --example f4_f5_hash_compare --example f4_f5_row_build_bench -j2
```

The baseline binary was copied before any production change to
`/workspace/scratch/f4-f5-fused-readback/a/bench`. Each candidate was copied to
the equally long `b/bench` path; prior binaries were saved in scratch before
replacement. SHA-256 and committed source revisions permit rebuilding them.
The thin `run-profiles.sh` refuses existing destinations and retains timestamps,
commands, exit codes, JSONL and all client-request dumps. It uses the locally
extracted Valgrind package and `VALGRIND_LIB`; package hash/version are recorded.

```sh
bash research/f4_f5_fused_readback_20261012/run-profiles.sh baseline
bash research/f4_f5_fused_readback_20261012/run-profiles.sh
bash research/f4_f5_fused_readback_20261012/run-profiles.sh legacy
```

For later cohorts, set `F4_F5_RECORD_ROOT` to their fresh record directory;
reuse the unchanged primary/A/A baseline records byte for byte and run new
candidate primary and all five paired processes. Rotated process order and
the one-worker Rust driver remain fixed. Corrected `legacy` mode records both
`KIC_F5_FUSED_BUILD=0` and `KIC_F4_HASH_COLUMN_INDEX=1`; the first candidate's
legacy invocation only used the F5 flag (the new F4 flag did not exist then).
Shared baseline records are not counted as independent new samples.

The native helper was compiled directly while Cargo held its build lock,
against the same already compiled locked JSON/BLAKE3 artifacts:

```sh
clippy-driver --edition=2021 -O -L dependency=target/release/deps \
  --extern serde_json=target/release/deps/libserde_json-d00c0fe8816b70e5.rlib \
  --extern blake3=target/release/deps/libblake3-e60f6d89c18cbd96.rlib \
  -D warnings examples/f4_f5_hash_compare.rs -o /workspace/scratch/f4-f5-fused-readback/compare
/workspace/scratch/f4-f5-fused-readback/compare research/f4_f5_fused_readback_20261012/records/attempt-2-fixed
/workspace/scratch/f4-f5-fused-readback/compare check research/f4_f5_fused_readback_20261012/records/attempt-2-fixed primary-a 3 legacy-b 1
```

`inspect RECORDS_DIRECTORY RUN REPEATS [f4|f5|criterion]` writes worker-inclusive
top self-costs; its default retains the old criterion fields. `check` validates
two runs without admitting any instruction/memory verdict. The original
all-33-case gate and the explicit complete-call gate coexist; all phase
regressions are retained. The native JSON recomputation was byte-identical.
An initial type-inference compile error and a version query without the
required `VALGRIND_LIB` were preserved and corrected.

The complete fixture repair is `844cb756`. The initial library failure and
the first incomplete fixture repair both remain recorded; the completed repair
passed 24 native coordinate-module tests, eight focused Cargo tests and all
required standard-library/example checks. Standard-suite ignores come from
the recorded current-main baseline. The exact
non-weak isogeny check is executed separately in a copied release-test binary:

```sh
timeout -k 5 1200 /workspace/scratch/f4-f5-fused-readback/jv-validation-bin \
  cryptanalysis::jv_isogeny_walk::tests::end_to_end_from_a_non_weak_curve_recovers_and_verifies \
  --exact --ignored --test-threads=1
```

Its source was not changed; its copied binary predates only the independent
coordinate-fixture repair. The command, binary hash, exact status and log stay
under `verification-final/`. Failed and author-interrupted earlier library
runs are never replaced by final receipts.

Report rendering uses the existing pinned scratch virtualenv:

```sh
/tmp/f4-f5-report-venv/bin/python research/f4_f5_fused_readback_20261012/render.py \
  --out research/f4_f5_fused_readback_20261012/reports/run-20261012T0543Z
```

Render a second time to scratch and compare every generated SVG, Typst source
and PDF byte for byte. Rasterize all six combined pages and all six one-page
PDFs with `pdftoppm`, inspect every page, and retain the final inspection/rebuild
receipt. The final validation strings come from actual completed native checks.
The earlier draft was labeled pending and is replaced before publication.

The normal fixture-only rebuild has a different whole-file hash; preserve its
failed byte comparison. The native `verify_elf.rs` checker compares loader
headers and every allocated section's layout/content, skipping only GNU
build-ID contents and zero-initialized NOBITS contents. Its successful receipt
confirms that runtime code/data and layout are unchanged. The differing symbol
metadata is not promoted as a code difference or an identical file. The original
full byte-difference listing is retained compressed, with its uncompressed hash.
