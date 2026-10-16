# Validation receipts

| Gate | Result | Raw receipt |
|---|---|---|
| Root `cargo test --release --lib --locked --offline -- --test-threads=2` | 4,083 passed, zero failed, 117 existing ignored; 1,347.49 seconds test interval | `evidence/root-release.log` |
| Standalone release library/binaries/integration/examples | PASS, one existing ignored exploratory library test; includes direct-division threshold checks and pivoted known solutions | `evidence/algorithm-release.log` |
| Native CLI and parser tests | 3 passed, zero failed | `evidence/perf-cli-tests.log` |
| Standalone `cargo clippy --all-targets -- -D warnings` | PASS | `evidence/clippy-final.log` |
| `cargo fmt --check` | PASS | `evidence/format-final.log` |
| Workflow lint with the documented Intel runner label configured | PASS | `evidence/actionlint.log` |
| Baseline/rung 1/2/3 constructor equality | 96/96 at each rung; separate incremental comparisons | `evidence/rung*/COMPARISON.json` |
| Final repeated Apple M4 Pro campaign | 20/20 valid; outputs equal; counts repeatable; zero qualified wall samples | `evidence/apple-m4-final-repeats/SUMMARY.json` |

The standalone package does not depend on the root library, and no root
library implementation changed. The root gate nevertheless ran unfiltered
as required. Before its final pass, the sparse checkout lacked compile-time
and runtime fixtures; a later attempt failed while writing to the full SSD.
Those failed logs remain beside the passing receipt. No test was skipped,
disabled or quarantined to obtain the pass.

No Python research code, harness or test suite changed. The relevant native
algorithm and example suites ran; unrelated historical Python autolab suites
are outside this change.

Actual local compiler: Homebrew rustc 1.93.1, LLVM 21.1.8, as retained in each
Cargo compiler fingerprint. Cargo was invoked through rustup 1.98 but found
the Homebrew compiler on PATH. Formatting/Clippy used their rustup drivers.
The CPU is Apple M4 Pro, 14 logical cores, 48 GiB, macOS arm64. Raw native
timing is excluded on this shared host. All measured field-call models,
binary hashes and constructor inputs are in the JSONL headers and rows.

CI checks for the eventual PR head and native Linux isolation execution are
separate from these local passes; they must be verified before merging.

## Post-merge validation, retained on 2026-10-15

At source `0191f5b924765a0f6d7882f13787b7391413c685`, the full standalone
release library/binary/integration/example suite passed again. The root suite
now contains twelve additional tests inherited from main. Its final unfiltered
run, with `RAYON_NUM_THREADS=1` and `--test-threads=1`, passed **4,095 tests,
zero failures, 117 existing ignored**. No test or assertion was changed.

The first post-merge full root attempt had 4,094 passes and one failure: the
unchanged Nagao six-sum test exhausted its existing solver deadline and
reported `timed_out: true`. That exact test passed alone with one Rayon worker.
A concurrent full retry was interrupted while the Nagao and registry fixtures
were still running; it has no complete-suite verdict. The successful serial
run does not replace either earlier receipt. These validation intervals are
diagnostic test-run data, excluded from algorithm speedup claims.

The fresh logs and command provenance are in `publication-20261015/`.
Earlier frozen archives and their inventories remain unchanged. The current
coordinator no longer contains old task T-2844 (404); resumed publication is
claimed as T-15 through a registered shell that maintains its lease.

## Final current-root and fourth-candidate validation

The unfiltered root release gate after merging main at `69c84f13a` passed
**4,104 tests, zero failed, 117 existing ignored**, test interval 3,454.89 s.
The registered shell returned exit zero. It used `RAYON_NUM_THREADS=1`,
`cargo test --release --lib --locked --offline -- --test-threads=1`, and the
build-SSD target/TMPDIR. Later commits change no root library bytes; a Git
comparison verifies the same `src`, root Cargo.toml and lockfile. Raw output
is sealed in `FINAL_VALIDATION.tar.gz`; its inventory is in
`validation-final-20261015/FILES.json`. Earlier failed receipts remain intact.

At fourth-candidate source `9a55ce6d22b9473c1d2afbf398bfdf51ae22b544`, the full
standalone release library/binary/integration/example suite, strict release
all-target Clippy, formatting and workflow lint pass. The five native CLI
tests pass. Both complete 96-case cumulative/incremental comparisons pass;
all four field counters match the previous rung in all 96 cases. The frozen
archive `RUNG4_RUNS.tar.gz` passes all thirteen byte/SHA256 records after fresh
extraction. The tested old checkout and continuation crate compare byte for
byte. These controls establish correctness, not a timing gain for rung four.

The first CPU wave independently audits all four full output matrices. AMD
has twenty valid native isolated timing commands and 640 recomputed medians;
the two macOS profiles have twenty valid count-only runs each. The ARM full
campaign and a deliberately wrong-source audit are rejected. The final CPU
packet passes ten byte/SHA256 checks, all four original ZIP CRC checks, and
byte-for-byte extraction comparison for the audit archive. The CPU and
fourth-candidate PDFs were visually inspected; PDF/SVG/150-dpi PNG browsing
copies match their study originals byte for byte. Native helper sources also
compile and pass formatting. This validation does not make unrun CPU waves,
other algorithm families or release publication complete.

The native legacy census compiles, passes formatting, and rejects empty
outputs and records missing their boolean verification flag. AMD, ARM and M4
each retain 395 records, including the same false path-search row. The local
M4 shell exits zero; a fresh archive extraction passes all eleven byte/SHA256
records, including executable, source, compiler, startup sample and resource
receipts. Its one-page report was visually reviewed and browsing copies match.
This is an exploratory coverage check, not a qualified timing comparison.
