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
