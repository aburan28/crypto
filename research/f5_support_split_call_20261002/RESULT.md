# Support-separated matrix-F5 complete-call candidate

## Decision

The exact support-split certificate is implemented and passed the frozen
four-seed local correctness and counted-work screen. The required further
2× **isolated Linux complete-call result is still unmeasured**. These Apple
ARM64 timings are nonpromoting; each seed's exact five-pair bootstrap lower
bound is below 2×. Two smaller-case timing controls also fall below their
respective A/A minima. The physical Linux workflow is prepared in
`.github/workflows/f5-support-split-call.yml` and has not run.

This is a matrix-F5 solver-stage experiment. It is not an index-calculus
one-target online result or an IC versus rho comparison.

## Change and invariant

For a 24-variable degree-4 matrix with at least 4,096 rows and 8,192
columns, form 5 partitions degree-4 rows by whether they touch a degree-4
column supported entirely on variables 0–19. It proves full rank for the
inner and outer diagonal projections, then for each lower-degree diagonal
projection. The outer degree-4 rows are zero in the inner columns, and
lower-degree rows are zero in all higher-degree columns. Full row rank in
each diagonal block therefore proves the original full matrix has full row
rank. Only then are its unchanged original rows returned. Any failed
check falls back to exact full-matrix echelon reduction. Default reduced
output and existing selective echelon behavior are unchanged.

The primary returns a different row representation from selective
echelon. Every measured call matched exact rank, canonical row-space
fingerprint, F5 criterion work, built/pruned row counts, and column count.
All six smaller cases matched every nontiming output and route field.

## Matched local screen

One release binary ran both forms in separate processes, one Rayon
thread, with identical frozen input per pair. Each seed had two warmups,
five reference/reference A/A pairs, and five alternating
reference/candidate pairs. All 22 processes per seed completed. Ratios
below are selective-echelon / support-split; counted-work ratios are
support-split / selective-echelon. Binary and source hashes, every raw
JSON case, stderr, status, and phase measurement are in the receipts.

| Seed XOR | Complete-call median ratio | Exact bootstrap 95% lower | A/A ratio range | Reduction word operations: reference → candidate | Candidate work / reference |
| --- | ---: | ---: | ---: | ---: | ---: |
| `00000000` | 2.825× | 0.966× | 0.804–1.290× | 100,213,183 → 33,812,101 | 33.74% |
| `badc0de1` | 3.511× | 1.676× | 0.539–1.543× | 100,342,346 → 33,875,172 | 33.76% |
| `5eed2026` | 2.113× | 1.686× | 0.392–1.947× | 100,270,257 → 33,814,010 | 33.72% |
| `f5c02a28` | 3.868× | 1.360× | 0.232–1.541× | 100,365,756 → 33,809,701 | 33.69% |

The local counted-work promotion threshold was at most 55% on all four
seeds, and it passed. The local complete-call medians exceed 2×, but the
A/A ranges are wide. On smaller cases, frozen `n20_m20_d4` measured a
0.680× reference/candidate ratio against a 0.742× A/A minimum;
holdout `5eed2026` `n20_m20_d3` measured 0.426× against 0.799×.
Neither case used the support route or changed its nontiming output.
These control failures remain in the receipts and require the isolated
Linux gate before any speed conclusion.

## Validation and artifacts

- Release `cryptanalysis::matrix_f5_f2::tests`: 14 passed, including
  independent mixed-degree projection, rank-deficient fallback, zero-row
  fallback, and a cross-word projection.
- Release `cryptanalysis::gf2_elim::tests`: 9 passed.
- `bash -n research/f5_support_split_call_20261002/run_linux.sh` and
  `git diff --check` passed. The full repository `cargo fmt --check`
  reports extensive formatting differences in unrelated existing files;
  the three new Rust examples were formatted directly.
- The native Linux analyzer was smoke-tested against all four local
  receipts and synthetic clean-isolation records. It rejected them with
  the expected Apple-host, bootstrap, and smaller-case reasons; this is
  a parser and negative-gate check, not a Linux measurement.
- The Linux-only CPU chooser passed `cargo check --offline --target
  x86_64-unknown-linux-gnu --example f5_support_split_cpus` with the
  installed rustup stable compiler. This is a compilation check, not a
  physical x86-64 execution.
- The exact-source local receipts are `local_frozen.json.gz`,
  `local_holdout_a.json.gz`, `local_holdout_b.json.gz`, and
  `local_holdout_c.json.gz`. The first preformat screen is preserved in the
  matching `local_*_preformat.json.gz` files. `RECEIPTS.sha256` hashes the
  compressed archives; `RAW_SHA256.txt` hashes the original JSON bytes.
  All eight archives were extracted and matched their original hashes.
  Each run's benchmark source and
  binary hashes remain embedded in the receipt. The native benchmark
  runner's exact-source SHA-256 is
  `f94a0bbbd8867514bcbfa6513f6fc0ecc81a96df7755fa20b5bd460f488243b8`;
  the matrix-F5 source SHA-256 is
  `b7bca61827c6f3b73b14bd9a1522ba271246a03e0ca1c4bc890ea666a4a9934c`.

The Linux gate preserves the first clean isolated attempt on each seed,
with up to three attempts and all failed/contended attempts retained. A
one-thread pass requires the median **and** exact bootstrap lower bound
to exceed 2.00× on every seed. Smaller cases must avoid regression beyond
their A/A lower bound; a separate two-thread job enforces the same
nonregression control. A timeout, contention, failed output check, or
unavailable host is not a pass.

For the local replay, build the release `f4_f2_bench` and
`f5_support_split_pair` examples, then run the latter under
`python3 tools/isolated_bench.py busy --` with arguments
`target/release/examples/f4_f2_bench`, a fresh receipt path, one of
`frozen`, `holdout_a`, `holdout_b`, `holdout_c`, and `1`. The Python
isolation tool is the repository-mandated host reservation wrapper;
the benchmark driver, output checks, and gate are native Rust.
