#set page(paper: "a4", margin: 17mm, numbering: "1")
#set text(font: "DejaVu Sans", size: 10pt)
#set par(justify: false, leading: 0.65em)
#set heading(numbering: none)
= Boolean F4/F5 construction allocation optimization

Repository: `aburan28/crypto`. Study: `f4_f5_row_build_20261011`.
Engineering result for the existing Rust Boolean matrix implementation.

Complete-call cumulative requested allocation bytes fell from 11,667,696 to
6,712,608 for F4 (42.47\%) and from 12,400,920 to 6,762,696 for F5 (45.47\%).
F5 criterion allocation requests fell from 8,748 to 5,676 (35.12\%).
These are recorded values across 24 calls per arm, eight systems, three repetitions.
Percentages are native derived values: 100 x (1 - candidate/baseline).
Bytes are cumulative requests, not peak memory. Runtime remains unmeasured.

#figure(
  image("reports/run-20261011T1715Z/allocation-audit.svg", width: 100%),
  caption: [Recorded candidate audit, `run-20261011T1715Z`. No wall-time values.]
)

#table(
  columns: (1fr, 1fr, 1fr, 1.1fr, 1.1fr), inset: 5pt,
  align: (left, right, right, right, right),
  table.header([Arm], [Baseline requests], [Candidate requests], [Baseline bytes], [Candidate bytes]),
  [F4], [48,414], [48,294], [11,667,696], [6,712,608],
  [F5], [42,801], [39,507], [12,400,920], [6,762,696],
  [Criterion], [8,748], [5,676], [1,329,024], [660,480],
)

All matched output digests, ranks, pruning reports and XOR counters agree.
No workload regressed on either allocation measure. The formatted-driver A/A
control reproduced the first baseline exactly; five alternating baseline/candidate
pairs reproduced the primary JSONL records byte for byte.

== Implementation and verification

The shared Macaulay builder compacts product parity in place and collects
unique columns before the original monomial-order sort, using the existing
high-bit-safe mask hasher. The F5 criterion XOR-packs products directly into
reusable rows and retains buffers only for new pivots. Its fallback preserves
cancellation of missing columns. A monotone word cursor avoids rescanning the
zero prefix during reduction.

Five added tests cover independent Boolean multiplication, zero and constant
polynomials, bits through 63, cancellation of absent columns, and reduction
across prefix and temporary pivots in multiple machine words. Existing randomized
F4/F5 row-space and certificate tests remain active.

Validation: 17 F5, 43 Macaulay/Groebner and 13 F4 tests passed. The example test
target passed (zero unit tests); the actual audit driver executed every recorded
case. Across the full run and independent remainder, 3,838 distinct library tests
passed and 110 were ignored. One unchanged isogeny test did not finish before
the 20-minute budget; the original process was stopped with SIGTERM at 1,257
seconds (Cargo exit 101). Full-suite verification remains incomplete.

#pagebreak()
== Frozen conditions and reproducibility

The workloads are planted quadratic Boolean systems with (n,m,degree) =
(8,8,3), (8,8,4), (10,10,4), (10,10,5), each with seeds
`0x2545f4914f6cdd1d` and `0x9e3779b97f4a7c15`. The second seed is a holdout
relative to the existing F4 benchmark. Warmup precedes three repetitions.
Fixture generation, worker startup, hashing and verification are excluded;
complete packing, elimination and output construction are included.
Criterion construction is a separate arm.

Host: Linux x86-64, AMD EPYC 9V74, 5 logical CPUs, 18,882,699,264 bytes RAM,
Rust 1.99.0 (b940084d7), Cargo 1.99.0, default release build, one Rayon worker.
No solver-tuning environment overrides or RUSTFLAGS were set. No exclusive CPU
reservation was earned. No wall-time, peak-memory, other-field, larger-system,
full-DLP, matched-rho or m=83/m=131 scaling conclusion is established.

Original main `c95a14de9bfa2711643a9318f3728d5d1109367b` could not compile.
Separate prerequisite commits removed duplicated constructor, helper, assertion,
return and statistics fragments while preserving intact existing implementations.
The test build also required restoring an undefined parser fixture and a
test-only historical sort helper; the original assertions remain. Failed builds
(exit 101) and successful builds are retained; failures produced no accepted
audit values.

Baseline production: `e957e48c2b8cf1f84a7fc44a1d673103964a789d`.
Test-only repair and formatted driver: `ee344374`.
Optimized kernels: `2e3dffda681a4771be511653180da8354de3bce3`.
All source and executable hashes are in `records/*sha256.txt`;
the saved baseline executables and rlib remain under `/workspace/scratch/`.

Commands and exact baseline-library reconstruction are in `REPORT.md`.
Native execution:
```sh
cargo build --locked --release --example f4_f5_row_build_bench -j 2
target/release/examples/f4_f5_row_build_bench 3 > candidate.jsonl
target/release/examples/f4_f5_row_build_bench compare baseline.jsonl candidate.jsonl
RAYON_NUM_THREADS=1 cargo test --locked --release --lib -- --test-threads=1
cargo test --locked --release --example f4_f5_row_build_bench
```

== Evidence and delivery status

Implementation and matched allocation checks are verified; full library tests
are partial: 3,838 passed, 110 ignored, one unchanged test unfinished. Its name
and the full/remainder/termination logs are retained in `records/`; its sources
do not call the optimized kernels. Budget enforcement was delayed because
the checker initially could not identify the executable through `/proc`.
Work is committed locally on `codex/optimize-f4-f5-row-build`.
Initial GitHub authentication failed with HTTP 401. Subsequent browser
authorization restored authentication as `aburan28`; repository push access
and the branch upload dry run were verified. The original failure records and
manifest workflow state are retained as historical evidence.
Conductor is not installed and this checkout has no workflow file.

Native evidence: `records/comparison.json`, both primary JSONL files,
`records/per-system-counts.csv`, `records/paired/`, host/compiler/hash files,
and focused/full test logs. Protocol: `PROTOCOL.md`; provenance: `manifest.json`.
Run reports and inspected PDFs:
`reports/run-20261011T1655Z/`, `reports/run-20261011T1700Z/`,
`reports/run-20261011T1715Z/`, `reports/run-20261011T1716Z/`.
SVG data and source are retained in these directories and `render_reports.py`.
Python only renders recorded native results.

Canonical graphs checked and unchanged:
`docs/index-calculus-scoreboard.html`, `docs/ic/progress-timeline.json`,
`docs/performance-gains.html`, `docs/performance-gains/`,
`docs/curves/`, `figures/`. No time, IC/rho, curve or scaling fact changes.
