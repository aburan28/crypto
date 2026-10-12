# Secondary descriptive export

This package is separate from the frozen runner, protocol and primary paired analysis. It prepares renderer data after timing has ended. It does not solve a target, create new timing, bootstrap, generate a new confidence interval or replace the primary geometric mean of paired online ratios.

**Do not build or run it during the timed panel.** Source preparation alone leaves the active benchmark unchanged.

After timing, build from the repository root with a local Rust toolchain:

```sh
CARGO_TARGET_DIR=target/n83-descriptive-summary cargo build --release --locked --manifest-path research/n83_evaluation_20261012/descriptive_summary/Cargo.toml
CARGO_TARGET_DIR=target/n83-descriptive-summary cargo test --release --locked --manifest-path research/n83_evaluation_20261012/descriptive_summary/Cargo.toml
```

After the frozen primary analysis has produced `statistics.json`:

```sh
target/n83-descriptive-summary/release/n83-descriptive-summary --pairs /absolute/path/to/panel/pairs.jsonl --primary-statistics /absolute/path/to/panel/statistics.json --out /absolute/path/to/new/report_data_secondary.json
```

The output explicitly labels its min/max/arithmetic means **SECONDARY_DESCRIPTIVE_NOT_PREREGISTERED_PRIMARY_ESTIMATOR**. It embeds the unmodified primary statistics, protocol and pair rows and retains source, binary, input, receipt and raw-file SHA256 digests. It refuses to overwrite output, requires the five formal matched pairs, checks the supplied primary analysis's pairs hash, and validates native receipt schema/arm identity. No frozen input is written.

The separate `primary_reciprocal_presentation` field shows the frozen matched primary results in the opposite direction: `IC/strong = 1 / (strong/IC)`, for online and whole-process costs. Its interval is the original interval with bounds reversed and inverted, `[1 / original_upper, 1 / original_lower]`. This is an algebraic presentation of the original primary estimator, not another estimator or confidence interval. Only positive finite source values are transformed; missing/invalid values remain null, and no interval is supplied for fewer than two retained pairs.

Each of the five formal rows per method includes seed, status, online and preparation milliseconds, complete-process wall and user/system/total CPU seconds, peak RSS bytes, phase residual, resource envelope, and the original full receipt. A missing receipt produces null measurements and an explicit missing status. No zero imputation occurs.

Descriptions distinguish successful verified runs from all charged attempts. Online and preparation summaries use successful runs only. Wall, CPU and memory additionally have summaries over every observed formal attempt, including timeout and failure charges; those costs are never presented as successful-solve costs. Every summary carries its observed/missing counts and original per-pair values. Warmups and strong/strong A/A arms remain in the embedded source indexes but do not enter these per-method formal summaries.

Whole-process wall runs from before child spawn until the parent's `wait4(WNOHANG)` observes exit; the frozen runner sleeps 100 ms between polls. Exit observation can therefore add roughly 0–100 ms plus scheduler and polling-body delays, with no strict 100 ms bound under arbitrary scheduling. Retain this original wall measurement and use child CPU totals and producer process timestamps as corroborating measurements; do not subtract an inferred delay or revise the frozen primary analysis. The producer's `process_ns_before_receipt_io` has an earlier endpoint and excludes final receipt I/O.

The delivery unit tests cover reciprocal bounds, absence of a one-pair interval, invalid/missing ratios and intervals, and the separation of timeout charges from successful-solve costs without zero imputation. The exact source compiled for the original export is retained in `../evidence/source/descriptive_main_compiled.rs`; delivery-only unit tests do not alter that recorded export.
