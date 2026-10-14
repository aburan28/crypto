# Post-run descriptive supplement

This native Rust analysis adds medians and minima beside the frozen campaign's
paired means. It leaves the measured records, runner sources and primary
statistics untouched. It creates a new JSON output and refuses overwrites.

Every one of the thirty verified held-out rows must be present, with three
matching seeds per target and arm, the frozen five target/seed assignments,
one fixed configuration per arm, matching paired scalars, and positive finite
wall, CPU, RSS, preparation, online and walk-step values. Missing or failed
values reject the complete-panel supplement; none is replaced by zero.

The output contains fifteen-draw summaries for each arm, three-draw summaries
for each target and arm, and all fifteen matched current/selected process-wall
ratios. Each summary reports arithmetic mean, median, minimum and maximum.
These launch-ratio descriptions do not replace the frozen five-target primary
estimand and carry no population confidence interval.

```sh
cargo test --release --locked --manifest-path research/n83_speed_20261012/analysis/Cargo.toml
cargo build --release --locked --manifest-path research/n83_speed_20261012/analysis/Cargo.toml
n83-speed-descriptive research/n83_speed_20261012/records/panel/holdout/records.jsonl research/n83_speed_20261012/analysis/descriptive_summary.json
```

The output pins the input, compiled source, compiled Cargo manifest and lock,
executable and exact invocation. `test.log`, `build.log` and `run.log` retain
the executed checks and command disposition. Tests cover odd and even medians,
input ordering, a singleton, missing metrics, zero values and an incomplete
panel. This is analysis of an existing completed experiment; execution cost is
not presented as a new performance measurement.
