# Post-campaign native presentation supplement

This separately versioned Rust exporter adds descriptive fields omitted from
the frozen primary analysis: all 104 IC five-phase rows and signed residuals,
114 rho integer phase/residual rows, and all 32 process CPU, memory and timing
boundary records. The frozen workers, controller and analysis remain unchanged.
No primary estimator, bootstrap, policy or seed is recalculated or selected.

It runs only after the full 218-answer independent PASS receipt, five A/A
semantic pairs and the unchanged successful frozen analysis are available.
The input and frozen analysis rows must match. Observation/producer source
pointers and hashes accompany each result. Hashing uses native `sha256sum`.
Values used for this supplement come from original producer `targets.jsonl`,
`stdout.jsonl`, and `execution.json`. Original versus normalized copies retain
exact JSON structure, integer types, strings, identities and provenance. Only
explicit diagnostic binary64 clock/CPU fields permit one representable step
of parser/serializer drift; two steps or a changed integer/type is rejected.
The frozen primary analysis remains read-only and is not recomputed.

```sh
cargo test --release --offline --manifest-path research/n83_matched_20261014/postprocessing/Cargo.toml
cargo build --release --offline --manifest-path research/n83_matched_20261014/postprocessing/Cargo.toml
n83-matched-supplement INPUT_JSON FROZEN_STATISTICS_JSON CAMPAIGN_ROOT NEW_OUTPUT_JSON research/n83_matched_20261014/records/campaign-20261014T1940Z
```

Output is create-new. IC residual is exactly the signed floating calculation
`target_ms - target_phase_sum_ms`; it is not clipped, assigned to a sixth
phase, or added back into the primary process metric. The exported five-phase
sum is checked against its five components. The small comparison tolerance is
16 times binary64 epsilon times the largest magnitude (or one), for decimal
round-trips and addition only; raw residual values are preserved.

Rho phase identities use exact integer nanoseconds. CPU user plus system is
checked against total, and native KiB peak RSS must convert exactly to bytes.
Setup totals already include nested index/load/basis phases and are never
added to those phases twice. Process peak RSS is not per-target memory.
Medians/p90/p95 use the same conventional linear interpolation as the frozen
analysis. Descriptive tails have no new population or host uncertainty claim.

Rho work counters remain exact integers. The exporter retains each outcome's
walk steps, walks, fruitless cycles/escapes, capped walks, wasted merges and
table entries, then checks sums for integer overflow. It emits all-114 totals,
main-104 totals excluding the ten A/A rows, per-process totals (including
each shared cohort), and cohort/cold/A/A group totals. A zero process-timeout
count does not mean zero internal capped walks/restarts. Summed table entries
describe recorded work; they are not simultaneous peak memory.
