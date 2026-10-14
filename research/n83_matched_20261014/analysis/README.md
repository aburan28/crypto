# Native matched-curve accounting

This prospective Rust crate summarizes the frozen 32-target / three-cohort
degree-83 experiment. It performs all statistical calculations; the PDF and
plots read its recorded JSON values. Every ratio is **IC/rho**, so a value
below one favors IC.

```sh
cargo test --release --manifest-path research/n83_matched_20261014/analysis/Cargo.toml
cargo run --release --manifest-path research/n83_matched_20261014/analysis/Cargo.toml -- INPUT_JSON NEW_OUTPUT_JSON
```

The input has schema `n83_matched_statistics_input_v1`, `metadata`,
`process_rows`, `target_rows`, and `independent_verification`.
The verification object must say `status: PASS`, `checked_target_rows: 218`,
`invalid_rows: 0`, `checked_aa_semantic_pairs: 5`, `aa_semantic_gate: true`, and
contain the 64-digit `campaign_manifest_sha256` and
`protocol_sha256`. Metadata must have `independently_verified: true`, matching
protocol/manifest hashes, and `independent_verification_receipt_sha256`.
This is a pipeline binding gate, not proof against fabricated arbitrary JSON.
The five AA semantic pairs are checked by independent native C++: identical
public point, seed, policy, outcome/counters and per-seed preparation charges.
Different work cannot qualify as a pure timing-noise observation.
`metadata.target_ids` contains the exact
ordered 32 frozen IDs; `metadata.cold_target_ids` contains the first eight.
The controller emits the following after the independent correctness gate:

- Process rows: `kind` (`cohort`, `cold`, `aa`), `arm` (`ic`, `rho`),
  `run_id`, `cohort_id` (1–3 for cohorts), `target_id` for singleton/AA,
  `target_count` (32 or 1), `status`, `all_verified`,
  `elapsed_process_wall_s`, `cpu_user_s`, `cpu_system_s`, `cpu_total_s`,
  `peak_rss_bytes`, and `producer_preparation_ms`. Rho cohort rows additionally
  have `producer_shared_preparation_ms` and `producer_per_target_preparation_ms`:
  the first is `setup_total - all_G_only_jumps`, and the second is
  `all_G_only_jumps / 32`, taken from the frozen producer summary. AA rows additionally
  have `pair_id` (1–5), `leg` (`A`, `B`), and matching `seed`.
- Exactly 192 cohort target rows: `target_id`, `cohort_id` (1–3), `arm`,
  `status`, `verified`, `online_wall_s`. Cold/AA target rows remain in the
  raw controller evidence and are not duplicated in `target_rows`.
- Accepted completion statuses are `VERIFIED`, `VERIFIED_SOLVE`,
  `VERIFIED_SOLVED`, `COMPLETED_VERIFIED`, `COMPLETED`, or `SUCCESS`, accompanied by the required
  true verification flag. A true timeout/wrong-answer flag or a nonzero
  supplied exit code rejects the row. Verification flags must reflect the
  independent checker; the analysis does not replay scalars itself.

Output is create-new. Source identities, binary identities, input/output hashes,
and the exact command belong to the controller's provenance receipt.

The primary batch estimator is the ratio of arithmetic mean actual whole-process
costs across three paired 32-target processes. Its paired-cohort percentile
bootstrap has only three setup/timing clusters. The cold estimator uses eight
actual paired one-target fresh processes. Online times are secondary: stock IC
and rho have different input-validation/serialization clock boundaries.

All bootstrap calculations use 20,000 draws, SplitMix64 and unbiased rejection
sampling of indices. Seeds are 202610146799 for batch, 800 for cold, 801 for
target-block online, 802 for crossed target/cohort sensitivity, and 803 for the
break-even model. These are the respective final three digits; full seeds are
202610146799 through 202610146803. Target-block resampling retains all three
observations of each selected target. Crossed sensitivity independently
resamples paired targets and paired cohort indices. It does not create 96
independent setups or claim independent-host uncertainty.

Arithmetic means, midpoint medians, interpolated p90/p95, extrema and every
individual target remain visible. No fastest seed or favorable estimator is
selected after timing. Incomplete/duplicated/unverified required matrices
withhold their predeclared comparison. An overall winner requires all 208 main
answers. A locally valid eight-target cold comparison does not imply completion
of the 32-target cohort comparison. The independent 218-row audit gate remains
required for either primary. All observed consumed failure costs remain in
descriptive summaries.

The break-even curve is a linear **model**, using observed producer preparation
and mean online clocks. IC prep is invariant across the retained-base batch.
Rho invariant prep is the observed total minus all seed-specific jump-table
prep; the mean jump-table cost is added to each target's rho online cost, so
that extrapolation to N targets charges N fresh seeds rather than treating
all 32 jump tables as free fixed work. All preparation is still charged once
in the actually observed batch process cost.

The model does not include residual launch/output/audit costs and does not
establish execution above 32 targets. Crossed bootstrap samples paired cohort
preparation and target observations using the same sampled cohort indices.
Draws with no eventual IC advantage retain their nonfinite mass; finite-only
quantiles are explicitly conditional and are not substituted for an
unconditional interval. If IC has cheaper setup but slower per-target work,
its possible initial advantage and last favorable N are classified separately
from an eventual amortized crossover. Physical
host replication, common operation units and original retained-base construction
remain outside the evidence boundary.
