# Direct-fill N53 C2800 feasibility

The [direct-fill source](../../examples/s3_partial_prepared_fullrank_c2800direct_targetspan.rs)
constructs the 415,520,000 N53 C2800 regular index states in one preallocated
40-byte-per-state array. Fourteen workers fill disjoint consecutive left-row
ranges. This preserves serial state order and first-witness root insertion
while avoiding the full worker-chunk-to-final-array merge. Other geometries
retain the general variable-yield state builder. The bounded root-scatter
path remains the same as the preceding variant. Ten release tests pass,
including direct-fill order across worker counts and the small serial versus
parallel root-index control. Frozen source commit: `54737edf4`.

[`FEASIBILITY_FREEZE.json`](FEASIBILITY_FREEZE.json) fixes the N53 curve,
T025 public point, factor-base construction, 2,800 prepared columns,
415,520,000 expected regular states, 1,073,741,824 expected root-table slots,
14 preparation workers, one target worker, a 900-second cap, and a 32 GiB
peak-RSS envelope. It binds source, binary, runner, point, and C2800
reference-report hashes before execution. `/usr/bin/time -l` records the
whole process's peak RSS and operation counters. This control measures rank,
scalar correctness, report semantics, and the resource gate. Its target time
is a diagnostic until fresh paired targets and host isolation are measured.

Run the committed protocol with:

```sh
python3 experiments/koblitz-s3-c2800-direct-20261009/run_feasibility.py \
  --binary target/release/examples/s3_partial_prepared_fullrank_c2800direct_targetspan \
  --receipt /Volumes/SSD990/llm/tmp/s3-c2800direct-feasibility-receipt.json \
  --raw /Volumes/SSD990/llm/tmp/s3-c2800direct-feasibility-report.json.gz
```

Preserve every outcome. Advance to a fresh paired target panel only if full
rank, scalar replay, exact report semantics, 32 GiB peak RSS, and the frozen
root-table geometry all pass. An exceeded resource cap is a recorded result.
