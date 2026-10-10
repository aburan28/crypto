# Bounded-scatter N53 C2800 feasibility

The [bounded-scatter source](../../examples/s3_partial_prepared_fullrank_c2800stream_targetspan.rs)
inserts canonical roots into 64 open-addressed shards in consecutive chunks
of at most 4,194,304 index states. Each chunk retains at most 8,388,608
`(key, witness)` entries before insertion. Worker ranges and state blocks are
consumed in their original order per shard, preserving first-witness choices.
Shard tables start from a uniform capacity and grow exactly if an input is
more skewed than that initial estimate. The source is committed at
`2862f1561`; nine release tests pass, including shard growth and serial versus
parallel index equivalence.

[`FEASIBILITY_FREEZE.json`](FEASIBILITY_FREEZE.json) fixes the N53 curve,
T025 public point, factor-base construction, 2,800 prepared columns,
415,520,000 expected regular states, 1,073,741,824 expected table slots,
14 preparation workers, one target worker, a 900-second cap, and a 32 GiB
peak-RSS envelope. It fixes source, binary, runner, point, and prior C2800
reference-report hashes before execution. `/usr/bin/time -l` measures the
whole child process's maximum resident set size and operation counters.
The run is a correctness and resource control; online time is a diagnostic.

Run the committed protocol with:

```sh
python3 experiments/koblitz-s3-c2800-stream-20261009/run_feasibility.py \
  --binary target/release/examples/s3_partial_prepared_fullrank_c2800stream_targetspan \
  --receipt /Volumes/SSD990/llm/tmp/s3-c2800stream-feasibility-receipt.json \
  --raw /Volumes/SSD990/llm/tmp/s3-c2800stream-feasibility-report.json.gz
```

[`record_feasibility.py`](record_feasibility.py) preserves the receipt and
raw report before promotion. For a passing result, it compares every
non-timing report field with the earlier C2800 T025 control, except the
declared capacity and scatter-policy fields. Only exact semantic equality,
full rank, scalar verification, zero fallback, the expected table geometry,
and peak RSS within 32 GiB permit a candidate ID and fresh T157–T164
workload records.
