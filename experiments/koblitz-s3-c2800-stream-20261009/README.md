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
raw report before promotion. It compares the algorithmic report fields with
the earlier C2800 T025 control, excluding timing, RSS, and the declared
capacity and scatter-policy fields. Exact semantic equality, full rank,
scalar verification, zero fallback, the expected table geometry, and peak
RSS within 32 GiB are required before assigning a candidate ID and running
fresh T157–T164 workloads.

## T025 feasibility result

The frozen control completed in 368.604 seconds of whole-process wall time.
It constructed 415,520,000 regular states, 1,073,741,824 root-table slots,
and 296,800 usable factor-base points. Preparation reached rank 2,800. The
target required 498 S3 calls, used zero fallback attempts, and recovered the
fixture scalar with independent group replay. The [semantic comparison](local/semantic_comparison.json)
reports equality for every compared algorithmic field.

[`/usr/bin/time -l` receipt](local/feasibility_receipt.json) recorded peak RSS
of **37,855,100,928 bytes**, exceeding the frozen 34,359,738,368-byte cap by
3,495,362,560 bytes. The compressed [raw report](local/feasibility_report.json.gz)
is preserved. The resource gate failed, so `candidate_id` remains `null` and
the fresh-target panel was not launched. The next variant fills the final
state allocation directly to remove the worker-chunk/merge overlap.
