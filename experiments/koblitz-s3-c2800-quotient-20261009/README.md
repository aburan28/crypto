# Exact eight-byte quotient root index for N53 C2800

The [source](../../examples/s3_partial_prepared_fullrank_c2800quotient_targetspan.rs)
retains the 415,520,000 ordered regular states and inserts each canonical
root into a 64-shard table of eight-byte slots. Multiplication by an odd
constant modulo `2^53` permutes the complete canonical key. A slot stores
its 23-bit key remainder, six-bit initial probe displacement, 29-bit state
ordinal, and six-bit Frobenius shift. The shard and original home slot supply
the other key bits. Shard growth recalculates the split exactly. Query-time
lookup verifies the complete permuted key through its split fields and
reconstructs the `left/right/relative/shift` witness arithmetically, without
reading the large state array. Duplicate roots preserve the first witness.
Thirteen release tests pass, including exhaustive 13-bit key probing and
forced shard growth. Frozen source commit: `a845983f1`.

[`FEASIBILITY_FREEZE.json`](FEASIBILITY_FREEZE.json) binds the source,
release binary, runner, T025 point, and full-key direct reference. It fixes
rank 2,800, 1,073,741,824 table slots, 14 preparation workers, one target
worker, 32 GiB maximum whole-process RSS, and a 1,200-second cap. The
[`run_feasibility.py`](run_feasibility.py) receipt captures the raw report,
scalar replay, memory, instructions, cycles, and scheduling counters. The
[`record_feasibility.py`](record_feasibility.py) comparison excludes only
timings, peak RSS, the declared root layout, and per-shard counts changed by
the new hash partition. The single-target online interval is a diagnostic
until an auditable host isolation record is available.

Run this control after any other large C2800 index job has finished:

```sh
python3 experiments/koblitz-s3-c2800-quotient-20261009/run_feasibility.py \
  --binary target/release/examples/s3_partial_prepared_fullrank_c2800quotient_targetspan \
  --receipt /Volumes/SSD990/llm/tmp/s3-c2800quotient-feasibility-receipt.json \
  --raw /Volumes/SSD990/llm/tmp/s3-c2800quotient-feasibility-report.json.gz
```

If the rank, scalar, table geometry, semantic, and memory gates pass, the
recorder assigns the content-hashed candidate ID and freezes fresh T173–T180
workloads. Retain any failed or capped run as a receipt with its raw status.

The frozen T025 control completed in 117.312 s, reached rank 2,800, and
recovered scalar 429728658711 after 498 S3 calls with group replay. The
[semantic comparison](local/semantic_comparison.json) matches the full-key
reference after excluding timing and the declared table geometry. Peak RSS
was 27,307,655,168 bytes. The low-bit shard hash was skewed: input counts
ranged from zero to 29,937,272 per shard, and 29 of 64 shards doubled. The
table reached 1,560,281,088 slots, so the frozen slot-count gate failed.
The [receipt](local/feasibility_receipt.json) and
[compressed report](local/feasibility_report.json.gz) preserve this result;
`candidate_id` remains null and the T173–T180 panel was not started. The
target online interval was 0.701 ms on this contended host and remains a
diagnostic.
