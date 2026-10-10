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
