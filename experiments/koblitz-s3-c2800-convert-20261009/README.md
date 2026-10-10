# Full-key construction with exact seven-byte C2800 lookup

The [converted-root source](../../examples/s3_partial_prepared_fullrank_c2800convert_targetspan.rs)
builds each root shard with a full canonical key and its compact 31-bit
state/root witness. After all canonical roots have been inserted in first-
witness order, it converts each shard slot by slot to the seven-byte layout:
24 hash-tag bits and a 30-bit state index plus one root selector. Query-time
tag matches reconstruct and compare the complete key. The conversion retains
hash-table slot positions, duplicate semantics, and missing-key probe chains.
It removes random state reconstruction from the large index build. Twelve
release tests pass, including forced compact/full-key small-index equality.
Frozen source commit: `c28cd00a2`.

[`FEASIBILITY_FREEZE.json`](FEASIBILITY_FREEZE.json) binds source, binary,
runner, T025 point, and the previously verified direct-fill report. It fixes
rank 2,800, 1,073,741,824 table slots, 14 preparation workers, one target
worker, 32 GiB maximum whole-process RSS, and a 1,200-second cap. The
compressed raw report and `/usr/bin/time -l` receipt preserve scalar replay,
semantic equality, memory, instructions, cycles, and scheduling counters.
Online time is a diagnostic in this control.

Run only after the direct-compact control has finished:

```sh
python3 experiments/koblitz-s3-c2800-convert-20261009/run_feasibility.py \
  --binary target/release/examples/s3_partial_prepared_fullrank_c2800convert_targetspan \
  --receipt /Volumes/SSD990/llm/tmp/s3-c2800convert-feasibility-receipt.json \
  --raw /Volumes/SSD990/llm/tmp/s3-c2800convert-feasibility-report.json.gz
```

[`record_feasibility.py`](record_feasibility.py) archives the result and
compares all algorithmic fields with the direct-fill reference, except the
declared slot format and construction policy. It assigns a content-hashed
candidate ID and freezes unused T165–T172 workloads only after the rank,
scalar, table geometry, semantic, and memory gates pass.
