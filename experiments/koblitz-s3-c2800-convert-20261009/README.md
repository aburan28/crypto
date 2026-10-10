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

The T025 control passed those gates. The 415,520,000 ordered regular states
produced 415,077,662 distinct roots in 1,073,741,824 slots; prepared rank
reached 2,800, and 498 target S3 calls recovered scalar 429728658711 with
group replay. The [semantic comparison](local/semantic_comparison.json) is
equal to the full-key direct reference for every field outside the declared
layout and timing fields. Peak whole-process RSS was 18,463,309,824 bytes,
within the frozen 32 GiB envelope. The whole-process interval was 471.875 s;
the target online interval was 9.499 ms on this contended host and is a
diagnostic. The candidate manifest is
[`IC1N53Ckb1fb296800PDP4rootRCsampleLAgaussTDdirectISO0h1bc6d6637960`](identities/candidates/IC1N53Ckb1fb296800PDP4rootRCsampleLAgaussTDdirectISO0h1bc6d6637960.json).
The [receipt](local/feasibility_receipt.json) and
[compressed report](local/feasibility_report.json.gz) retain the raw result.
