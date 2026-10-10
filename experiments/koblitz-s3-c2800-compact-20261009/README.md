# Exact seven-byte root slots for N53 C2800

The [compact-root source](../../examples/s3_partial_prepared_fullrank_c2800compact_targetspan.rs)
keeps the direct-filled 415,520,000-state array and replaces the C2800
13-byte root slots with seven-byte slots. Each slot stores a 24-bit hash tag,
a 30-bit state index, and one root selector. On a matching tag, the lookup
reconstructs the complete 53-bit canonical key from the indexed polynomial
root before accepting the witness. Insertions also check complete keys, so
hash-tag collisions preserve exact first-witness semantics. The frozen
1,073,741,824-slot layout has 7,516,192,768 bytes of slot payload, a
6,442,450,944-byte reduction from the 13-byte layout. Other base sizes keep
the established full-key root table. Twelve release tests pass, including
compact growth, missing-key rejection, and small-index equivalence to the
full-key table. Source commit: `257c2818a`.

[`FEASIBILITY_FREEZE.json`](FEASIBILITY_FREEZE.json) binds source, binary,
runner, T025 public point, and the direct-fill reference report before the
control. It fixes 14 preparation workers, one target worker, rank 2,800,
1,073,741,824 expected slots, a 32 GiB whole-process RSS envelope, and a
1,200-second cap. `/usr/bin/time -l` captures process RSS, instructions,
cycles, and scheduling counters. Rank, scalar replay, table geometry, and
exact report semantics must pass before a candidate ID or fresh comparison
is assigned.

Run the committed control after the ongoing direct-fill paired screen has
finished, so the two large index builds do not contend:

```sh
python3 experiments/koblitz-s3-c2800-compact-20261009/run_feasibility.py \
  --binary target/release/examples/s3_partial_prepared_fullrank_c2800compact_targetspan \
  --receipt /Volumes/SSD990/llm/tmp/s3-c2800compact-feasibility-receipt.json \
  --raw /Volumes/SSD990/llm/tmp/s3-c2800compact-feasibility-report.json.gz
```

[`record_feasibility.py`](record_feasibility.py) archives every completed
outcome and compares the algorithmic fields with the direct-fill T025 report.
It excludes timing, RSS, and the declared slot representation only. A passing
control generates a content-hashed candidate manifest and previously unused
T165–T172 workload manifests for the next paired screen.
