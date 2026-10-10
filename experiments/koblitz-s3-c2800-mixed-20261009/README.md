# Mixed eight-byte root index for N53 C2800

The [source](../../examples/s3_partial_prepared_fullrank_c2800mixed_targetspan.rs)
uses a reversible 53-bit permutation: xor with a 26-bit right shift,
multiplication by `0x9e3779b97f4a7c15` modulo `2^53`, and the same right-xor
again. The low six mixed bits select one of 64 shards. Each eight-byte slot
retains the rest of the exact mixed key through its home position, remainder,
and 12-bit probe displacement, plus a 29-bit ordered state index. On a
complete key match, the first matching root in that retained state supplies
the Frobenius shift. This keeps first-witness semantics and leaves up to
4,095 probe positions before a shard must grow.

The previous [quotient control](../koblitz-s3-c2800-quotient-20261009/README.md)
completed rank and scalar replay, but a low-bit-only hash skewed the input
among shards and exceeded its frozen table-size gate. Fifteen release tests
pass for this variant, including exhaustive 13-bit key lookup, first-witness
and growth checks, exact shift reconstruction, and a mixed-shard balance
check on N53 C244 roots. Frozen source commit: `fbcccb303`.

[`FEASIBILITY_FREEZE.json`](FEASIBILITY_FREEZE.json) binds the source,
release binary, runner, T025 point, and full-key direct reference. It fixes
rank 2,800, 1,073,741,824 table slots, 14 preparation workers, one target
worker, 32 GiB whole-process RSS, and a 1,200-second cap. The
[`run_feasibility.py`](run_feasibility.py) receipt keeps the raw report,
scalar replay, memory, instructions, cycles, and scheduling counters. The
[`record_feasibility.py`](record_feasibility.py) comparison excludes timings,
peak RSS, and declared table geometry while checking all other algorithmic
fields. The one-target online interval is a diagnostic without host-wide
isolation evidence.

```sh
python3 experiments/koblitz-s3-c2800-mixed-20261009/run_feasibility.py \
  --binary target/release/examples/s3_partial_prepared_fullrank_c2800mixed_targetspan \
  --receipt /Volumes/SSD990/llm/tmp/s3-c2800mixed-feasibility-receipt.json \
  --raw /Volumes/SSD990/llm/tmp/s3-c2800mixed-feasibility-report.json.gz
```

If rank, scalar, semantic, memory, and table-geometry gates pass, the
recorder assigns the content-hashed candidate ID and freezes fresh
T173–T180 one-target workloads. Preserve any failed or capped row.

The T025 control passed. Its 415,520,000 regular states produced 415,077,662
distinct canonical roots in exactly 1,073,741,824 slots, with every shard
remaining at 16,777,216 slots. Input counts ranged from 12,971,694 to
12,998,628 per shard. Prepared rank reached 2,800; 498 target S3 calls
recovered scalar 429728658711 with group replay. Peak whole-process RSS was
22,526,918,656 bytes, within the frozen 32 GiB envelope. The
[semantic comparison](local/semantic_comparison.json) matches the full-key
reference outside the declared layout and timing fields. The complete
process ran 216.328 s; the 5.110 ms target online interval was measured on a
contended host and remains a diagnostic. The [receipt](local/feasibility_receipt.json)
and [compressed report](local/feasibility_report.json.gz) preserve the raw
result. Candidate ID:
[`IC1N53Ckb1fb296800PDP4rootRCsampleLAgaussTDdirectISO0hbbd5b5a111fb`](identities/candidates/IC1N53Ckb1fb296800PDP4rootRCsampleLAgaussTDdirectISO0hbbd5b5a111fb.json).

After that gate, [`prepare_freeze.py`](prepare_freeze.py) binds the C1200 and
C2800 candidate manifests, the same eight public points, and the paired
[`run_paired.py`](run_paired.py) runner. Its screen orders each point's two
arms AB or BA, caps each process at 600 seconds, rejects a run whose measured
peak RSS exceeds 32 GiB, and stops
after the first invalid or over-budget row. A complete verified screen needs
at least a twofold geometric mean reduction in S3 calls before further
measurement; online wall time and its target-level uncertainty remain
separate. [`preserve_screen.py`](preserve_screen.py) and
[`analyze.py`](analyze.py) retain all launched rows and leave aggregates
unknown when the panel is incomplete.

The [archived local screen](analysis.json) completed all 16 runs on T173–T180
with native group and known-scalar checks, full-rank preparation, and peak
whole-process RSS below 32 GiB in every row. The geometric mean target S3
call ratio was **3.220×** (target bootstrap 95% interval **1.826–5.542×**):
five of eight points needed at least twofold fewer calls with C2800. The
corresponding verified one-target online wall ratio was **1.121×** (target
bootstrap interval **0.621–1.967×**), with one of eight points at or above
twofold. The largest measured process RSS was 25,814,581,248 bytes. The
baseline and expanded median reusable setup times were 30.976 and 179.971 s,
respectively; these are recorded separately from the online interval.
Median charged PDP time per target S3 call was 0.805 and 2.523 microseconds
for C1200 and C2800, respectively, on the same contended host.
This Mac screen lacks a host-level isolation receipt, so its timing ratio is
an exploratory diagnostic. The [ledger](local/raw_paired.jsonl.gz) and
[raw reports](local/raw_reports.tar.gz) retain every measured row. The
call-count gate passed. The next optimization step is to profile the wider
index's per-query work and memory access, which must fall for a twofold wall
gain on this workload.
