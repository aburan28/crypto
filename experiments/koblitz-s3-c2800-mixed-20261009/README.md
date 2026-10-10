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

After that gate, [`prepare_freeze.py`](prepare_freeze.py) binds the C1200 and
C2800 candidate manifests, the same eight public points, and the paired
[`run_paired.py`](run_paired.py) runner. Its screen orders each point's two
arms AB or BA, caps each process at 600 seconds and 32 GiB RSS, and stops
after the first invalid or over-budget row. A complete verified screen needs
at least a twofold geometric mean reduction in S3 calls before further
measurement; online wall time and its target-level uncertainty remain
separate. [`preserve_screen.py`](preserve_screen.py) and
[`analyze.py`](analyze.py) retain all launched rows and leave aggregates
unknown when the panel is incomplete.
