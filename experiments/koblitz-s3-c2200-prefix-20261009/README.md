# N53 C2200 root-index feasibility

The frozen C2200 control reached rank 2200 and recovered the T025 scalar
with no target fallback. Its exact factor base has **233,200** usable points;
the index contains **256,520,000** states and **256,347,555** distinct root
keys in **536,870,912** packed slots. The observed child peak RSS was
**17,543,626,752 bytes**, below the predeclared 32 GiB cap. This is one
macOS process observation and will be checked again in subsequent runs.
The exact candidate is
[`IC1N53Ckb1fb233200PDP4rootRCsampleLAgaussTDdirectISO0hd2accd16723d`](identities/candidates/IC1N53Ckb1fb233200PDP4rootRCsampleLAgaussTDdirectISO0hd2accd16723d.json).
Eight new T125–T132 [one-target workloads](identities/workloads/) are
reserved for the next paired screen.

The C2200 candidate retains the C1200 point-seeded target-state prefix and
queries a larger full root table. Its frozen feasibility control checks
whether 256,520,000 indexed states, 513,040,000 root insertions, full rank
2200, and verified one-target scalar recovery fit within a 32 GiB peak RSS
envelope. The root-table gate is at most 536,870,912 packed slots, the same
allocation step observed for C1800. The frozen pre-run record kept the
actual factor-base count and candidate ID null until enumeration and
verification completed.

[`FEASIBILITY_FREEZE.json`](FEASIBILITY_FREEZE.json) fixes the previously used
T025 correctness target, source commit, release binary, target point, relation
seed, 14 preparation workers, one target worker, 600-second cap, and memory
and slot gates. This is a correctness and memory control. Fresh timing points
remain unused until a complete candidate manifest and a separate screen
protocol have been committed.

Run the frozen control with:

```sh
python3 experiments/koblitz-s3-c2200-prefix-20261009/run_feasibility.py \
  --binary target/release/examples/s3_partial_prepared_fullrank_c2200prefix_targetspan \
  --receipt /private/tmp/s3-c2200prefix-feasibility-receipt.json \
  --raw /private/tmp/s3-c2200prefix-feasibility-report.json.gz
```

The runner retains unsuccessful statuses and checks the complete online
report, full-rank preparation, exact scalar, schedule policy, root-table
capacity, and observed child peak RSS.
[`local/feasibility_receipt.json`](local/feasibility_receipt.json) records the
observed gates, status, source and binary hashes, and child peak RSS;
[`local/feasibility_report.json.gz`](local/feasibility_report.json.gz)
retains the full native report and relation witnesses. The T025 target is a
feasibility control; its target S3-call count of 498 equals the C1800
control's count and does not estimate fresh-target yield.
