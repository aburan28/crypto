# N53 C2200 root-index feasibility

The C2200 candidate retains the C1200 point-seeded target-state prefix and
queries a larger full root table. Its frozen feasibility control checks
whether 256,520,000 indexed states, 513,040,000 root insertions, full rank
2200, and verified one-target scalar recovery fit within a 32 GiB peak RSS
envelope. The root-table gate is at most 536,870,912 packed slots, the same
allocation step observed for C1800. The pre-run actual factor-base count and
candidate ID remain null until the exact base has been enumerated and
verified.

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
