# Fresh N53 C1200/C2200 prefix call-count screen

[`SCREEN_FREEZE.json`](SCREEN_FREEZE.json) fixes eight previously unused
one-target workloads T125–T132 on the same N53 public curve and exact target
point in each pair. It compares C1200 with the C2200 root table and C1200
target-state prefix. Both arms use the same committed source and release
binary; candidate manifests fix the distinct factor bases, ranks, search
orders, and source hashes. Runs alternate AB/BA by point, use 14 preparation
workers and one target worker, and retain all target-dependent online phases
and scalar replay.

This is a **one-repetition screening gate** for deterministic target S3-call
counts. Advancement to a replicated fresh-target wall panel requires all 16
scalars verified and a target-level geometric mean C1200/C2200 call ratio of
at least 3. The threshold, targets, 600-second per-run cap, and 32 GiB
memory envelope were fixed before execution. A wall ratio from this local
screen is exploratory and does not establish the requested 2× online speedup.

Run the frozen screen with:

```sh
python3 experiments/koblitz-s3-c2200-prefix-screen-20261009/run_paired.py \
  --binary target/release/examples/s3_partial_prepared_fullrank_c2200prefix_targetspan \
  --output /private/tmp/s3-c2200-prefix-screen-20261009
```

The runner validates candidate/workload hashes, source and binary hashes,
public points, full-rank preparation, root-table sizes, search policies,
exclusive online phase sum, and independently replayed scalars. Its
append-only ledger retains timeouts, failures, and complete raw reports.
