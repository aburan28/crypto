# N53 C1800 compact-state feasibility control

This frozen control tests whether a full-rank 1,800-column N53 factor base fits a 32 GiB observed peak-memory envelope using the 40-byte retained-state implementation. A 1,800-column base would index 171,720,000 regular states and 343,440,000 input roots if every pair lifts. Its table crosses the current capacity step to at most 536,870,912 13-byte slots. The actual base size, digest, shard counts, rank, and memory are measured in the control.

[`FEASIBILITY_FREEZE.json`](FEASIBILITY_FREEZE.json) fixes the committed source and release binary hashes, public T025 point, relation seed `20260928`, 14 preparation workers, one target worker, 360-second timeout, and 32 GiB observed RSS gate. This is a design-feasibility run: candidate ID and measured IC costs remain null until the exact base and execution gates are established. The raw result is preserved on success, failure, timeout, or OOM.

Run after the active C1200/C1580 paired screen completes, to avoid resource contention:

```sh
python3 experiments/koblitz-s3-c1800-20261009/run_feasibility.py \
  --binary target/release/examples/s3_partial_prepared_fullrank_compactstate_targetspan \
  --receipt experiments/koblitz-s3-c1800-20261009/local/feasibility_receipt.json \
  --raw experiments/koblitz-s3-c1800-20261009/local/feasibility_report.json.gz
```

If full rank and resource gates pass, assign the exact `IC1` identity and freeze fresh one-target point comparisons before interpreting online root-call or wall-time ratios. The local CPU host does not have a host-wide isolation receipt for a controlled wall-time claim.
