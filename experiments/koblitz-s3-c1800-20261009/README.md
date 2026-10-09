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

## Feasibility result and exact identity

The [receipt](local/feasibility_receipt.json) and [raw report](local/feasibility_report.json.gz) show full rank 1,800 before target input, zero online fallback, and independently checkable scalar replay. The base has **190,800** usable points, digest `0c3597d24ca943d7bb081ee30078596ffc433766356a8042ba8180c2ced78030`, and last scanned `x` 3,623. Its 171,720,000 states use a 6,868,800,000-byte retained-state vector; the root table has 171,640,231 distinct keys in 536,870,912 slots. Observed peak RSS was **22,582,263,808 bytes**, below the frozen 32 GiB envelope. T025 used 1,270 target S3 root calls, a known-point diagnostic.

The exact candidate is [`IC1N53Ckb1fb190800PDP4rootRCsampleLAgaussTDdirectISO0h904e60fe9a23`](identities/candidates/IC1N53Ckb1fb190800PDP4rootRCsampleLAgaussTDdirectISO0h904e60fe9a23.json). Fresh N53 T109–T124 [one-target workloads](identities/workloads/) are reserved for a later paired panel. The C1580 fresh-target screen found strong point-dependent reversals, so target search order should be fixed and tested on this holdout before expanding it to C1800.
