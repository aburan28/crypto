# N53 C1580 compact-state feasibility control

This control tests whether the 40-byte retained-state S3 implementation can prepare a full-rank 1,580-column N53 factor base within a 32 GiB observed peak-memory envelope. It uses the already checked T025 public point to establish the actual base size and digest, table capacity, relation rank, scalar recovery, and resource use. The [freeze](FEASIBILITY_FREEZE.json) fixes source and binary hashes, input point, relation seed, 14 preparation workers, one target worker, and a 300-second limit before execution.

The candidate ID and measured IC costs remain null until the exact base and full algorithm gates are established. This control is a design-feasibility measurement, not a paired one-target speed comparison. Its raw solver output is preserved even if rank, memory, or scalar verification fails.

If all 1,580² × 53 indexed states are regular, the state count is 132,309,200 and the two-root input count is 264,618,400. The present dense table has a 268,435,456-slot capacity step so long as every one of 64 shards receives at most 4,194,303 inputs. The actual shard counts and capacity are measured, not assumed.

After the frozen controls and the [compact-state N41/N53 T025 checks](../koblitz-s3-compact-state-20261009/README.md), run this from the repository root:

```sh
python3 experiments/koblitz-s3-c1580-20261009/run_feasibility.py \
  --binary target/release/examples/s3_partial_prepared_fullrank_compactstate_targetspan \
  --receipt experiments/koblitz-s3-c1580-20261009/local/feasibility_receipt.json \
  --raw experiments/koblitz-s3-c1580-20261009/local/feasibility_report.json.gz
```

If the control passes, fix its actual base, candidate ID, and fresh one-target workloads before a matched comparison with C1200. Report all target-dependent attempts, exclusive phase costs, scalar replay, resource use, and failures; the local CPU host still requires an isolation receipt for a controlled wall-time claim.
