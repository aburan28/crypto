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

## Feasibility result and exact identity

The [receipt](local/feasibility_receipt.json) and [raw report](local/feasibility_report.json.gz) show verified scalar recovery with full rank 1,580 before target input, zero online fallback, 167,480 distinct subgroup-usable base points, and 132,309,200 indexed states. The actual base digest is `5be73166ac72caba26766e1f8a62f957fb54adb88aabac784f483c8b8610aa4f`; its last scanned `x` is 3,180. Every shard stayed within the declared capacity step: 268,435,456 total 13-byte root slots. The observed peak RSS was **20,308,869,120 bytes**, below the 32 GiB envelope. T025 required 1,048 target S3 root calls, versus 3,272 in the compact-state C1200 control; this known-point count is a feasibility diagnostic.

The exact candidate is [`IC1N53Ckb1fb167480PDP4rootRCsampleLAgaussTDdirectISO0h2090cce30506`](identities/candidates/IC1N53Ckb1fb167480PDP4rootRCsampleLAgaussTDdirectISO0h2090cce30506.json). The 24 previously unused N53 targets T101–T124 retain their [one-target workload manifests](identities/workloads/); the first eight are reserved for a paired screen, with the remainder held for later confirmation.
