# Sixteen-point N83 panel bounds the K1104–K1120 cold difference

The sealed Tier 1 panel measured K1120/root11 **145.454 ms faster in mean
complete cold process time** than K1104/root11 across sixteen fresh public
points. The paired two-sided 95% interval for K1104 minus K1120 was
**[-92.776, 383.684] ms**. The frozen selector therefore records
`BOUNDED_UNKNOWN` and a null exact configuration. All 64 IC arms and 16
same-point rho solves passed independent arithmetic and resource replay.

| Complete pipeline | Usable base points | Folded columns | Mean cold process, s | Mean one-target online, ms | Mean paired rho online, ms |
| --- | ---: | ---: | ---: | ---: | ---: |
| K1104/root11 dynamic | 183,264 | 1,104 | 119.847 | 176.455 | 2,553.156 |
| K1120/root11 dynamic | 185,920 | 1,120 | 119.701 | 158.121 | 2,553.156 |

Each process built its base, root index and logs before solving one public
target. The online interval begins with target-dependent work after that
preparation. One independent rho process solved the same point for each of
the sixteen point workloads. The table reports descriptive means; the frozen
statistical test concerns paired **cold** process differences. The
[paired-difference graph](precision.svg), [run PDF](RUN.pdf), and
[results PDF](RESULT.pdf) retain these timing boundaries.

The [GitHub Actions run](https://github.com/aburan28/crypto/actions/runs/38970960376)
executed source head `d13e28e76d3709c19f15cf2c352d4ab912812335` on a
Tier 1 cloud instance. Its original provider artifact ID is `12083922301`;
the archived [provider receipt](provider_38970960376.json) records ZIP
SHA-256
`559578cd386693ce3681827ba3a9f192bd99be7ca72920fc35734cf515fd844d`.
The complete ZIP remains at
`/Volumes/SSD990/crypto/artifacts/n83_e2e_precision_20261017/raw.zip`.
The [cloud analysis](analysis_38970960376.json) and
[local replay](local_replay_38970960376.jsonl) match as parsed JSON with
empty stderr and exit status 0; the [replay receipt](replay_receipt_38970960376.json)
pins their hashes. The [IC1 record](ic1_rows_38970960376.json) binds all 645 ZIP result
members, sixteen workload IDs, 64 run rows, two exact candidate manifests,
the paired interval, phase timings, memory, host receipt, and all sixteen rho
receipts. `materialize.py` checks the provider digest and every extracted ZIP
member before writing the IC1 record. The independent verifier checks each
base, relation, scalar, source, schedule, topology and telemetry receipt.

The next separately sealed finalist panel compares nearby base sizes and
root layouts on fresh points under two independent Tier 1 strata. Its
candidate selection remains separate from this two-size precision result.
