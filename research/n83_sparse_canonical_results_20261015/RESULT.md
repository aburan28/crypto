# N83 sparse canonical selector accelerates root11 decomposition

The sparse linear-map canonical selector reduced K3300/root11 mean target
online time from **39.722 ms** to **31.963 ms** on 1,000 paired N83 points.
The paired sparse-minus-83-rotation-scan mean was **−7.759 ms** with the
predeclared 95% contiguous-block bootstrap interval **[−8.379, −7.125]
ms**. Both arms used the same 547,800-point base and 903,852,075-key exact
root table. All 2,000 target records and both complete bases passed
independent scalar and four-point relation replay under the 14 GiB guard.

| Selector | Mean target online ms | Setup s | Peak RSS GiB | Status |
| :--- | ---: | ---: | ---: | :--- |
| 83-rotation scan | 39.722 | 651.971 | 10.052 | `PASS_INDEPENDENT_REPLAY` |
| Sparse linear map | 31.963 | 446.585 | 10.053 | `PASS_INDEPENDENT_REPLAY` |

The mean online saving is **19.53%** relative to the scan arm; the measured
setup saving is **205.387 s**. Both arms used fresh processes with empty
base, root index and factor-log state. The target probe count was identical
at 55,787.783 mean; the selector changes canonicalization cost per probe.
The target interval includes query, four-summand lookup, relation check,
descent and recovery check. The predeclared statistic uses 20 contiguous
blocks of 50 targets, 10,000 bootstrap resamples and seed `2026101508`;
the upper endpoint below zero meets the frozen sparse-selector gate.

This is a cloud-instance stage comparison on the reused 1,000-point
refinement panel, recorded as `cloud_instance_exploratory`. The exact sparse
candidate is
`IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0h064d6bb97962`.
The curve is `EC1N83Ckb1a1haacb69ee9f37`, polynomial-basis `GF(2^83)`
with modulus bitmask `9671406556917033397649557`, model
`y²+xy=x³+x²+1`, subgroup order `8569786107849059` and cofactor
`1128547018`. The base hash is
`a4ba33098d2776520b98a326585a2301f8ccd5bc005df7c1a70fbe8496c61d8b`.

## Receipts

Executed source commit `39f8c5e74b2916333c6be60f2e653d7651eb7308` has
IC source SHA-256
`c9824acebda0360416c3c507594bb07c17ca798c32cbf8f84bba95ca2c44d5b1`.
[GitHub Actions run 38676235613](https://github.com/aburan28/crypto/actions/runs/38676235613)
completed the source tests, producer and independent verifier. Raw artifact
`11947421783` has SHA-256
`be12dbe701fbe33092b2bb104d2e4cd07d321374596b1f8c69b83c332ac137fc`;
its analysis has SHA-256
`6d555f4d6f8596fda5451b18e9f6cc0451fc63b2f792fd7baf1f33cf1fcf501c`.
The [selection](selection.json), two canonical candidate manifests,
[measurement rows](measurements.json.gz) and [digest index](records_index.json)
seal the exact result. The shared 1,000-point
[workload](../n83_factor_base_refine_results_20261014/workload.json.gz) has ID
`0d7b05d9fca6` and public-point SHA-256
`35e34151dcc1ad6048939560adf28fc7c1b1017833fe005a656c3f4e2bc07de8`.
