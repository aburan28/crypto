# N83 root22 index improves the K1104 cold process on the exposed point

At K1104, the dynamic root22 index reduced the complete fresh-process IC
interval by **1.640 s** (1.37%) against the dynamic root11 index on the same
public point and cloud allocation. The symmetric A/B schedule measured
root11 at 119.446 s and root22 at 117.806 s. Both within-layout controls
passed the frozen 5% spread gate. K1120/root22 had the lowest observed cold
interval in this screen, **117.258 s**; its fresh-point comparison remains the
next selection experiment.

| Pipeline | Usable base points | Folded columns | Cold process, s | Setup, s | Rank stage, s | One-target online, ms | Peak RSS, MiB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| K1104/root11 dynamic, A/B mean | 183,264 | 1,104 | 119.446 | 48.431 | 70.959 | 52.122 | 1,936.2 |
| K1104/root22 dynamic, A/B mean | 183,264 | 1,104 | 117.806 | 48.304 | 69.448 | 50.658 | 2,466.0 |
| K1120/root11 dynamic | 185,920 | 1,120 | 118.017 | 49.174 | 68.816 | 23.690 | 1,952.0 |
| **K1120/root22 dynamic** | **185,920** | **1,120** | **117.258** | **49.240** | **67.991** | **22.995** | **2,496.5** |

Root22 saved 1.511 s in the K1104 rank stage and 0.127 s in setup on the A/B
means. It used about 530 MiB more peak RSS than root11 at K1104. The two
layouts reconstructed the identical K1104 base digest
`62679ca766bbac029175604ef0af8368b69edaebdc8c160eba262d1413ca9cc9`;
the K1120 layouts shared digest
`4771fee7fd73c98e15198dbc997423b224f3e1f75c675f1cb2e985d8c11ec7f4`.
Each layout had the same root-table entry count at fixed K. The K1104 A/B
cold intervals were 119.584/119.308 s for root11 and 117.861/117.751 s
for root22, giving 0.231% and 0.094% within-layout spreads. The A/B pair
differences were 1.722 s and 1.557 s. K1120 had one process per layout;
its 0.759 s layout difference is a screen observation.

## Reproducibility and candidate identities

The [frozen protocol](PROTOCOL.md) ran six fresh IC processes in the order
K1104/root11 A, K1104/root22 A, K1120/root22, K1120/root11,
K1104/root22 B, K1104/root11 B, followed by one same-point rho process.
Every IC arm completed base, four-point relation, scalar, resource, topology,
telemetry, execution-order and exclusive-phase replay. The rho scalar replay
passed; its online interval was 721.786 ms. The cloud instance exposed four
AMD EPYC 7763 CPUs and 15.61 GiB host memory. The analysis status is
`SCREEN_COMPLETE` with `PASS_INDEPENDENT_REPLAY` and a null promoted selection.

The [IC1 row record](ic1_rows_38976053815.json) binds six run rows on workload
`eebf98ee1031` to four immutable candidate manifests. The root11 identities
remain the earlier K1104 and K1120 candidates; the new root22 identities are
`IC1N83Ckb1a1fb183264PDP4rootRCguidedLAgaussTDdirectISO0hce7cffedf9c0`
and
`IC1N83Ckb1a1fb185920PDP4rootRCguidedLAgaussTDdirectISO0h82b74c23b7bc`.
The source SHA-256 is
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`.

[GitHub Actions run 38976053815](https://github.com/aburan28/crypto/actions/runs/38976053815)
completed at exact execution head `cec93c657f49441f59d6408a3a81038f72a6ee97`.
The provider artifact `n83-e2e-crosslayout-38976053815-1` has ID
`12080884073` and ZIP SHA-256
`5556c506faaff4f54baac31010da3f0abde0747410045828acdcbb983ef6cb6b`.
The archived analysis SHA-256 is
`af48f1cc555b6c65bed99f5426dba04c0b5053f039a2003169110b31cd08a3e1`.
An independent local verifier exited `0` with empty stderr and reproduced
the parsed analysis; the local stdout SHA-256 is
`a18c1927b433500754ca9f156abfe42c7a7b3dacbe494421b5de098d74c2cf02`.
The materializer checks the provider ZIP, all 60 extracted result members,
host/source identity, replay, and canonical candidate hashes.
The raw ZIP and extracted receipts, including the independent replay logs,
are retained under
`/Volumes/SSD990-2/n83-archive-20261016/{raw,extracted}/n83-e2e-crosslayout-38976053815`;
the raw file has the `.zip` suffix.

Next, compare the K1120/root22 candidate with the best K1104/K1120 root11
candidate from the separate 16-point precision panel on fresh public points,
with repeated A/B controls and alternating same-point rho pairs under one
Tier 1 resource envelope. The adjacent-K screen may nominate another K for
that holdout; its finalist must be frozen before the fresh-point run.
