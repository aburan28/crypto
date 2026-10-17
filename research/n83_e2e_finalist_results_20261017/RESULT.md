# N83 factor-base frontier on sixteen fresh Tier 1 points

The K1107 factor base lowers one-target **cold** process cost against the
same-source K1120 control by 3.949 seconds in root11 layout (paired mean;
95% interval [3.404, 4.494] seconds; exact two-sided sign-flip
`p = 0.000030518`). Root22 lowers it by 4.332 seconds
([3.318, 5.346] seconds; the same `p`). All 128 IC processes, 128
same-point rho solves, four IC repeat controls and 256 before/after rho
calibrations passed the two Tier 1 host gates and independent local replay.
The frozen eight-configuration selector records `BOUNDED_UNKNOWN` for both
online and cold choice: the K1107 layouts are indistinguishable at its
Bonferroni threshold, and online differences among factor bases do not pass
that threshold.

K1120/root11 has the lowest descriptive root11 online median, 73.487 ms,
and a median paired same-point rho/IC online ratio of 22.355 over the
sixteen points. It is the predeclared fifth arm in the next [five-arm
challenge](../n83_smallk_challenge_20261017/PROTOCOL.md), alongside the new
K1248 and K1280 factor bases and the K2543 and K3072 controls. The primary
online comparison begins after reusable base, index and log preparation;
the cold process column separately charges that preparation.

| Exact configuration | Usable points | Folded columns | Online median, ms | Cold mean, s | Paired rho/IC online median | IC faster on points |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| K1104/root11 | 183,264 | 1,104 | 113.394 | 108.382 | 16.247 | 16/16 |
| K1104/root22 | 183,264 | 1,104 | 111.651 | 108.112 | 11.334 | 16/16 |
| K1107/root11 | 183,762 | 1,107 | 79.820 | 104.826 | 18.507 | 16/16 |
| K1107/root22 | 183,762 | 1,107 | 82.929 | 104.635 | 23.910 | 16/16 |
| K1112/root11 | 184,592 | 1,112 | 195.899 | 106.596 | 12.457 | 15/16 |
| K1112/root22 | 184,592 | 1,112 | 191.117 | 106.699 | 15.636 | 16/16 |
| K1120/root11 | 185,920 | 1,120 | 73.487 | 108.775 | 22.355 | 16/16 |
| K1120/root22 | 185,920 | 1,120 | 71.913 | 108.966 | 19.394 | 16/16 |

For the canonical progress chart, K1120/root11's median **paired**
IC/rho online ratio is `0.0450834882`. K1107/root11's mean cold process
time divided by its mean paired rho **online** time is `47.06115645`;
that hybrid cold/reference ratio is supplementary and names both timing
boundaries explicitly. Its mean paired rho online time is 2.227432 s.

The [graph](finalist.svg) and [two-page PDF report](RESULT.pdf) show every
point, median and interquartile range
for online IC time, cold process time and the independently paired rho/IC
ratio. Its [source](build_visual.py) reads only the provider-bound
[validated rows](../n83_e2e_final_holdout_20261016/validated_rows.json) and
[frozen selection](../n83_e2e_final_holdout_20261016/selection.json); the
underlying chart numbers are in [plot_data.json](plot_data.json). The exact
eight [candidate manifests and 132 IC1 rows](../n83_e2e_final_holdout_20261016/)
retain every per-point phase, rho and resource record. Each pairwise online
comparison uses the mean log time ratio over the same sixteen points;
cold comparisons use paired process milliseconds. The selector tests all
28 unordered pairs with exact `2^16` sign flips and
`0.05/28 = 0.001785714` per-pair alpha. Its smallest online `p` is
0.045258. Both K1107 layouts beat every other tested K in cold paired mean
at this threshold, while their mutual cold comparison has `p = 0.617950`.

## Instance, workload and evidence

The curve is `EC1N83Ckb1a1haacb69ee9f37`, polynomial-basis `GF(2^83)`
with defining-polynomial integer `9671406556917033397649557`, equation
`y² + xy = x³ + x² + 1`, prime subgroup order `8569786107849059` and
cofactor `1128547018`. The [sealed protocol](../n83_e2e_final_holdout_20261016/COMPARISON_PROTOCOL.md)
and [fixture](../n83_e2e_final_holdout_20261016/PROTOCOL.md) fixed sixteen
distinct public points before the finalists were selected. Target file
SHA-256 is `407cd248797100214ff039193bb04bf89b5c2708a8cd339ca98f7aca98174062`;
the verifier-only answer sidecar SHA-256 is
`fe771284a3b0e5c0301309f678751f0630f0fa23844dd0ce1dbc9fa47f243076`.
Each candidate-point arm used a fresh process and its own same-point rho
solve. Source SHA-256 is
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`;
the executed IC and rho binary SHA-256 values are
`3e0a5e99ba0fd7cf507bfa129b50ee2a39c80d338261413e87c6ded69e18f90f`
and `4ed1335f55240ff98470caf0096b9c0cbb3da0a01b157ce88bc83cf886f761e9`.

[GitHub Actions run 39005111289](https://github.com/aburan28/crypto/actions/runs/39005111289)
executed source commit `416489e09c3d8885b6c12f262935e6ba29373460`.
Half 0 used runner `GitHub Actions 1000326160`; half 1 used
`GitHub Actions 1000326147`. Each exposed four AMD EPYC 9V74 vCPUs, two
cores with sibling threads, and respectively 16,373,452 and 16,372,440
KiB of host memory. Both
recorded `PASS_TIER1`, `PASS_10_PERCENT_PER_PAIR`, and
`PASS_INDEPENDENT_REPLAY`. The host, cgroup, affinity, throttling, steal,
memory, source, command and event-order receipts are in the raw archive.
Their provider artifact IDs are `12095116702` and `12094374158`; exact ZIP
SHA-256 values are
`50db2c8e06c884ad9f6a2731aa350dc8bde708d8ea0c52fddf740404874908a4`
and `737df73989d1a2e663ef3ea9253bb088b6aa92d2734070cae3a3d4b641b75c20`.
Separate local verification reproduced both producer analyses as exact JSON
records.

The [deduplicated raw archive](raw_dedup.tar.gz) retains every one of 3,802
file paths and hashes, storing the 132 repeated base dumps as four distinct
blobs. Its [manifest](raw_dedup_manifest.json) and
[receipt](raw_dedup_receipt.json) bind the original provider ZIP hashes.
The 23,533,732-byte archive SHA-256 is
`b057954c12a969591a3d4451eb2a693b5f7ec6cbb611905a4c686733abedd422`.
The [reconstruction checker](verify_archive.py) restored and verified all
3,802 paths in a fresh directory. The two original provider ZIPs and
extractions remain under
`/Volumes/SSD990/crypto/artifacts/n83_e2e_finalist_20261017`.
The [validation record](VALIDATION.md) retains the complete serial Rust
suite outcome, earlier capped and sandboxed attempts, example tests, Python
replay gates, and generated-site checks.

## Next measurement

The [ten-point five-arm protocol](../n83_smallk_challenge_20261017/PROTOCOL.md)
compares K1248/root11 and K1280/root11 against K2543/root11,
K3072/root11 and the predeclared K1120/root11 fallback on the same source.
It keeps one public point per fresh process, an independent same-point rho
solve per arm, two Tier 1 host strata, alternating order and before/after
rho calibration. Its [dispatch](../n83_smallk_challenge_20261017/dispatch.json)
binds this result's validated rows, selection, materialized IC1 rows and
provider run before the cloud workflow can start.
