# Five-point N83 cold holdout confirms two smaller bases against K1152

K1120 saved **2.645 s** of one-target cold process time against K1152 on five
sealed public points (paired 95% interval **[2.085, 3.205] s**). K1104 saved
**2.244 s** (**[1.980, 2.508] s**). K1120 had the lowest observed mean,
**119.826 s**, but its **0.401 s** advantage over K1104 has a paired 95%
interval **[−0.236, 1.038] s**. The frozen winner rule therefore records
`BOUNDED_UNKNOWN` and `selected_candidate_id: null`. All 25 IC arms and five
same-point rho arms passed independent cloud replay, resource checks, and
the five K1152 A/A gates.

| Pipeline | Actual usable points | Folded columns | Mean cold process, s | Setup, s | Rank, s | Target, s | Mean target online, ms | Peak RSS, MiB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| K1104/root11 dynamic | 183,264 | 1,104 | 120.227 | 48.728 | 71.320 | 0.176 | 175.766 | 1,937.5 |
| K1108/root11 dynamic | 183,928 | 1,108 | 121.296 | 48.955 | 72.171 | 0.167 | 166.802 | 1,942.7 |
| K1120/root11 dynamic | 185,920 | 1,120 | **119.826** | 50.081 | 69.525 | 0.216 | 215.955 | 1,954.5 |
| K1152/root11 dynamic A/B mean | 191,232 | 1,152 | 122.471 | 53.067 | 69.255 | 0.144 | 144.047 | 1,990.8 |

The cold process interval starts with empty base, index, and log state in a
fresh child and includes construction, relation rank, one public target,
and producer verification. K1152's two child times are averaged separately
on each point before the five-point mean. The target-online interval starts
after reusable setup and logged base preparation. Same-point rho's five
target-online intervals were 2.132, 1.359, 2.892, 1.926, and 2.512 s;
their mean was 2.164 s. These intervals and every phase are retained in
the raw run rows.

| Paired cold saving | Mean, s | 95% interval, s | All five point differences positive? |
| --- | ---: | ---: | --- |
| K1120 versus K1152 | 2.645 | [2.085, 3.205] | yes |
| K1104 versus K1152 | 2.244 | [1.980, 2.508] | yes |
| K1104 versus K1108 | 1.069 | [0.734, 1.405] | yes |
| K1120 versus K1108 | 1.470 | [0.716, 2.225] | yes |
| K1120 versus K1104 | 0.401 | [−0.236, 1.038] | no |

The five K1152 A/A relative spreads were 0.588%, 0.536%, 0.038%, 0.651%,
and 0.009%; each passed the prespecified 10% gate. Topology, CPU counter,
execution-order, resource, and cross-point exact-base identity checks passed.
The allocation was a GitHub-hosted Linux x86-64 cloud instance exposing four
AMD EPYC 7763 CPUs and 16,373,452 KiB of memory. Its observed affinity and
cgroup envelopes are in `host.json`; the measurement tier is
`cloud_instance_tier1`.

## Frozen instance and evidence

The curve is `EC1N83Ckb1a1haacb69ee9f37` over polynomial-basis
`GF(2^83)`, defining-polynomial bitmask `9671406556917033397649557`, model
`y² + xy = x³ + x² + 1`, prime subgroup order `8569786107849059`, cofactor
`1128547018`. The five points were committed before screen timing, with
point-file SHA-256
`722bc60cd010e533cd715a9bfb321fba52d950352c99684cd422b23af7f93d26`
and verifier-only answer-file SHA-256
`b4c2ae858770b7c8842703558984326bc0002b0f1deb3c898d5ac30562de8b96`.
The candidate screen and holdout rules are in the [protocol](PROTOCOL.md).

[GitHub Actions run 38955665283](https://github.com/aburan28/crypto/actions/runs/38955665283)
executed commit `e44c5022519e190cfa72fad326e1ba399d49893d` with IC source
SHA-256 `8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`.
Provider artifact `12073412864` is archived at
`/Volumes/SSD990-2/n83-archive-20261016/raw/n83-e2e-holdout-38955665283.zip`;
its SHA-256 is
`d2e73fff94e3dc972c47be646667eef7ab2fa4a20f661f44f02323b25d8ff060`,
equal to GitHub's artifact digest. ZIP integrity and path checks passed for
all 255 members. The cloud `analysis.json` SHA-256 is
`b436fcad4ad5ae8af8b773c1dcbfb8c043fda12433fa509ae49639cb77f33b89`.
Separate local independent replay of all 25 IC arms and five rho arms exited
zero. Its regenerated `analysis.json` is byte-identical to the cloud ZIP;
its stdout SHA-256 is
`b3efa6240fe6188e4d4f746bcd6b596469b837ac2c935d16a10ec5e31f2aa2b3`
and stderr is empty. The exact [IC1 row record](ic1_rows_38955665283.json),
SHA-256 `4c296ed1a4e145cfaeff850b0a8882c5721971f607f2be2b4119d7e210c85600`,
contains four candidate manifests, five workload IDs, and 25 complete IC
rows. K1120 is
`IC1N83Ckb1a1fb185920PDP4rootRCguidedLAgaussTDdirectISO0hc922b499a234`;
K1104 is
`IC1N83Ckb1a1fb183264PDP4rootRCguidedLAgaussTDdirectISO0h6f30a77abf9e`.
The archive and local replay logs remain under
`/Volumes/SSD990-2/n83-archive-20261016/`.

The next precision experiment is a fresh, sealed K1104/K1120 matched panel
with mirrored order and an A/A bracket on each point. Its own frozen
confidence rule should resolve the 0.401 s comparison under a recorded
Tier 1 envelope. Other base geometries and index layouts remain separate
candidate families with their own matched cold comparisons.
