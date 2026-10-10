# N53 C1800 prefix plus prefetch and projective target checks

The combined S3 target path recovered all eight frozen N53 single-target
scalars in two fresh processes per arm and reduced the local geometric mean
online wall time by **3.247×** against the C1200 packed reference. The
target-resampled 95% bootstrap interval is **[1.868, 5.843]×**. The local
point-estimate gate of at least 2× passed; the stricter lower-bound gate of
at least 2× did not. This host has no exclusive CPU partition receipt, and
both arms' reusable setup times moved outside the frozen 10% noise band, so
the controlled speedup field remains `null`.

The method searches the existing C1200 state order first, then the complete
C1800 state index. It pairs that schedule with late root-slot prefetch and
portable projective scalar checks. The larger factor base has **190,800**
distinct subgroup-usable points before sign/Frobenius folding and **1,800**
relation columns; the reference has **127,200** points and **1,200** columns.
Both reached complete pre-target rank, used one target worker, and incurred
zero target-charged relation-collection fallback attempts. Factor-base size
and the online kernel both change in this comparison.

| Public target | C1200 online ms, two repeats | C1800 combined online ms, two repeats | C1200 / C1800 geometric ratio | Target S3 calls, C1200 / C1800 |
| --- | ---: | ---: | ---: | ---: |
| T205 | 2.271, 1.712 | 0.154, 0.106 | 15.409× | 3,104 / 128 |
| T206 | 0.559, 0.443 | 0.111, 0.086 | 5.103× | 508 / 128 |
| T207 | 0.255, 0.287 | 0.276, 0.093 | 1.690× | 116 / 128 |
| T208 | 0.330, 0.268 | 0.115, 0.882 | 0.933× | 188 / 160 |
| T209 | 1.583, 1.794 | 0.295, 0.311 | 5.565× | 2,670 / 560 |
| T210 | 0.373, 0.627 | 0.131, 0.723 | 1.574× | 212 / 224 |
| T211 | 0.650, 8.554 | 0.284, 2.233 | 2.959× | 656 / 624 |
| T212 | 2.042, 2.463 | 1.055, 0.322 | 3.850× | 3,142 / 672 |

The online interval starts with the public target available and the reusable
index and base logs ready. It includes target query, charged target PDP,
relation check, descent, and scalar recovery check. Point generation, input
loading, reusable setup, and checked Sage replay are outside that interval.
The geometric mean C1200/C1800 target S3-call ratio is **2.628×**, with a
target-bootstrap interval of **[1.307, 5.818]×**. S3-call counts are a
deterministic stage diagnostic; wall time includes the complete online path.

The [repeat-aware robustness analysis](robustness.json) resamples both target
identities and the two timings per arm without dropping any run. Its 95%
interval is **[1.595, 6.429]×**; omitting any one target leaves a geometric
ratio between **2.599× and 3.881×**. The empirical log-ratio variance has
0.683 between-target and 0.332 within-target-repeat components. If the
observed effect and both variance components persisted on a qualified host,
a normal approximation calls for about **17 distinct targets** for an
expected two-sided 95% lower bound above 2× and **34 distinct targets** for
80% power, with two runs per arm and target. These are planning estimates
from exploratory data; the Modal pilot and host-isolation preflight must
establish the variance and fidelity for the actual confirmation host.

[`CONFIRM_FREEZE.json`](CONFIRM_FREEZE.json) fixes the two candidate IDs,
source and binary hashes, all eight public points and workload IDs, the
SHA-256 scalar input law inherited from the parent freeze, ABBA/BAAB order,
14 setup workers, one target worker, two repetitions per arm and point,
180-second per-run cap, and 32 GiB observed RSS envelope. The combined
source's T025 control reached rank 1,800, recovered its frozen scalar, and
observed 23,356,014,592 bytes peak RSS before the fresh panel. Across the
panel, all **32/32** native runs passed phase accounting, exact scalar,
full-rank, and memory checks; the largest observed RSS was
25,172,295,680 bytes. The checked Sage launcher independently replayed all
32 scalars against the public points and verified the generator order.

The raw paired ledger, compressed reports, stdout, stderr, host record,
resource counters, and independent Sage receipt are in
[`local/paired_confirm`](local/paired_confirm). The ledger SHA-256 is
`21e3548722e1b442b7c878cac063bb2f8a1babe0cd5e63a5c7bd6887a4923344`.
[`analysis.json`](analysis.json) retains every paired phase and uncertainty
calculation. The baseline and combined reusable setup medians were 14.355 s
and 36.406 s; neither arm passed the 10% setup-drift check. For example,
T208's 0.882 ms combined repeat had six major page faults and substantially
more involuntary context switches than its 0.115 ms repeat. All such runs
remain in the analysis. The exact sources and workload can now be taken to
an auditable isolated host for the >=2× controlled confirmation.

The [bounded Modal T205 pilot](modal/README.md) now has exact Linux binaries,
an image-transfer manifest, a public-only point bundle, guest fidelity
telemetry, and a passing offline preview. It awaits authorization for the
source-derived transfer and Sandbox creation.
