# N83 same-point Tier 1 online result for the two finalists

Both complete N83 pipelines solved the same previously unseen public subgroup point on one Tier 1 cloud allocation and passed independent scalar, factor-base, four-point relation, and noise replay. K3425/root4/standard took **10.792847 ms** online against **4,461.115133 ms** for its matched rho solve, a **413.340×** ratio. K3450/root11/wide3 took **35.623502 ms** online against **2,276.915256 ms** for its matched rho solve, a **63.916×** ratio. The two rho solves used independently frozen seeds; each ratio belongs to its own same-point pair.

| Complete pipeline | Actual usable base points | Folded columns | IC online ms | Matched rho online ms | Rho / IC | IC setup s | IC peak GiB | Order |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| K3425/root4/standard | 568,550 | 3,425 | **10.792847** | 4,461.115133 | **413.340×** | 614.681 | 14.425 | IC, rho |
| K3450/root11/wide3 | 572,700 | 3,450 | **35.623502** | 2,276.915256 | **63.916×** | 552.886 | 11.503 | rho, IC |

The single public point is `[4738094823042181344100879, 7217856775769869847588318]` on `EC1N83Ckb1a1haacb69ee9f37`, whose prime subgroup order is `8569786107849059`. Its public-point SHA-256 is `dfea95ee4e423dd627adbd116828e9f478779c51476f0fa172cf4c160e645d12`; the verifier-only scalar sidecar SHA-256 is `ef4d939a30fbf26809f97acc0607fcd4fe04661e658e8b353c1ec77f5075d335`. The workload ID is `d20e46bcef1e`. Each IC process constructed an empty index and factor-log system before its timed target interval. Those cold setup costs are retained separately. The online clock starts at the first target-dependent operation and stops after producer scalar verification; independent replay follows outside that clock.

| IC online phase, ms | K3425/root4/standard | K3450/root11/wide3 |
| --- | ---: | ---: |
| Target query | 0.000030 | 0.000040 |
| Point decomposition | 10.655241 | 35.502725 |
| Relation check | 0.022070 | 0.005161 |
| Target descent | 0.000220 | 0.000310 |
| Scalar recovery check | 0.114786 | 0.114745 |
| Timing boundary overhead | 0.000500 | 0.000521 |
| **Charged online total** | **10.792847** | **35.623502** |

The GitHub-hosted Ubuntu 24.04 allocation exposed four logical AMD EPYC 7763 CPUs on two cores and 16 GB host memory. The runner pinned both pairs to those CPUs. Both calibration triples passed: after/before median ratios were 0.967789 and 1.025547; the largest within-triple spread was 1.106353. The selected CPUs recorded zero steal ticks, and each arm recorded zero new cgroup CPU-throttled periods, OOM events, and OOM kills. CPU frequency was not exposed. This is a **Tier 1 cloud-instance** measurement under its observed resource envelope.

The one-target IC times are descriptive for factor-base choice. The earlier fresh 4,000-point confirmation placed K3425 0.398616 ms ahead by mirrored mean, with a two-sided 99% paired interval of [-2.103728, +1.322782] ms. Its preregistered factor-base decision remains **`BOUNDED_UNKNOWN`** between K3425/root4/standard and K3450/root11/wide3. The same-point result supplies a verified online IC/rho comparison for each finalist; a single point cannot resolve the factor-base population ranking.

## Frozen evidence and replay

[Cloud run 38818664674](https://github.com/aburan28/crypto/actions/runs/38818664674) executed source commit `e1a77618587533364cf2236f09cc41b5bc816429`. The IC binary SHA-256 is `decd62935e518ff769cb1198b9da4d368f7b98ec8169d4f9e73c15d0e902ece9`; the rho binary SHA-256 is `7d17330ab29e753d5b3aad20b886951a9515d64ac78276a10401cef69f226380`. Artifact `12013488144` ZIP SHA-256 `73547413c41f05450e60274e695b94a38b574e10c1be9cffdaf98cf0e6810778` matches GitHub's artifact digest. The cloud analysis SHA-256 was `c6127202584a3410dc81382af9fbb2f0088dd1a7c4ab111ca62e63c6a5174f26`; separate local replay of the archived raw data regenerated that same SHA-256 and returned `PASS_TIER1_SAME_POINT_PAIR_SET`.

The raw ZIP and extracted receipts are archived at `/Volumes/SSD990-2/n83-archive-20261015/raw/n83-unresolved-tier1-38818664674(.zip)`. `records_index.json` pins that artifact, the source and host receipts, one compressed workload, and both compressed paired measurement rows. The current-head replay workflow downloads the exact artifact, independently verifies the frozen executed source and point, and regenerates the sealed records.
