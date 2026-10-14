# N83 guarded root-only factor-base frontier

The optional `KIC_ROOT_ONLY22=1` mode stores only two 83-bit roots per
symmetric S3 state. The state labels are recovered from the deterministic
`(left,right,relative)` traversal ordinal after an exact root-key hit. It
asserts that every retained pair produced S3 roots and that the final state
count equals the full loop bound. A missing root stops the arm with a
recorded failure before any target claim. K1798, K2200, and K2543 prior
receipts each recorded the full bound; this gate tests the larger domain.

## Frozen instance and source

- Curve: `y^2+xy=x^3+x^2+1` over `F_2^83`, polynomial low terms
  `[0,2,4,7]`; subgroup order `8569786107849059`, cofactor `1128547018`,
  generator `(6614244651398184860292015,8868260410029249391539529)`.
- Public point: `research/n83_indexed5_20261013/primary_target.jsonl`,
  SHA-256 `21c8dfc7d488909e5c4a42c839fcc01bdefa9821c11ba001626a1d15697240e8`.
  The independent scalar sidecar SHA-256 is
  `651eb7d0ef08a1655036746311b19b6989e340ca4752b6d675606055e59c12a1`.
- Producer source SHA-256:
  `02ca5368352267fe4defca4e32a68ee5bd09968a1cb6e75d294c5f75d8599c60`.
  `Cargo.lock` SHA-256:
  `ac5cd7a3459be776180afe1eb2975c76b1a88f84bd1e99aa20475927ec260a90`.
- Every arm uses symmetric S3, exact five-byte root slots, the
  `ceil(4S/3)` request, no Bloom filter, no sparse canonicalization, rank
  seed 29, four rank threads, one target thread, and a fresh process.
  K and the root-only versus 26-byte layout vary as shown below.

| Order | Arm | K | Predicted full-orbit usable points | State bound | State plus slot payload GiB | Child cap |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 1 | root-only 22-byte | 3114 | 516,924 | 402,426,891 | 13.245 | 2400 s |
| 2 | root-only 22-byte | 3050 | 506,300 | 386,055,275 | 12.910 | 2100 s |
| 3 | root-only 22-byte | 2940 | 488,040 | 358,710,870 | 12.350 | 1800 s |
| 4 | root-only 22-byte | 2543 | 422,138 | 268,375,505 | 10.499 | 900 s |
| 5 | encoded-label 26-byte | 2543 | 422,138 | 268,375,505 | 11.499 | 900 s |

The bound is `83*K*(K-1)/2+42*K`. Lean checks that K3114 fits the
`2^30`-slot table while K3115 crosses to `2^31` slots, and that the K3114
two-array payload is below 13.3 GiB. This excludes allocator and other
process memory. The K2543 prior receipt measured approximately 0.297 GiB
of overhead above state and slot arrays; actual RSS is measured here. The
runner stops a child above 14 GiB RSS, keeps its stdout/stderr, `wait4`,
cgroup, wall-cap, and resource-stop records, and runs all later arms.

The same-point strong-rho arm uses archived `fold-L64-DP6-E0` config SHA-256
`cf76ebb5372acfc1355b88374f875980be66cdffba1a82cb6e4c2f28bcfa723b`,
seed `2026101301`, one thread, and a 120-second cap. The independent
verifier checks every completed base, subgroup, four-point relation,
scalar, exact state count, layout, and five exclusive online phases. At
K2543 it compares base hash, exact root-key count, and scalar across the
root-only and encoded-label layouts. The source test compares all six
real-field representatives and 2,047 probes for exact first witnesses.

This one-run cloud-instance gate tests correctness and capacity. Its
same-point timing ratios are exploratory until five alternating pairs pass
the Tier 1 noise gate. Actual usable counts and canonical `IC1` manifests
are materialized only for independently replayed arms. A new public point
will be frozen after the base is selected for the controlled online study.
