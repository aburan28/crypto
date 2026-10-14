# N83 one-root state capacity and online gate

The S3 solver returns two normal-basis roots per symmetric state. Its regular
branch sets the second root to the first plus `q`, where
`q = xy/(x+y)^2` for the two shifted representative coordinates `x,y`.
When `x+y=0` or `xy=0`, its exceptional branch returns equal roots.
The optional `KIC_ROOT_ONLY11=1` layout retains the first 83-bit root in
11 bytes and reconstructs the second root from the state ordinal and the
shifted representatives on an exact-key comparison or query. The index
inserts both roots while each state is built, so ordinary construction does
not require a second pass of reconstructed roots. Exact full-key equality
and first-inserted witnesses remain the table policy.

## Frozen instance and variants

- Curve: `y^2+xy=x^3+x^2+1` over `F_2^83`, polynomial low terms
  `[0,2,4,7]`; subgroup order `8569786107849059`, cofactor `1128547018`,
  generator `(6614244651398184860292015,8868260410029249391539529)`.
- Public point: `research/n83_indexed5_20261013/primary_target.jsonl`,
  SHA-256 `21c8dfc7d488909e5c4a42c839fcc01bdefa9821c11ba001626a1d15697240e8`.
  The scalar sidecar SHA-256 is
  `651eb7d0ef08a1655036746311b19b6989e340ca4752b6d675606055e59c12a1`.
- Producer source SHA-256:
  `07ce5094f0e561a75bd65ca143f2301830447742268936805c0b920f7a93f6ae`.
  `Cargo.lock` SHA-256:
  `ac5cd7a3459be776180afe1eb2975c76b1a88f84bd1e99aa20475927ec260a90`.
- All arms use symmetric S3, exact five-byte root slots, a fixed `2^30`-slot
  table, no Bloom filter, no sparse canonicalization, rank seed 29, four
  rank threads, one target thread, and a fresh process. The K2543 arms
  compare 11-byte reconstructed roots with 21-byte retained root pairs on
  the same source, factor base, root keys, target, and table capacity.

| Order | Arm | K | Predicted usable points | Complete state bound | Array payload GiB | Child cap |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 1 | one root | 3400 | 564,400 | 479,741,700 | 9.915 | 1800 s |
| 2 | one root | 3300 | 547,800 | 451,936,650 | 9.630 | 1800 s |
| 3 | one root | 3200 | 531,200 | 424,961,600 | 9.354 | 1500 s |
| 4 | one root | 2543 | 422,138 | 268,375,505 | 7.749 | 900 s |
| 5 | retained pair | 2543 | 422,138 | 268,375,505 | 10.248 | 900 s |

The complete state bound is `83*K*(K-1)/2+42*K`. At K3400 twice that
bound is `959,483,400`, below the fixed slot count `1,073,741,824`;
the table may nevertheless suffer long probes near its 89.36% worst-case
load. The array payload is `11*S+5*2^30` bytes. It excludes allocations,
shifted representative tables, and runtime memory, which the runner measures
with sampled RSS and final `wait4` peak. The runner stops a child above
14 GiB sampled RSS and preserves the final peak if it exceeds the cap
between samples. Each arm retains stdout, stderr, usage, cgroup, and wall-cap
receipts, including failures and subsequent fallback arms.

The same-point strong-rho arm uses archived `fold-L64-DP6-E0` config SHA-256
`cf76ebb5372acfc1355b88374f875980be66cdffba1a82cb6e4c2f28bcfa723b`,
seed `2026101301`, one thread, and a 120-second cap. The independent
verifier replays each complete base, subgroup, four-point relation, scalar,
exact state count, allocation, and five exclusive online phases. It checks
K2543 base hash, root-key count, and scalar parity between layouts.

The primary cloud run measures a complete single-target online interval
after setup on a public point used by previous capacity gates; its timing
is exploratory. The archived five-point panel is a secondary selection
experiment. A new point and five alternating same-point candidate/rho pairs
are the Tier 1 timing gate after factor-base selection.
