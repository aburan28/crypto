# N83 symmetric S3 root-index factor-base gate

The optional symmetric-S3 state quotient keeps one representative of each
pair of indexed S3 inputs related by exchange and Frobenius rotation. It
retains the earlier source-loop state, so the exact root table can still
return the first witness for every canonical key. This gate tests whether
the resulting memory reduction admits a larger complete N83 factor base on
one 16 GB cloud instance and measures the same public target against the
full-state index and strong rho.

## Frozen mathematics, inputs, and source

- Curve: Koblitz `y^2+xy=x^3+x^2+1` over `F_2^83`, polynomial low terms
  `[0,2,4,7]`, subgroup order `8569786107849059`, cofactor `1128547018`,
  generator `(6614244651398184860292015,8868260410029249391539529)`.
- Public point: `research/n83_indexed5_20261013/primary_target.jsonl`,
  SHA-256 `21c8dfc7d488909e5c4a42c839fcc01bdefa9821c11ba001626a1d15697240e8`.
  The worker receives only this point. The separate correctness sidecar is
  `answer.jsonl`, SHA-256
  `651eb7d0ef08a1655036746311b19b6989e340ca4752b6d675606055e59c12a1`.
- Source file: `examples/koblitz_orbit_dlp_fast.rs`, SHA-256
  `e3a37208fc72fe809a8bbc2a1604f0d59f8a6582a9deaa3fd58a077083bf2c3f`.
  `Cargo.lock` SHA-256 is
  `ac5cd7a3459be776180afe1eb2975c76b1a88f84bd1e99aa20475927ec260a90`.
  The workflow checks all source and input hashes before the run; its host
  receipt records the executed commit and binary hashes.
- Shared flags: exact five-byte root index, 28-byte N83 state layout,
  twofold root-table request, Bloom/dense/load75/sparse-canonical off,
  rank seed 29, four rank threads, one target thread. Only
  `KIC_SYMMETRIC_S3` and K vary between arms. Each arm starts a fresh
  process with an empty base, index, and factor-log table.

## Ordered one-target arms and bounds

| Order | Arm | K | Predicted usable points if all orbits are full | State-count upper bound | Root-slot bucket from bound | Wall cap |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 1 | symmetric S3 | 2,543 | 422,138 | 268,375,505 | 5 GiB | 1,800 s |
| 2 | symmetric S3 | 2,200 | 365,200 | 200,861,100 | 5 GiB | 1,500 s |
| 3 | symmetric S3 | 1,798 | 298,468 | 134,162,265 | 2.5 GiB | 1,200 s |
| 4 | full states | 1,798 | 298,468 | 268,322,732 | 5 GiB | 1,800 s |

The quotient bound is `83*K*(K-1)/2 + 42*K`; the full bound is `83*K^2`.
Both are source-loop limits, not measured regular-state counts. Lean
checks that K2543 stays in the 5-GiB slot bucket, K2544 crosses it, and
the K2543 state plus slot array payload is below 12 GiB. Actual usable
points, root keys, allocated bytes, peak RSS, and online cost come only
from the run receipt and independent verifier.

The runner stops a child if sampled RSS exceeds 14 GiB or its wall cap
expires. It records wall/user/system CPU time, RSS, major faults,
involuntary context switches, cgroup memory and CPU counters, every raw
stdout/stderr/exit status, base dump, and target row. It attempts all arms
even if an earlier arm stops; timeout, OOM, guard stop, or verification
failure remains a censored or failed row with its exact limit.
The workflow has a 135-minute overall cap, above the 107 minutes of summed
child caps, leaving time for compilation, independent replay, and artifact
upload even when every child reaches its cap.

The same-point strong-rho reference uses the archived
`fold-L64-DP6-E0` policy, seed `2026101301`, one OS thread, and a 120 s
cap. Its config SHA-256 is
`cf76ebb5372acfc1355b88374f875980be66cdffba1a82cb6e4c2f28bcfa723b`.

The independent verifier replays the complete point-defined factor base,
subgroup membership, four-point relation, recovered scalar, and five
exclusive online phases. It checks the quotient state bound and layout;
for K1798 it compares the quotient's exact root-key count with the full
index. The full K1798 arm is also compared with the archived first-witness
row. The source unit test directly compares every real-field canonical
key and first witness in both root-table layouts for six representatives.

Each verified arm gets a distinct canonical `IC1` candidate manifest and
one-target workload/run ID. The primary online interval begins with
target-dependent query work after reusable base, index, and factor logs are
ready, and ends after scalar verification. Setup and rank are recorded
separately. One run per arm is an instance-scoped Tier 1 capacity and
correctness gate; timing ratios from this gate are exploratory. A
controlled Tier 1 wall-time claim requires five alternating same-point
pairs and the tiered service's noise gates. A later 64-point panel and a
256-point finalist panel are separate coverage questions, each retaining
one-target rows rather than shared-target averages.

The branch-scoped workflow uploads its raw artifact on success or failure.
It does not alter repository-wide Actions settings or expose credentials
to the benchmark children.
