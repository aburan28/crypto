# N83 indexed five-byte dense-request factor-base gate

The optional `KIC_INDEXED5_LOAD75=1` mode requests `ceil(2S/3)` entries for
`S` retained symmetric-S3 states. The five-byte constructor allocates the
next power of two above twice that request. It keeps the same state and
canonical-key stream, exact full-key verification, and first witness. This
gate tests how far the larger usable factor base improves one-target online
cost within one 16-GB cloud instance.

## Frozen inputs and source

- Curve: Koblitz `y^2+xy=x^3+x^2+1` over `F_2^83`, polynomial low terms
  `[0,2,4,7]`, subgroup order `8569786107849059`, cofactor `1128547018`,
  generator `(6614244651398184860292015,8868260410029249391539529)`.
- Public target: `research/n83_indexed5_20261013/primary_target.jsonl`,
  SHA-256 `21c8dfc7d488909e5c4a42c839fcc01bdefa9821c11ba001626a1d15697240e8`.
  The worker receives only the point. The answer sidecar SHA-256 is
  `651eb7d0ef08a1655036746311b19b6989e340ca4752b6d675606055e59c12a1`.
- Producer source: `examples/koblitz_orbit_dlp_fast.rs`, SHA-256
  `c97ec7d8caef4c2a93b82378c5961a05d874e754f3873de64e7bcf2afcfd8e21`.
  `Cargo.lock` SHA-256 is
  `ac5cd7a3459be776180afe1eb2975c76b1a88f84bd1e99aa20475927ec260a90`.
- Shared flags: symmetric S3 state quotient, exact five-byte root index,
  28-byte N83 states, Bloom/dense/load75 packed and sparse-canonical off,
  rank seed 29, four rank threads, one target thread. Each arm starts a new
  process with an empty base, index, and factor-log table. Only K and
  `KIC_INDEXED5_LOAD75` vary.

## Ordered one-target arms

| Order | Arm | K | Predicted usable points if all orbits are full | State bound | Slot bucket | Child wall cap |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 1 | dense request | 3,114 | 516,924 | 402,426,891 | 2.5 GiB | 2,400 s |
| 2 | dense request | 3,000 | 498,000 | 373,501,500 | 2.5 GiB | 2,000 s |
| 3 | dense request | 2,543 | 422,138 | 268,375,505 | 2.5 GiB | 1,400 s |
| 4 | regular request | 2,543 | 422,138 | 268,375,505 | 5 GiB | 1,400 s |

The bound is `83*K*(K-1)/2+42*K`; actual usable points and root keys come
from the producer and independent replay. Lean checks K3114 as the last K
in this slot bucket and its two-array payload below 13 GiB. The runner
records sampled child RSS and stops a child above 14 GiB. It retains wall,
CPU, RSS, page faults, context switches, cgroup counters, raw output,
resource stops, and timeouts. All four arms run even if an earlier one is
censored. The overall workflow cap is 160 minutes, above the 120 minutes
of summed arm caps plus the 120-second rho cap.

The same-point strong-rho reference uses archived `fold-L64-DP6-E0`, seed
`2026101301`, one OS thread, and a 120-second cap. Its config SHA-256 is
`cf76ebb5372acfc1355b88374f875980be66cdffba1a82cb6e4c2f28bcfa723b`.

The independent verifier replays the complete point-defined base, subgroup,
four-point relation, recovered scalar, and five exclusive online phases.
It checks the state bound and table layout. At K2543 it compares dense and
regular root-key counts, base hashes, and first target witness. The N83
source test compares every canonical key and first witness at both table
densities for six representatives. Every completed arm receives a distinct
canonical `IC1` manifest and one-target workload/run ID.

One execution per arm is an exploratory Tier 1 capacity gate. A controlled
wall-time ratio requires five alternating same-point pairs and the tiered
noise gates. The branch-scoped workflow uploads all raw receipts on success
or failure and does not alter repository-wide Actions settings.
