# N83 symmetric S3 index: 422,138-point factor base

The symmetric-S3 state quotient built a **422,138-point N83 factor base**
(`K=2543`) within 12.30 GiB peak child RSS and recovered the frozen public
target scalar in **12.075402 ms online**. The same-point strong-rho run took
1,880.391384 ms online on the same cloud instance, a 155.72 ratio for this
exploratory execution. Complete base, subgroup, four-point relation, scalar,
and phase-sum checks passed both on the runner and in an independent local
replay. The local replay reproduced the runner's `analysis.json` SHA-256
`a379888bb02f501b1debf331bfba5f928ef814154eaa47fd5967544f3e9ed5e0`.

The quotient keeps `(left,right,relative)` only for `left<right`, or for
`left=right` with `relative≤41`. S3 exchange and Frobenius transport are the
mathematical basis. K1798 used 134,162,265 retained states versus
268,322,732 full states, and its 2.5-GiB exact root table used the same
268,316,606 entry count as the full-state table. Its peak RSS fell from
12.20 GiB to 6.21 GiB. The original source test checked key and first-witness
equality for six N83 representatives; the large K1798 receipt checks equal
key **counts** and base hashes, with the scope correction in [AMENDMENT.md](AMENDMENT.md).

## One-point measured rows

All rows use curve `EC1N83Ckb1a1haacb69ee9f37`, workload
`64d976711e0e`, the same public target, four visible CPUs, seed 29,
one target thread, and a fresh process with empty base/index/log tables.
The factor-base column is the actual count before sign/Frobenius folding.
Every listed run has status `PASS_INDEPENDENT_REPLAY`.

| Index | K | Factor-base points | Exact root keys | Target probes | Online IC ms | Setup ms | Peak RSS GiB | Rho/IC online |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Symmetric S3 | 2,543 | 422,138 | 536,737,856 | 16,174 | 12.075402 | 430,337.169 | 12.30 | 155.72 |
| Symmetric S3 | 2,200 | 365,200 | 401,711,801 | 119,402 | 62.602255 | 276,589.412 | 10.49 | 30.04 |
| Symmetric S3 | 1,798 | 298,468 | 268,316,606 | 119,686 | 63.461915 | 184,520.302 | 6.21 | 29.63 |
| Full states | 1,798 | 298,468 | 268,316,606 | 111,004 | 61.057458 | 440,056.218 | 12.20 | 30.80 |

The K2543 online interval comprises 0.000030 ms target query generation,
11.920432 ms point decomposition, 0.022351 ms relation check, 0.000361 ms
target descent, and 0.131427 ms scalar recovery check. Setup and rank are
recorded separately in [measurements.json](measurements.json). The
K2543 candidate is
`IC1N83Ckb1a1fb422138PDP4rootRCguidedLAgaussTDdirectISO0hba518cf0de57`;
its full immutable configuration is [candidate_K2543_symmetric.json](candidate_K2543_symmetric.json).
The other three manifests retain their distinct bases and state policies.

## Hardware and receipt boundary

Run [38501277471](https://github.com/aburan28/crypto/actions/runs/38501277471)
used a GitHub Actions Ubuntu 24.04 runner with AMD EPYC 7763, four visible
logical CPUs, 16,373,452 KiB reported memory, and image
`20261004.327.1`. The executed commit was
`05dcac9d1ced7e33d4a3f5bb7cb5871202ca06db`; source SHA-256 was
`e3a37208fc72fe809a8bbc2a1604f0d59f8a6582a9deaa3fd58a077083bf2c3f`.
All four children exited zero within their caps; cgroup memory OOM and CPU
throttling counters stayed at zero. The child involuntary context-switch
counts were 38,808 to 40,566, so the one-run timing ratios remain Tier 1
exploratory observations. The exact resource and failure fields are in the
artifact's `usage.json` files.

The raw [artifact ZIP](raw-artifact-38501277471.zip) has SHA-256
`4012697a128029113b7795325ec4c1ade709f0d837d9e91c3ac735d5c3eeea83`,
matching GitHub's artifact digest. `materialize.py` binds its independently
replayed rows to the four canonical `IC1` manifests and the one-point
[workload.json](workload.json). The next online comparison freezes a newly
generated point after base selection, then alternates five candidate/rho
pairs on the same resource envelope and applies the Tier 1 noise gates.
