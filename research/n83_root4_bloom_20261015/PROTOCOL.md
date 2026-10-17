# N83 root4 with a 1 GiB Bloom filter

This panel tests whether a root4 factor base below the current 29-bit state
index ceiling is faster when the exact root table is preceded by a Bloom
filter. It uses one fixed executable source snapshot and compares the full
target-dependent online path against K3300/root11/wide3.

| Pipeline | Usable base points | Folded columns | Exact root index | Bloom |
| :--- | ---: | ---: | :--- | :--- |
| `root4k3350std` | 556,100 | 3,350 | four-byte, 2^31 slots | standard: 1 GiB, three hashes |
| `root4k3400std` | 564,400 | 3,400 | four-byte, 2^31 slots | standard: 1 GiB, three hashes |
| `root4k3425std` | 568,550 | 3,425 | four-byte, 2^31 slots | standard: 1 GiB, three hashes |
| `root11k3300wide3` | 547,800 | 3,300 | five-byte, 2^30 slots | wide3: 1 GiB, three hashes |

All four candidates receive the same **14.5 GiB observed peak RSS guard**
and 1,200-second child wall cap. This is a different resource envelope from
the preceding 14 GiB panels; no cross-panel wall-time or selection ratio is
formed. A host must expose at least 256 MiB above the guard, 4 CPUs, and
2 GiB free disk before execution. An arm beyond the guard or terminated by
the host remains an explicit bounded result. Each process starts with empty
reusable base, index, and logs. The fixed order is **A B C D D C B A** on
one cloud allocation.

The common 1,000-point panel uses SHA-256 domain
`n83-root4-bloom-fresh-1000-20261015-v1` and excludes earlier factor-base,
Bloom, ABBA, high-K, six-arm, and root22/root21 challenger panels. Known
scalars are a separate sidecar outside producer inputs and timing. The
targets SHA-256 is
`71b93fd3a397ee2c1204dc14fc451894eee7537b99f50d75c98d65dee6b1e3d1`;
the scalar sidecar SHA-256 is
`5417c1a10a4c2240eaec304627eddfbd287e1d005f674b1a2b17b3b957aa53e5`.

The verifier replays every base orbit, target scalar, and four-point
relation, checks exact index layout and source receipts, and requires each
pipeline's repeated base digest and root-key count to agree. Eligible arms
must have zero new cgroup throttled periods and OOM kills. A completed panel
averages the two online intervals for each point and pipeline, then computes
ordered pairwise 99% intervals from 40 contiguous blocks of 25 targets and
10,000 resamples with seeds `2026101580+i*4+j`. Select the lowest mean only
when its upper paired endpoint is below zero against each alternative.
Keep first-minus-second contrasts with seeds `2026101590+i`; otherwise record
`BOUNDED_UNKNOWN`. This is a secondary pipeline screen. A selected root4
configuration advances to fresh same-point Tier 1 IC/rho validation under
the same 14.5 GiB envelope.
