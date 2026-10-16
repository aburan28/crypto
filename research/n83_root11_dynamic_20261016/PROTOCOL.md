# N83 dynamic one-root exact index, one-point cold screen

The one-root state stores one 83-bit S3 root in 11 bytes and reconstructs its
second root when an exact table collision or lookup needs it. The existing
implementation inserts both exact roots during state construction but
allocates a fixed `2^30` five-byte-slot table. This experiment sizes that
early table from the known symmetric state count. It tests a potentially
smaller cold index without changing the factor base, relation law, or
first-inserted witness semantics.

Use the single previously unseen public point frozen in
[`n83_cold_factor_base_20261016`](../n83_cold_factor_base_20261016/PROTOCOL.md),
SHA-256 `2c068d2f659588edcdbaa685fc73bb146e2bfee3be123d377d537658f8021c6d`.
The verifier-only scalar sidecar has SHA-256
`b71c5a8bfdae603931fa6c03fc6f6faa3b261288cf787ac682bf8a3e3afe262d`.
The curve is `EC1N83Ckb1a1haacb69ee9f37`, subgroup order
`8569786107849059`, model `y²+xy=x³+x²+1` over polynomial-basis `GF(2^83)`
with modulus bitmask `9671406556917033397649557`. Rank seed is 29.
The executed source SHA-256 is
`d094345474605075566d646ddecef8f694cdd644b0f9190a6b1ef1551fd06888`.

| Order | Arm | Usable base points | State and table | Child cap |
| ---: | --- | ---: | --- | ---: |
| 1 | K1271 fixed A | 210,986 | 11-byte; `2^30` five-byte slots | 900 s |
| 2 | K1271 dynamic A | 210,986 | 11-byte; `2^28` five-byte slots | 900 s |
| 3 | K1271 root22 | 210,986 | 22-byte; dynamic five-byte slots | 900 s |
| 4 | K1024 dynamic | 169,984 | 11-byte; dynamic five-byte slots | 900 s |
| 5 | K1536 dynamic | 254,976 | 11-byte; dynamic five-byte slots | 900 s |
| 6 | K1271 dynamic B | 210,986 | repeated dynamic control | 900 s |
| 7 | K1271 fixed B | 210,986 | repeated fixed control | 900 s |

All bases scan polynomial x from one and retain the first K signed-Frobenius
orbits after subgroup projection: exactly `166*K` distinct usable points and
K effective columns. Symmetric S3 retains
`S = 83*K*(K-1)/2 + 42*K` states. The dynamic table requests
`ceil(4*S/3)` entries; its constructor doubles the request and rounds the
slot count to a power of two. At K1271, `S=67,041,437`, so the request is
89,388,583 and the table allocates `2^28` five-byte slots, 1.25 GiB. The
fixed control allocates 5 GiB. Both retain exact root equality on lookup;
Bloom is off. Four rank threads, one target thread, sparse canonicalization,
and the target point are identical in all arms.

Run the seven arms in order on one Ubuntu 24.04 allocation with at least four
visible CPUs and 15.5 GiB host memory. Enforce observed peak RSS at most
15.25 GiB, zero new cgroup throttle, OOM and OOM-kill events, and at least
2 GiB free disk at admission. Preserve exact source and binary hashes,
topology, affinity, cgroup limits, per-arm wall/CPU and fault receipts, and
every timeout or failed producer. Run the frozen strong-rho binary once on
the same public point after the IC arms, seed `2026101601`, cap 120 s.

Independently replay every complete base, four-point relation, recovered IC
scalar, and rho scalar. For each IC arm, record the producer's cold
`process_total`, exclusive setup/rank/target stages, and the five exclusive
target-online phases. The outer child wall and replay time remain separate.
The fixed and dynamic K1271 A/B cold-process differences must each be at
most 10% of their means before descriptive comparisons are admitted.
The selected factor-base ID remains null until independent confirmation on
the four already frozen reserve points. This cloud-instance screen does not
replace alternating same-point Tier 1 IC/rho pairs for a speedup claim.
