# N83 sparse K3250 versus K3300 on 20,000 new public points

The preceding root11 neighborhood screen found K3250's mean target online
time 0.692 ms below K3300's on 4,000 paired points, with a 95% block
interval [−2.109, 0.717] ms. This study resolves that close factor-base
choice on 20,000 new subgroup points under the sparse canonical selector,
which is fixed **on in both arms**. Bloom is fixed off. The exact curve,
source kernel, point order, five-byte `2^30` root table, 11-byte states,
four-summand decomposition and resource envelope are common to both arms.

The public points derive from SHA-256 domain
`n83-k3250-sparse-fresh-20000-20261015-v1`; the generator checks
disjointness from earlier factor-base, refinement, reserve and confirmation
panels, including the five points frozen for Bloom Tier 1. It verifies the
subgroup generator's order once, computes each public point as `scalar*G`,
checks the curve equation and freezes the scalar sidecar separately. The
independent verifier checks every frozen scalar and every recovered
four-point target relation for each arm. Only public points enter the IC
producer. The hashes below are filled from the frozen generator output
before cloud execution:

| Frozen input | SHA-256 |
| :--- | :--- |
| Public points | `c80f5217350d3b89bbf252298295b66e81dea673c222ce36816c77c182e9cda7` |
| Answer sidecar | `9cfe15c36e06a42592344db6781a47c828e55bd0bfcb31239019ec4f6977f045` |
| Executed IC source | `c9824acebda0360416c3c507594bb07c17ca798c32cbf8f84bba95ca2c44d5b1` |

The exact point and answer files are retained in executed source commit
`402c3142fb3f1cfe70b9de58d54e51e50b2a40c5` and in the archived
producer artifact. To recreate the byte-identical files from this branch,
run `python3 research/n83_k3250_sparse_fresh_20261015/generate.py` from the
repository root. The generator checks the prior panels and prints the two
hashes above. Generated files are ignored in this branch to keep the review
diff focused on the method and results; the source-pinned replay uses the
executed commit's original fixtures.

The K3250 base has 539,500 usable subgroup points, 3,250 folded columns,
base hash `52ec996e34e9c2819cfb715898c9f605fcb896dd528c5dac01b6994b85900776`,
and 876,670,157 retained roots. K3300 has 547,800 usable points, 3,300
folded columns, base hash
`a4ba33098d2776520b98a326585a2301f8ccd5bc005df7c1a70fbe8496c61d8b`,
and 903,852,075 retained roots. Each arm starts a fresh process with an
empty base/index/log cache, a 14 GiB observed RSS guard and an 1,800-second
child cap. Target-independent preparation and complete process cost are
recorded separately from target-online query, decomposition, relation
check, descent and recovery check.

The frozen statistic is K3250 minus K3300 paired target online time over
400 contiguous 50-target blocks, 10,000 bootstrap resamples and seed
`2026101513`. K3250 is promoted only if the upper 95% endpoint is below
zero after both complete bases and all 40,000 target records pass independent
replay under the resource guard. Capped, failed, noisy and incomplete arms
remain in the raw artifact; their aggregate selection stays incomplete.
The study is a factor-base stage comparison on one cloud instance. Any new
selection requires fresh Tier 1 same-point IC/rho validation.
