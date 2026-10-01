# Decision: fixed degree-263 leaves break the source m10 slot-count transfer

The preregistered hypothesis was that the source curve's one-slot m10 count
cannot be copied to every slot of a **fixed** descending leaf by coordinate
squaring. The exact census confirms this on both preselected degree-263 leaves.
Each source 13-bit slot has 7,977 physical points and 3,988 nonzero `[4]`
sign classes. On leaf `[1,0]`, the ten physical counts range from 8,063 to
8,339; on leaf `[1,4]`, they range from 8,087 to 8,323. The complete ordered
slot vectors are in [ANALYSIS.json](ANALYSIS.json). The 14-bit high slot has
16,125 source points, 16,259 on `[1,0]`, and 16,323 on `[1,4]`.

| Curve | m10 arm | Physical tuple product / q | Product / source | Distinct uncompressed `[4]` sign classes across slots | Classes / source |
| --- | --- | ---: | ---: | ---: | ---: |
| Source | Balanced | 1.53294467 | 1.000000 | 39,880 | 1.000000 |
| `[1,0]` | Balanced | 1.96048190 | 1.278899 | 40,875 | 1.024950 |
| `[1,4]` | Balanced | 2.01826995 | 1.316597 | 40,994 | 1.027934 |
| Source | Unequal | 3.09875051 | 1.000000 | 43,954 | 1.000000 |
| `[1,0]` | Unequal | 3.91350217 | 1.262929 | 44,932 | 1.022251 |
| `[1,4]` | Unequal | 4.03678721 | 1.302714 | 45,075 | 1.025504 |

These are **exact enumerated representation counts**, not measured PDP hit
rates. The product/q ratio is a necessary counting ceiling with no coset,
relation-rank, solver, or target-support certificate. In these particular
domains, each nonzero liftable x has a distinct nonzero projected x within
its slot, and the ten projected slot sets are disjoint; there are no observed
projection collisions to discount. This does not create a cheap fixed-leaf
Frobenius action. The source curve can use its geometric Frobenius between
rotated slots; coordinate squaring sends either fixed leaf to a conjugate
curve. The 26–32% larger raw tuple products therefore do **not** establish
an equal-useful-size factor-base or attack advantage.

The frozen source and inputs were merged in [PR #1059](https://github.com/aburan28/crypto/pull/1059)
after the protocol in [PR #1057](https://github.com/aburan28/crypto/pull/1057)
and before any leaf scan. The source lock names commit
`1f69d4b034a50abb2be4d724454f66d0327a4f56`; the one-shot hosted
[run 36707051160](https://github.com/aburan28/crypto/actions/runs/36707051160)
used merged main `3813d1879ca03483ac3445ce6d431f0319d25f94`. Its artifact
has 33 complete scans, 294,912 raw masks and SHA-256 result
`71c6a748a3d65f61e26d139713498aa6731cf64b53c066886ed663f2489ea1d5`.
The hosted independent verifier and a separate macOS arm64/Python 3.13.1
replay each checked every raw mask, all row/chunk hashes, both source
positive controls, and 264 sampled full-point `[4]` images. Their `PASS`
receipts are byte-identical. [MANIFEST.json](MANIFEST.json) gives each
archived file's byte count and SHA-256; the raw compressed rows and all scan
summaries live in the adjacent `leaf_m10_support_20260930/evidence/` tree.
[analyze.py](analyze.py) recomputes every number above from that archive.

**Gate decision.** Reject source one-slot count transfer for fixed leaves;
retain the actual per-slot counts as inputs to any leaf-native m10 design.
This is an accounting correction, not an index-calculus advance. Next, hold
the original, transported, native, and exact-pullback arms to the same
**useful projected physical size** and same Q streams; freeze deterministic
leaf selection and the cost of transport, inverse maps, cofactor projection,
any leaf action, PDP solving, relation rank, and scalar recovery before
measurement. The review/label-gated [m10 capacity PR #937](https://github.com/aburan28/crypto/pull/937)
remains an independent prerequisite for its n131 representation run; this
census does not release that gate. Natural-target PDP yield, full ECDLP
cost, and matched-rho crossover remain unset.
