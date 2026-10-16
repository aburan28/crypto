# N83 root11 factor-base neighborhood on 4,000 paired points

K3300 remains the selected factor base under the frozen 95% paired rule. All
four root11 bases completed and independently replayed the same 4,000 public
subgroup points under the 14 GiB RSS guard. K3250 had the lowest observed
mean, **33.546 ms** per target versus **34.239 ms** for K3300, but its paired
K3250-minus-K3300 mean was **−0.692 ms** with a contiguous-block bootstrap
interval of **[−2.109, 0.717] ms**. The upper endpoint exceeds zero, so the
predeclared selection rule retains K3300. K3250 is the next local factor-base
lead for a fresh, independently allocated confirmation panel.

## Four-arm result

| Root11 arm | Usable base points | Folded columns | Mean target online ms | Setup s | Peak RSS GiB | Replay status |
| :--- | ---: | ---: | ---: | ---: | ---: | :--- |
| K3200 | 531,200 | 3,200 | 35.466 | 576.687 | 9.762 | `PASS_INDEPENDENT_REPLAY` |
| K3250 | 539,500 | 3,250 | 33.546 | 580.051 | 9.905 | `PASS_INDEPENDENT_REPLAY` |
| K3300 | 547,800 | 3,300 | 34.239 | 598.758 | 10.052 | `PASS_INDEPENDENT_REPLAY` |
| K3350 | 556,100 | 3,350 | 34.314 | 619.775 | 10.200 | `PASS_INDEPENDENT_REPLAY` |

| Paired contrast, arm minus K3300 | Mean ms | 95% block interval ms | Confirmed faster |
| :--- | ---: | :--- | :--- |
| K3200 | 1.227 | [−0.368, 2.783] | No |
| K3250 | −0.692 | [−2.109, 0.717] | No |
| K3350 | 0.076 | [−1.432, 1.547] | No |

The comparison uses 80 contiguous blocks of 50 target points and 10,000
bootstrap resamples, with seed `2026101599+K` for each contrast. It promotes
the fastest alternative only when the upper 95% endpoint for its paired
difference from K3300 is below zero. Every arm ran in a fresh process with an
empty base, root index and log cache. Preparation time is reported separately
from the target online interval. The target interval contains query,
four-summand lookup, relation check, descent and recovery check. This is a
secondary factor-base stage comparison on one cloud instance using points
from the earlier 4,000-point confirmation panel. The selected K3300 base has
a separate fresh-target Tier 1 same-point IC/rho receipt.

The exact curve is `EC1N83Ckb1a1haacb69ee9f37`, with polynomial-basis
`GF(2^83)` modulus bitmask `9671406556917033397649557`, model
`y²+xy=x³+x²+1`, subgroup order `8569786107849059`, and cofactor
`1128547018`. The K3300 geometric base hash matches the prior confirmed
base: `a4ba33098d2776520b98a326585a2301f8ccd5bc005df7c1a70fbe8496c61d8b`.
The four canonical candidate manifests record their exact constructions and
executed source; the K3300 candidate ID in this source snapshot is
`IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0hb4769a96f33e`.
The common workload ID is `12e14f274f8c`.

## Verification and custody

The producer was frozen at commit `85ad48ca4ab668e850623416b814a350b162beef`,
with exact IC source SHA-256
`a1eb3b937106bb7cfdee6c54d6422ceb126b628462361c9ff7a127127f230c89`.
[GitHub Actions run 38670919975](https://github.com/aburan28/crypto/actions/runs/38670919975)
completed on Linux x86-64 with an exposed AMD EPYC 7763. Raw artifact
`11946288154` has SHA-256
`f4e1edb15db6e3a4306f5f208e8f0714ef1a1e81d738ed9b30d3d28d3c30649b`;
the cloud analysis has SHA-256
`fe9c87c6730befc67fa8338204e48edc35dd3c66aa1cbd21b55fa29036b8eb88`.
The verifier replayed the complete factor base and every target scalar and
four-point relation for all arms. The source-pinned sealer checks these
digests and emits [selection.json](selection.json), four canonical candidate
manifests, lossless compressed [workload](workload.json.gz) and
[measurement rows](measurements.json.gz), and their
[digest index](records_index.json). The raw ZIP is retained in the durable
N83 artifact archive.
