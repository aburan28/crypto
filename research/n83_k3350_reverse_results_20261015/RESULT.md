# Reversed process order changes the N83 factor-base ranking

The K3350-first replication selected K3300 on the same 20,000 public
subgroup points. K3350 took **38.314261 ms** per target and K3300 took
**33.967423 ms**. The paired K3350-minus-K3300 online mean was
**+4.346838 ms**, with a 95% contiguous 50-target block bootstrap interval
of **[+3.654838, +5.043167] ms** (10,000 resamples, seed `2026101519`).
Both arms completed all target answers and independently replayed every
scalar and four-point relation. The frozen replication rule therefore records
`REPLICATION_NOT_CONFIRMED`.

| Process order | Cloud CPU | K3300 mean online ms | K3350 mean online ms | K3350 − K3300 ms | Paired 95% interval ms |
| :--- | :--- | ---: | ---: | ---: | :--- |
| K3300, K3350 | AMD EPYC 9V45 | 28.957926 | 26.877454 | −2.080472 | [−2.600900, −1.546586] |
| **K3350, K3300** | **AMD EPYC 9V74** | **33.967423** | **38.314261** | **+4.346838** | **[+3.654838, +5.043167]** |

The two process orderings favor different bases, so the combined evidence
does not rank K3300 and K3350 by wall time under a stable host envelope.
The next comparison should balance order on one allocation, using fresh
empty-cache processes in a K3300/K3350/K3350/K3300 sequence and calibrated
noise receipts. The first K3300-first result remains a separate verified
stage record at
[`../n83_k3350_sparse_fresh_results_20261015/RESULT.md`](../n83_k3350_sparse_fresh_results_20261015/RESULT.md).

| Reversed-order arm | Usable base points | Folded columns | Setup s | Full process wall s | Peak RSS GiB | Replay |
| :--- | ---: | ---: | ---: | ---: | ---: | :--- |
| K3350 first | 556,100 | 3,350 | 521.917 | 1,338.925 | 10.201 | `PASS_INDEPENDENT_REPLAY` |
| K3300 second | 547,800 | 3,300 | 448.315 | 1,178.318 | 10.054 | `PASS_INDEPENDENT_REPLAY` |

The run used four exposed CPUs, separate empty-cache processes, a 14 GiB
peak RSS gate, and a 1,800-second wall cap per arm. Both cgroup records show
zero CPU-throttling and OOM-kill events. Setup and full-process walls are
reported separately from the online interval. The common workload ID is
`8453444a4631`; the K3350 and K3300 candidate IDs in this executed snapshot
are `IC1N83Ckb1a1fb556100PDP4rootRCguidedLAgaussTDdirectISO0h58b539381e30`
and `IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0h76d2d77dd6c3`.
Both base hashes match the first study.

The producer executed source commit
`53339d82a372e09958dd43bb7391e283de8933f5`, with the same Rust IC
source SHA-256
`c9824acebda0360416c3c507594bb07c17ca798c32cbf8f84bba95ca2c44d5b1`.
[Actions run 38705055043](https://github.com/aburan28/crypto/actions/runs/38705055043)
passed source preflight, both capped arms, full independent replay, and raw
upload. Artifact `11960333030` has ZIP SHA-256
`8ad1eb8759a7d6d9872b93697990992e0a91eb2b728bdfd211d4a3c744ac8c76`;
its `analysis.json` SHA-256 is
`a7d94c1a201868f8084d414edd2a1d99baf853ac92393cf4fe616c6ba8732db7`.
The raw ZIP and a separate local verifier receipt are archived on SSD2 at
`/Volumes/SSD990-2/n83-archive-20261015/raw/`. The local verifier replayed
both bases and all 40,000 answers and reproduced the analysis digest byte
for byte. `records_index.json` pins the compressed workload and measurements;
the current-head CI workflow regenerates both and requires a clean diff.
