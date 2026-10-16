# K3350 improves the N83 sparse root11 online stage by 7.18%

> **Superseded ranking (annotation added 2026-10-16, claim audit; the frozen
> record below is unchanged).** The 7.18% K3350 win did not replicate. With
> the process order reversed (K3350 first), K3350 was **slower** by
> +4.346838 ms, 95% interval [+3.654838, +5.043167]
> (`REPLICATION_NOT_CONFIRMED`,
> [`../n83_k3350_reverse_results_20261015/RESULT.md`](../n83_k3350_reverse_results_20261015/RESULT.md)).
> A balanced K3300/K3350/K3350/K3300 run on one allocation gave −0.429354 ms,
> interval [−2.057542, +1.227720] (`BOUNDED_UNKNOWN`,
> [`../n83_k_abba_results_20261015/RESULT.md`](../n83_k_abba_results_20261015/RESULT.md)),
> and K3300's first-versus-second-process shift alone (+3.39 ms) exceeded
> the balanced difference. The 7.18% is a process-order effect within this
> run, not a factor-base improvement; K3300 and K3350 are unranked. Read the
> title and "improves" below as the claim made on 2026-10-15, now withdrawn.

K3350 is the selected factor base in this fresh sparse root11 comparison.
Both arms solved and independently replayed all 20,000 new public subgroup
points and their four-summand relations. The paired K3350-minus-K3300 mean
online time was **−2.080472 ms** per point, with a 95% contiguous 50-point
block bootstrap interval of **[−2.600900, −1.546586] ms** (10,000 resamples,
seed `2026101517`). Its upper endpoint is below zero, satisfying the frozen
promotion rule.

| Root11 arm | Usable base points | Folded columns | Mean target online ms | Setup s | Peak RSS GiB | Replay |
| :--- | ---: | ---: | ---: | ---: | ---: | :--- |
| K3300 | 547,800 | 3,300 | 28.957926 | 424.495 | 10.055 | `PASS_INDEPENDENT_REPLAY` |
| **K3350** | **556,100** | **3,350** | **26.877454** | **394.673** | **10.200** | `PASS_INDEPENDENT_REPLAY` |

The observed online reduction is **7.18%** against K3300 on this panel.
Each arm ran in a separate empty-cache process, under a 1,800-second wall
cap and a 14 GiB peak RSS gate, on four exposed CPUs of a GitHub Actions
AMD EPYC 9V45 Linux instance. Producer process walls were 1,050.242 s and
976.680 s, respectively. Setup is recorded separately from the online
target interval. The online interval contains query, four-summand lookup,
relation check, descent, and recovery check. This is a cloud-instance
factor-base stage comparison; the existing Tier 1 same-point IC/rho receipt
uses K3300. The next end-to-end measurement is K3350 with the selected
sparse plus Bloom kernel on unseen same-point targets, paired to rho.

The curve is `EC1N83Ckb1a1haacb69ee9f37`, with polynomial-basis
`GF(2^83)` modulus bitmask `9671406556917033397649557`, model
`y²+xy=x³+x²+1`, subgroup order `8569786107849059`, and cofactor
`1128547018`. K3350 has exact base hash
`dc2e0c5185b88b2b0c74ede0bbc7b889d756f2f21558e9b932850d97a1bdd07b`.
Its canonical candidate is
`IC1N83Ckb1a1fb556100PDP4rootRCguidedLAgaussTDdirectISO0h05cb9f088227`;
the paired workload ID is `8453444a4631`. The candidate manifest records
the exact base construction, source, flags, and stage policies.

The producer executed commit `01a93f468067bc5d7599c089a0e0cc53293b3b4b`
with Rust IC source SHA-256
`c9824acebda0360416c3c507594bb07c17ca798c32cbf8f84bba95ca2c44d5b1`.
[Actions run 38697369233](https://github.com/aburan28/crypto/actions/runs/38697369233)
passed the exact-source build, both capped arms, independent replay of all
40,000 answers and four-point relations, and raw upload. Artifact
`11956213859` has ZIP SHA-256
`7b83f5c57764e265ee26cfab2f76e9dffd79eabd6c9f65bec66afd4dbe5c4897`;
its `analysis.json` SHA-256 is
`0ce80ea1a49564ce0e8bc4d797f4b7d32fef763e7a702c640e5d0308694c9006`.
The ZIP is archived at
`/Volumes/SSD990-2/n83-archive-20261015/raw/n83-k3350-sparse-fresh-38697369233.zip`.
A separate local verifier replay passed and reproduced the analysis digest
byte for byte. `records_index.json` pins the compressed workload and both
measurement rows. The current-head workflow independently regenerates and
replays the executed snapshot and checks canonical record equality.
