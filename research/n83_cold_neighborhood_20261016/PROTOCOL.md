# N83 cold neighborhood around K1024

This screen tests whether a nearby contiguous-orbit factor base improves the
complete cold one-point process interval. The public target is the already
exposed point in `research/n83_cold_factor_base_20261016/targets.jsonl`, SHA-256
`2c068d2f659588edcdbaa685fc73bb146e2bfee3be123d377d537658f8021c6d`.
The four pre-frozen reserve points remain untouched. The producer receives
only the public point; the answer sidecar is used by the independent verifier.

The curve is `EC1N83Ckb1a1haacb69ee9f37`, polynomial-basis `GF(2^83)`
with modulus bitmask `9671406556917033397649557`, model
`y²+xy=x³+x²+1`, subgroup order `8569786107849059`, cofactor
`1128547018`, and rank seed 29. The IC source is the exact dynamic-index
revision with SHA-256
`d094345474605075566d646ddecef8f694cdd644b0f9190a6b1ef1551fd06888`.
Each K gives exactly `166*K` subgroup-usable factor-base points before orbit
folding and K relation columns. Every arm starts with an empty base, index,
and factor-log system. Root11 uses 11-byte states and dynamic exact five-byte
slots; root22 uses 22-byte states and dynamic exact five-byte slots. Neither
uses a Bloom filter. All arms use symmetric S3, sparse canonicalization,
four rank threads and one target thread.

| Order | Arm | K | Usable base points | Layout | Child cap |
| ---: | --- | ---: | ---: | --- | ---: |
| 1 | k1024dynamic_a | 1024 | 169,984 | root11 dynamic | 900 s |
| 2 | k896dynamic | 896 | 148,736 | root11 dynamic | 900 s |
| 3 | k1152dynamic | 1152 | 191,232 | root11 dynamic | 900 s |
| 4 | k1271dynamic | 1271 | 210,986 | root11 dynamic | 900 s |
| 5 | k1024root22 | 1024 | 169,984 | root22 dynamic | 900 s |
| 6 | k1152root22 | 1152 | 191,232 | root22 dynamic | 900 s |
| 7 | k1024dynamic_b | 1024 | 169,984 | repeated base/index A/A | 900 s |

Run all arms in this order on one Ubuntu 24.04 allocation exposing at least
four CPUs and 15.5 GiB memory, with a 15.25 GiB RSS guard and at least 2 GiB
free disk at admission. Record source and binary hashes, visible affinity,
cgroup limits, major faults, wall/CPU time, throttling and OOM deltas. After
the IC arms, run strong rho once on the same point with seed `2026101601`
and a 120-second cap. Independently replay every completed factor base,
four-point relation, IC scalar and rho scalar. Preserve timeout and failed
rows with their resource receipts.

The cold clock is producer `process_total`, from initial setup through its
one-target verification. Retain exclusive setup, rank, target and online
phase costs, plus the separately measured process outer wall and independent
replay. The A/B K1024 dynamic cold-process difference must be at most 10% of
their mean, with matching base hash and exact index entries, before ranking
resource-admitted completed arms. A capped arm cannot enter the ranking.
Advance the two lowest complete pipelines to the already frozen reserve-point
comparison, alongside the prior K1024/root22 and K1536/root22 finalists.
This exposed-point neighborhood is an exploratory cold screen. The selected
factor-base ID remains null until reserve-point confirmation; an online
IC/rho speedup requires a separate Tier 1 panel of at least five fresh
same-point alternating pairs under `docs/IC_BENCHMARK_TIERS.md` in the
cryptanalysis repository.
