# Fresh 4,000-point comparison of N83 K3300 and K3100 factor bases

The preceding 1,000-point refinement selected K3300/root11, with a paired
runner-up-minus-selected mean of 1.875 ms and a target-block interval that
included zero. This confirmation freezes 4,000 new public subgroup points
before either arm runs and compares the exact selected K3300/root11 base
with the exact K3100/root22 runner-up on one cloud allocation. Their
refinement candidate IDs are
`IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0h0f8ab8bfcff7`
and
`IC1N83Ckb1a1fb514600PDP4rootRCguidedLAgaussTDdirectISO0haa5b1411e6a0`.
The source base digests must reproduce
`a4ba33098d2776520b98a326585a2301f8ccd5bc005df7c1a70fbe8496c61d8b`
and
`bea9c9cda2826a16c5e7309f2fd5267085190e91978a535eba34d254b00554e6`.

Targets use the deterministic input law
`scalar_i = 1 + (little_endian(SHA256("n83-factor-base-confirm-4000-20261015-v1" || u32_le(i))) mod (r-1))`
for `i = 0,...,3999` in the prime subgroup of order
`r = 8569786107849059` on `EC1N83Ckb1a1haacb69ee9f37`. The point
file contains only public coordinates; a separate scalar sidecar supports
independent replay. Generation rejects duplicates against the prior N83
selection and fresh-validation fixtures. Exact file digests are pinned in
`generate.py`, `run.py`, and `verify.py`: public targets SHA-256
`451df8524b2c356dd4a523260620c4dff406975a122762b98ddddb7987b75c2c`,
independent scalar sidecar SHA-256
`8e18c9c97ba9f57e5eb43f1b39b9e9575b5cfab87840c25c428dac150c868101`.

| Order | Arm | Usable points | State/table layout | Child cap |
| ---: | --- | ---: | --- | ---: |
| 1 | K3300/root11 | 547,800 | 11-byte states, exact five-byte `2^30` slots | 1,200 s |
| 2 | K3100/root22 | 514,600 | 22-byte states, exact five-byte load-75 slots | 1,200 s |

Each arm uses a fresh process and empty factor-log/index setup, rank seed
29, symmetric S3 traversal, four rank threads, one target thread, and a
14-GiB observed RSS guard. The independent verifier checks every base
representative, its 166-point orbit, subgroup membership, every recovered
scalar and four-point relation, all exclusive target-online phase sums,
the two frozen base digests, and all 4,000 same public points. Setup, final
rank, final matrix solve, memory, process wall, attempts, and incomplete
arms stay in the raw receipt.

The predeclared contrast is
`Δ = mean_i(online_ms(K3100, i) − online_ms(K3300, i))`.
The 95% interval resamples 80 contiguous 50-target block means 10,000
times with seed `2026101599`; the lower endpoint above zero confirms
K3300's online advantage under this workload and allocation. An upper
endpoint below zero confirms K3100's advantage; any interval crossing zero
is inconclusive. This is a secondary paired factor-base confirmation. The
five-point same-target IC/rho study in PR #1855 separately measures the
single-target Tier 1 online result.
