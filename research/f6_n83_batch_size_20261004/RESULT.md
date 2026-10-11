# n=83 F6 wide-pair residual batch size: rejected

The [preregistered protocol](PROTOCOL.md) changed only the residual query
batch in `F6WidePairIndex::solve4` from 256 to 2,048. Both arms used the
registered K0 n=83 curve, its actual cofactor-projected bases, and the same
public T001 point. All six processes exited 0; all small and full queries
returned the same exact no-witness answer with the expected pair counts.
The full-base query contains 8,219,485 unordered pairs over 4,054 usable
points. The frozen [runner](run.sh), [status rows](status.tsv), and the six
`*.jsonl` files preserve every observed run; stderr files are empty.

| Base points | Metric | Baseline 256 | Candidate 2,048 | Baseline / candidate |
| ---: | --- | ---: | ---: | ---: |
| 258 | Index build, three-run median | 19.765 ms | 21.170 ms | 0.934× |
| 258 | Exact query, three-run median | 16.936 ms | 16.703 ms | 1.014× |
| 1,048 | Index build, three-run median | 459.870 ms | 398.066 ms | 1.155× |
| 1,048 | Exact query, three-run median | 346.977 ms | 332.105 ms | 1.045× |
| 4,054 | Index build, two-run median | 12.448 s | 9.099 s | 1.368× |
| 4,054 | Exact query, two-run median | 6.981 s | 6.376 s | 1.095× |
| 4,054 | Peak process RSS, maximum | 2,678,407,168 B | 2,673,967,104 B | 1.002× |

The full runs were interleaved baseline, candidate, candidate, baseline.
Their query times were respectively 8.337, 6.782, 5.969, and 5.625 s;
their build times were 18.585, 11.271, 6.928, and 6.311 s. This broad
drift makes the unisolated timing evidence especially weak. The paired
two-run median query is **8.7% lower**, short of the registered 10% reduction
gate. The code change was reverted. The exact one-line
[candidate patch](candidate.patch) and its compiled binary hashes remain
for reproducibility; final source retains batch size 256. The patch's
SHA-256 is `d17ee94dfa9fcdfc21981fce42126e04c35cc707c905fc9ef990df44cc3090ec`.

The baseline source was commit `2a0eb9b2e3fc6c24319754b2a4d77c2e01e1202a`,
`f6_wide_geometry.rs` SHA-256
`7cae9b6c9610bfe17cb0c16e053b616971fdff3a91748e1fc6a41c12b2499a60`.
The candidate source SHA-256 was
`ce9ed373095cba098d8fd1efd4e2e987f6ae2a978898aa96011415919949da44`.
The exact small and full probe source hashes were respectively
`263a0b043c23f890a695370e66aaebbcd8ecd6a053c9413ed821201f5a2c5f39`
and `d17d2e1952d101b435295728248be061de04c4bd305db66f52dc6c890e8648ad`.
The protocol SHA-256 was
`0297c5cba04e1a1aa7a928c3ed8c66459d72501dca6413ddc00ef6c147420e06`;
the runner SHA-256 was
`1e6170ad37d8b40204c814a8af42e37c3939fd7f7eac271242fbc04a7931b9c2`.
The release full-probe binaries were SHA-256
`81aad73752ec48e391b95be66d3281834cd1799f17440cb5a3b8253eef523bc9`
(baseline) and
`8be57913e56bda6790e759a48182149409d36bcdfb2c1e274765baf025066d45`
(candidate). The small-probe binaries were SHA-256
`d41e451c46c3012506fd7160fceaee941e02447826792be8611cd8a3eae035e6`
and `c2a04a410021b0d93b5e6faf605fb08850fffec1556a2e807f4d7113afd1ba59`.
Both were built with Rust 1.93.1 on arm64 macOS (Darwin 25.6.0), release
profile and `--offline`. The host lacked an isolation receipt.

These observations are exact four-summand stage diagnostics, not a
controlled CPU speedup or an F6/F4/F5 end-to-end comparison. The dimension-12
base's uniform-target four-summand coverage ceiling is `4.662e-12`.
No ordinary relation yield, recovered discrete logarithm, one-target IC
online interval, or paired rho run was measured in this experiment.
