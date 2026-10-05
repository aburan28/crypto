# n=83 F6 sign-folded pair index: retained component

The [preregistered protocol](PROTOCOL.md) compared the existing exact
four-summand pair index with a sign-folded index on the same registered K0
n=83 curve, actual cofactor-projected bases, and public T001 point. The new
index requires a distinct base closed under point negation. It stores one
representative for each pair-sum x-coordinate, derives the negative pair
from the negated source points, and computes both target residuals with one
batch of denominator inverses. A hit is replayed in the full curve group.

All six probe processes exited 0 and produced the same exact no-witness
outcome, with expected pair counts and on-curve/subgroup checks. The
[status table](status.tsv), [runner](run.sh), six JSONL records, and six
empty stderr files preserve every run. The final-source
[correctness log](correctness.log.gz) covers the dual-sign group law and
exhaustive small-base four-sum search, including witnesses and misses.
The [geometry test log](geometry_tests.log.gz) passed all six module tests.
The two logs are compressed losslessly with `gzip -n` and can be read with
`gunzip -c`.

| Usable base points | Unordered pairs / signed sums | Index build baseline / candidate | Exact query baseline / candidate | Query ratio |
| ---: | ---: | ---: | ---: | ---: |
| 258 | 33,411 / 16,640 | 23.702 / 29.655 ms | 22.953 / 18.204 ms | 1.261× |
| 1,048 | 549,676 / 274,570 | 457.465 / 347.898 ms | 433.599 / 268.015 ms | 1.618× |
| 4,054 | 8,219,485 / 4,108,723 | 8.417 / 6.561 s | 6.705 / 5.040 s | 1.330× |

Small-base values are three-run medians within each arm. Full-base values
are two-run medians; the arms ran baseline, candidate, candidate, baseline.
The four full query times in that order were 6.860, 4.270, 5.811, and
6.550 s. The 4,054-point candidate median was **24.8% lower** than baseline,
passing the registered 20% exploratory query gate. Its full index-build
median was 22.1% lower, within the 2x build-cost limit. The small-base
query medians also improved. Maximum process RSS fell from
2,677,833,728 B (2.49 GiB) to 1,276,133,376 B (1.19 GiB), a 52.3%
reduction. The sign-folded component is retained alongside the original
index, which still supports bases without sign closure.

The baseline was commit
`2a0eb9b2e3fc6c24319754b2a4d77c2e01e1202a`;
its `f6_wide_geometry.rs` SHA-256 was
`7cae9b6c9610bfe17cb0c16e053b616971fdff3a91748e1fc6a41c12b2499a60`.
The candidate final source SHA-256 was
`9e5ce9ccc3be6a5427a6782ca5fbc65c0f4d6247037776b74969fba9a0af3ec9`.
The small and full candidate probe sources were SHA-256
`efd2bd31623b6d2ea4f3c7712eeaaf8873d74e09dc10bc4266d88ba5a6ead502`
and `d6b85dcc0be6382eac600aa974299617170d7d1a1d210a305396fa551c42bc61`.
The protocol and runner were SHA-256
`70f8e932f4fe807333ebd659573b4174a7e77fb02437e8c272c5b4c1db45c2e0`
and `7aaa4a7137f1ce4d3a760d08689d3d7c2e6f03dd3610efb408b84363e66648aa`.
The candidate full and small release binaries were SHA-256
`99ccb5c75fbf29798a261b837914ba7b004e94946e3148980eb8338e40db8b81`
and `adcddfbe351f6b4b4501f21b9e614ea8dc0169104ceecdfd9ea97ba4cc99d612`.
The baseline full and small binaries were SHA-256
`81aad73752ec48e391b95be66d3281834cd1799f17440cb5a3b8253eef523bc9`
and `d41e451c46c3012506fd7160fceaee941e02447826792be8611cd8a3eae035e6`.
Both arms used Rust 1.93.1 release builds on arm64 macOS (Darwin 25.6.0).

This is an **unisolated stage diagnostic**, not a controlled CPU speedup.
It is a four-summand meet-in-the-middle component, not a complete n=83 F6
point-decomposition or one-target IC pipeline. The dimension-12 base's
uniform-target four-summand coverage ceiling is `4.662e-12`; T001 was an
exact miss in every run. The build is target-independent and was reported
separately from the target-dependent query. Natural relation yield, a
recovered discrete logarithm, F4/F5 comparisons, and paired one-target rho
online timing remain unmeasured.
