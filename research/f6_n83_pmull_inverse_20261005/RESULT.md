# n=83 F6 packed inversion: rejected

The [preregistered protocol](PROTOCOL.md) tested replacing the last
portable field inversion in the packed degree-83 F6 pair-index builder
and x-only query with a PMULL Itoh–Tsujii inverse. The
[rejected patch](REJECTED.patch) preserves the exact implementation;
runtime source was reverted after it missed the registered 15% full-base
query reduction gate.

The [release geometry test log](geometry_tests.log.gz) passed all ten
module tests. The new test compared zero and 10,000 deterministic
nonzero inverses with the reference implementation and independently
checked `a * inverse(a) = 1`. Existing tests covered ordinary and
exceptional point sums and exhaustive small-base four-sum search. The
portable fallback source was unchanged; no physical x86 run was made.
The full-base K0 [planted control](planted_control.jsonl)
returned `[0,2,4,6]` under both reference and packed queries, and
both witnesses replayed in the curve group. This is a correctness
control, not ordinary relation yield. All six benchmark processes
exited 0, and all ordinary T001 queries reported exact `no_witness`.
The [status table](status.tsv), [runner](run.sh), six JSONL outputs,
and empty stderr files preserve failures as well as successes.

| Actual base points | Metric | #1399 baseline | Packed inversion | Baseline / candidate |
| ---: | --- | ---: | ---: | ---: |
| 258 | Index build, three-run median | 6.800 ms | 5.501 ms | 1.236× |
| 258 | Exact query, three-run median | 2.241 ms | 1.589 ms | 1.410× |
| 1,048 | Index build, three-run median | 71.427 ms | 65.984 ms | 1.082× |
| 1,048 | Exact query, three-run median | 33.182 ms | 28.716 ms | 1.156× |
| 4,054 | Index build, two-run median | 1.423 s | 1.436 s | 0.991× |
| 4,054 | Exact query, two-run median | 0.675 s | 0.664 s | 1.017× |
| 4,054 | Maximum process RSS | 1,276,231,680 B | 1,276,264,448 B | ~1.000× |

The full-base order was baseline, candidate, candidate, baseline.
Query times were 0.627, 0.659, 0.669, and 0.724 s; adjacent paired
baseline/candidate ratios were 0.951× and 1.082×. The median query
reduction was only **1.66%**, below the registered 15% gate and within
this unisolated host's observed variation. Index-build median rose by
0.95%; peak RSS was essentially flat. Every full run used 4,054
distinct usable points, 8,219,485 unordered pairs, 4,108,723
signed-sum representatives, and 2,027 signed columns. The smaller
query improvements are diagnostics and do not override the frozen
full-base gate.

The exact commands were `cargo test --offline --release --lib
cryptanalysis::f6_wide_geometry::tests::`, `cargo build --offline
--release --example f6_wide_n83_pmull_index_probe --example
f6_wide_n83_pmull_full_probe --example f6_wide_n83_pmull_planted_probe`,
and `sh research/f6_n83_pmull_inverse_20261005/run.sh`, with the same
SSD-backed `TMPDIR` and shared `CARGO_TARGET_DIR` in both arms. The
[candidate build log](candidate_build.log.gz) and test log use lossless
`gzip -n` compression. Both arms used Rust 1.93.1 release builds on
physical arm64 macOS. The #1399 baseline F6 source SHA-256 was
`887677cbc6c91bd996caf92e79fb575d0af6295fbb1e4fb0eb4c9a7cefb2e619`;
the candidate source was
`26a9472523b95147a0c4b466bc42da563beda607ed8981a3f1544fc89d9d6be3`.
The protocol and runner hashes were
`fffd4bf936e915204fa793761b0d335093d7f6ac5521e0017d25e49461694af3`
and `2acb0e76da4bba38074f834f2f76d2e7193dae04cfe231770c266cc4d5f91894`.
The baseline small/full binary SHA-256 values were
`1eb31927f4cfad1f8cf913d5adea846df999dfff2b215ff30c74d9826f00f20e`
and `dfca7d30489c898cfe695ef4b5d2fbfa28dcc8aa71f42e711482b9f6c624b41d`;
candidate small/full/planted values were
`19b51f4d7b6f3c551beb076afc2a7c716b6f732594ffd18de51ae995ba3fbfb9`,
`50a798d42cfe9bcecc304bf9fec89bf7e64c52b5126e573aedbb15619aa1091c`,
and `47dcf5956762f82a06e3b1ce18439182241d063cb6dc99b11ab8c0ba8186e39a`.

These unisolated timings are exploratory four-summand stage diagnostics.
Index construction is target-independent; none of these costs is a
complete one-target IC online interval. The K0 base's uniform-target
four-summand coverage ceiling is `4.662e-12`, and all ordinary queries
were exact misses. No natural relation yield, recovered logarithm,
complete F6/F4/F5 comparison, or paired rho run was measured. The
end-to-end IC effect remains unknown.
