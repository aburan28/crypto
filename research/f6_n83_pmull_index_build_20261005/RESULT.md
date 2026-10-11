# n=83 F6 packed pair-index construction: retained

The [preregistered protocol](PROTOCOL.md) compared the retained x-only
PMULL query with a focused PMULL pair-index builder. The candidate
calculates exact full affine sums in packed degree-83 field arithmetic
while constructing `F6SignedPairIndex`. It keeps one reference inversion
per row and uses the reference group law for identity, same-x,
inverse, and repeated-point cases. The exact polynomial and ARM64
AES/PMULL feature are checked once; every other field and CPU uses the
existing portable builder. The target-query method, lookup, sign
check, and group-witness replay are unchanged.

The [geometry test log](geometry_tests.log.gz) passed all nine module
tests, including exceptional K1 full-point comparisons, selected K0
factor-base row comparisons, and exhaustive small-base witness replay.
The exact validation command was `cargo test --offline --release --lib
cryptanalysis::f6_wide_geometry::tests::` with the shared `CARGO_TARGET_DIR`
and SSD-backed `TMPDIR` used for this worktree. The native probe
binaries were compiled with `cargo build --offline --release --example`
for `f6_wide_n83_pmull_index_probe`,
`f6_wide_n83_pmull_full_probe`, and
`f6_wide_n83_pmull_planted_probe`, plus the corresponding parent
small/full examples, then copied to the frozen `target/bench_bins`
paths identified by their hashes below. The thin
orchestration command was `sh research/f6_n83_pmull_index_build_20261005/run.sh`.
A full-base K0 [planted control](planted_control.jsonl) returned
`[0,2,4,6]` from both query methods, and both sums replayed in the
curve group. This is a correctness control, **not** an ordinary-yield
observation. All six benchmark processes exited 0 with the same
exact T001 `no_witness` result. The [status table](status.tsv),
[runner](run.sh), six JSONL files, and six empty stderr files preserve
the raw observations. The [build log](candidate_build.log.gz) records
the offline release build and is losslessly compressed with `gzip -n`.

| Actual base points | Metric | Parent baseline | Packed builder | Baseline / candidate |
| ---: | --- | ---: | ---: | ---: |
| 258 | Index build, three-run median | 17.952 ms | 6.249 ms | 2.873× |
| 258 | Exact query, three-run median | 1.827 ms | 1.867 ms | 0.978× |
| 1,048 | Index build, three-run median | 327.176 ms | 104.210 ms | 3.139× |
| 1,048 | Exact query, three-run median | 43.459 ms | 53.933 ms | 0.806× |
| 4,054 | Index build, two-run median | 5.335 s | 1.733 s | 3.079× |
| 4,054 | Exact query, two-run median | 0.945 s | 0.884 s | 1.069× |
| 4,054 | Build plus query, supplementary | 6.280 s | 2.617 s | 2.400× |
| 4,054 | Maximum process RSS | 1,276,329,984 B | 1,275,772,928 B | ~1.000× |

The full-base order was baseline, candidate, candidate, baseline.
Index-build times were 5.495, 1.769, 1.696, and 5.175 s;
adjacent paired baseline/candidate build ratios were 3.105× and
3.052×. These two ratios describe observed spread, not a confidence
interval. Full-query times were 0.926, 0.933, 0.835, and 0.963 s.
Every full run used 4,054 distinct usable points, 8,219,485 unordered
pairs, 4,108,723 signed-sum representatives, and 2,027 signed columns.
The candidate full-index build median was **67.5% lower**, clearing
the registered 20% gate. The full-query median improved, peak RSS
fell slightly, and all exactness controls passed, so the guarded
builder is retained. The 1,048-point query median was 24.1% higher;
that small-base query is a diagnostic, not a gate, because its code
path is unchanged. This unisolated host cannot determine whether that
observation is repeatable. No small-base query timing is used as an
end-to-end IC claim.

Index construction is target-independent. The build-plus-query row
is a supplementary component cost and must not replace the primary
one-target online IC metric, which starts after reusable preparation.
The two full pairs and three small repeats do not supply a controlled
CPU speedup bound.

The baseline was #1394 head `8b6e6e395`, with F6 geometry source
SHA-256 `9ce12438b77d0a77fbaa414d15a2d5f59e8555a1e3dd383368eb6e22454a887c`
and field source SHA-256
`15d8f5f357e863e0f16d49f04e88bbbd3269a8ddbb930079f41dff8ef2c0f197`.
The candidate F6 and field source hashes were
`887677cbc6c91bd996caf92e79fb575d0af6295fbb1e4fb0eb4c9a7cefb2e619`
and `16ea7db9ff4792f54e3f1b21454a7a2d109b9a2e27245f93e27c1ec530b0559e`.
The protocol and runner hashes were
`76bcbca30c9e47bf19162fd6ad95c5cb7a723e7799c48da4ae421738f8ab5666`
and `7e46504b8af52009cfff72c99ad84c26adb2b888a04dd994e24b554217b83403`.
The baseline small and full binaries were SHA-256
`83b8f59d66422cfcbbe47dd840c6b476fac64cec6baea74d0c7debb6805a521c`
and `8923553ad79d2647383fd8a929b12fd4ddbfe0f880d95b9d641de295fc6e541b`.
The candidate small, full, and planted binaries were SHA-256
`1eb31927f4cfad1f8cf913d5adea846df999dfff2b215ff30c74d9826f00f20e`,
`dfca7d30489c898cfe695ef4b5d2fbfa28dcc8aa71f42e711482b9f6c624b41d`,
and `1e176fbeba3d7e1ae6cc2b23def456134c5a339613d088071282040fd9a0b94c`.
Both arms used Rust 1.93.1 release builds on physical arm64 macOS.

This host has no CPU-isolation receipt. All ratios are exploratory
four-summand stage diagnostics, not controlled complete-solver
speedups. The K0 base's uniform-target four-summand coverage ceiling
is `4.662e-12`, and T001 was an exact miss in every ordinary run.
No natural relation yield, recovered logarithm, complete F6/F4/F5
comparison, one-target IC online interval, or paired rho run was
measured. The end-to-end IC effect remains unknown.
