# Shared degree-83 PMULL field operations: rejected by small-query gate

The [preregistered protocol](PROTOCOL.md) tested a guarded ARM64 PMULL
fast path in shared `F2mElement::mul`, `square`, and `square_k_times` for
the exact polynomial `z^83 + z^45 + z^2 + z + 1`. Every other field and
CPU retained the portable path. The [rejected patch](rejected_candidate.patch.gz)
preserves the complete tested source and independent field test; use
`gzip -dc rejected_candidate.patch.gz | git apply -` against the parent
branch to restore it. The
final branch restores the baseline runtime source.

The [independent field test](field_tests.log.gz) passed edge values,
10,000 deterministic pseudorandom multiply/square comparisons against
schoolbook arithmetic, repeated squares, and inversion round trips.
The [general field suite](general_field_tests.log.gz) passed 17 tests;
the [F6 geometry suite](geometry_tests.log.gz) passed eight. A full-base
K0 [planted control](planted_control.jsonl) returned and group-replayed
`[0,2,4,6]` through both F6 query paths. Planted controls do **not**
estimate natural relation yield. All six paired benchmark processes
exited 0 with the same exact T001 `no_witness` outcome. The
[status table](status.tsv), [runner](run.sh), six JSONL files, and six
empty stderr files preserve the raw observations. The
[candidate build log](candidate_build.log.gz) records the offline
release build; logs are losslessly compressed with `gzip -n`.

| Actual base points | Metric | Parent baseline | Shared-field candidate | Baseline / candidate |
| ---: | --- | ---: | ---: | ---: |
| 258 | Index build, three-run median | 22.300 ms | 13.288 ms | 1.678× |
| 258 | Exact query, three-run median | 1.911 ms | 2.131 ms | 0.897× |
| 1,048 | Index build, three-run median | 361.491 ms | 210.759 ms | 1.715× |
| 1,048 | Exact query, three-run median | 57.118 ms | 45.666 ms | 1.251× |
| 4,054 | Index build, two-run median | 7.784 s | 3.858 s | 2.018× |
| 4,054 | Exact query, two-run median | 2.021 s | 1.144 s | 1.767× |
| 4,054 | Build plus query, supplementary | 9.805 s | 5.002 s | 1.960× |
| 4,054 | Maximum process RSS | 1,276,149,760 B | 1,276,133,376 B | ~1.000× |

The full-base order was baseline, candidate, candidate, baseline.
Index-build times were 8.763, 4.038, 3.678, and 6.805 s; query times
were 2.276, 1.054, 1.234, and 1.767 s. Every full run used 4,054
distinct usable points, 8,219,485 unordered pairs, 4,108,723 signed-sum
representatives, and 2,027 signed columns. The full index-build median
fell **50.4%**, clearing that registered 20% threshold; full query and
1,048-point query also improved, and RSS passed. The 258-point query
median rose **11.49%** (1.911 to 2.131 ms), missing the registered
maximum regression of 10%. Candidate small-query observations ranged
from 1.687 to 2.867 ms while baseline observations ranged from 1.894
to 1.920 ms. The tiny absolute difference and high candidate spread
make a genuine small-query slowdown uncertain, but the preregistered
retention gate failed. The shared runtime fast path was reverted.

Index construction is target-independent. The build-plus-query row is
a supplementary component cost and must not replace the one-target
online IC metric, which charges only target-dependent work after
reusable preparation. The two full pairs and three small repeats do
not supply a confidence interval or host-isolation proof.

The baseline was #1394 head `8b6e6e395`, with field source SHA-256
`15d8f5f357e863e0f16d49f04e88bbbd3269a8ddbb930079f41dff8ef2c0f197`.
The tested candidate commit was `36221a6ad`, with field source
SHA-256 `3302c9e93dddb09666243bf8fca01ba05c873e6ef5731ce6e15173efbb151f62`.
The uncompressed rejected patch SHA-256 is
`574911e616c1ad0e4d405ddcfd5b19227010723c2631530f75cf16e603f72836`.
The protocol and runner were SHA-256
`695257702f43af067f5ad8bb568050a4ae6da129113c5cb6ccfd546516adfa88`
and `d0a46776ccdc72e1c8400f486a1c27097371a4085759734aeddc0d81b015e5b7`.
The parent small and full probe binaries were SHA-256
`83b8f59d66422cfcbbe47dd840c6b476fac64cec6baea74d0c7debb6805a521c`
and `8923553ad79d2647383fd8a929b12fd4ddbfe0f880d95b9d641de295fc6e541b`.
The candidate small, full, and planted probe binaries were SHA-256
`d6384c1a0831991e31ddf90ed3171ba54cdf155711263241bf4e56ca8f31e933`,
`d13ac4ec6493ce54a44eb551ffd6238c6c83aed5f75368d57448761d638ed9b5`,
and `6f7b9628b39b3551aa2b41317b69c625fa52701379855a2128703cebdcdef45f`.
Both arms used Rust 1.93.1 release builds on physical arm64 macOS.

This host has no CPU-isolation receipt, so every timing ratio above is
an exploratory stage diagnostic. The four-summand K0 base's
uniform-target relation coverage ceiling is `4.662e-12`, and T001 was
an exact miss in every ordinary run. No natural relation yield,
recovered logarithm, complete F6/F4/F5 comparison, one-target IC online
interval, or paired rho solve was measured. The end-to-end effect
remains unknown.
