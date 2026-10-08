# n83 F6 compact signed-pair storage: retained component

The [preregistered gate](PROTOCOL.md) passes on the pinned K0 curve
`icv1-f2m83-tm6151469093347-debefd74`. The signed four-summand index now
stores each pair-sum representative as packed field coordinates and four
`u32` point indices. The n83 ARM64 PMULL residual query reads those words
directly; nonmatching curves and CPUs retain the portable group-law path.
An index with more points than `u32` can address is rejected. A matching x
key still triggers exact sign selection and independent full-group witness
replay.

The baseline is PR #1399 head `a6577e05039851bedd101dd3f38dba11b3ee4829`.
The frozen standard dimension-12 base has **4,054 actual distinct usable
subgroup points**, 2,027 sign-folded columns, 8,219,485 unordered pairs,
and 4,108,723 signed-sum representatives. Both arms queried the same
public T001 point
`(355fb5df7a905f16921eb,5900a390f42d290f1bbe)` and returned the
same exact `no_witness` outcome. This is expected: the uniform-target
four-summand coverage ceiling is `4.662e-12` for this base. The full-base
planted `[0,2,4,6]` control returned the same witness from the portable
and packed queries, and both full curve sums replayed to the target.

| Full-base run order | Arm | Index build (ms) | Exact query (ms) | Peak RSS (B) | Outcome |
| ---: | --- | ---: | ---: | ---: | --- |
| 1 | baseline | 12,419.392 | 1,888.020 | 1,270,218,752 | exact miss |
| 2 | compact | 2,389.316 | 1,004.642 | 684,359,680 | exact miss |
| 3 | compact | 5,903.493 | 1,374.830 | 678,559,744 | exact miss |
| 4 | baseline | 9,600.969 | 2,437.565 | 1,270,235,136 | exact miss |

All six panel processes exited zero; their individual
[JSONL files and stderr](status.tsv) preserve the raw outputs. The
two-point medians are **11,010.181 versus 4,146.405 ms build** and
**2,162.792 versus 1,189.736 ms query**, baseline versus compact.
Those give exploratory ratios of 2.655× build and 1.818× query. The
adjacent baseline/compact query ratios are 1.879× and 1.773×. Maximum
RSS fell **46.1%**. The [native layout receipt](layout.txt) reports a
208-byte old entry and a 64-byte packed entry on this build. Smaller
dimension-8 and dimension-10 three-repeat probes also produced identical
exact misses and representative counts; all individual times remain in
[`baseline_small.jsonl`](baseline_small.jsonl) and
[`candidate_small.jsonl`](candidate_small.jsonl). The wide timing drift
across runs is preserved, not smoothed away.

All nine focused release geometry tests passed, including ordinary,
same-x, inverse, repeated-point, and identity cases; exhaustive
small-base four-sum membership for 41 target scalars; and exact
stored-word residual keys against the reference. See
[`geometry_tests.log`](geometry_tests.log), the successful
[`build.log`](build.log), [`planted_build.log`](planted_build.log),
[`planted.jsonl`](planted.jsonl), and [`planted.status`](planted.status).
The frozen shell [runner](run.sh) used a 120-second cap per process.
Both release probes were built before the panel; no build ran between
arms. Baseline binary SHA-256 values are
`8923553ad79d2647383fd8a929b12fd4ddbfe0f880d95b9d641de295fc6e541b`
(full) and
`83b8f59d66422cfcbbe47dd840c6b476fac64cec6baea74d0c7debb6805a521c`
(small). Compact binaries are
`4b5a292164aec2e4584dec1f2470f694f78b71f510eabf47782e3f33afe093d7`
and `76400560b6481206d5e6934c7ce7fa7faf4ebc51bbe03f7b2e65f6715cff1ab8`.
The final source SHA-256 is
`2cd2962ed9375bc235291a3a43364e9e0971d62b0e24d22ce715e5745261377f`;
the protocol and runner hashes are
`a85a988c893c335dbc8d235627e893a27ef2785167a93e6d54d03092cf60ef83`
and `f909970606d326729c6deb45d9677662993a85999078aeebdbf649af9d774041`.
The [SHA-256 receipt](SHA256SUMS) covers each raw artifact.

The host was physical Apple M4 Pro arm64, Darwin macOS 26.6, Rust 1.93.1.
It had no auditable exclusive CPU partition; several other builds ran
during this work. The timing ratios are **exploratory component diagnostics**
and do not satisfy the isolated CPU speedup gate. The query is
target-dependent, while index construction is reusable preparation.
Neither row is a complete n83 F6 point-decomposition attempt. There is no
verified ordinary six-summand relation, F4/F5/F6 complete-call comparison,
one-target IC online interval, matched one-target rho result, or end-to-end
speedup. Those values remain unknown.

Decision: retain the compact index representation because exactness,
46.1% lower peak RSS, 45.0% lower median query time, and faster build
meet every registered gate. The next n83 search still requires a solver
for natural six-summand relations; this four-summand component alone has
negligible natural-target coverage.
