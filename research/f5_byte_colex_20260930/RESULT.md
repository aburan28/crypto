# Byte-colex F5 row packing: stage gain, full-call gate failed

The [protocol](PROTOCOL.md) and [unmodified baseline build](BASELINE.json)
preceded candidate code and timing in draft PR #1074. The candidate used a
small per-call byte lookup for colex column indices only when direct packed
row construction was selected. It remained opt-in throughout testing.

The exactness screen passed on both local seeds and all seven benchmark cases.
Every returned-row fingerprint, canonical row-space fingerprint, term count,
rank, matrix shape and counted operation matched. The exhaustive degree-four
column test passed; F5 release tests passed 10/10 and GF(2) elimination tests
6/6. The [local paired screen](LOCAL_RESULT.md) advanced the candidate to the
qualified x86 run but was exploratory because macOS could not reserve a CPU.

The [Linux x86 workflow](https://github.com/aburan28/crypto/actions/runs/36760909104)
passed both one- and two-thread jobs. Each selected seed block had five A/A
and five alternating A/B pairs on one binary, an empty eligible-user-thread
set on reserved CPUs, zero contended samples, and exact output on every call.
The one-thread runner used AMD EPYC 9V74; the two-thread runner used AMD EPYC
7763. Times below are medians of complete F5 calls in milliseconds, and ratios
are medians of paired reference/candidate ratios. They are compared only
within their own host and thread count.

| Threads | Seed XOR | Reference ms | Candidate ms | Build ratio | Full-call ratio | 95% paired-median interval |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 1 | `0` | 85.988 | 83.948 | 1.458 | 1.0243 | 1.0168–1.0298 |
| 1 | `badc0de1` | 86.143 | 83.940 | 1.463 | 1.0262 | 1.0239–1.0309 |
| 1 | `5eed2026` | 84.124 | 82.304 | 1.457 | 1.0229 | 1.0076–1.0308 |
| 1 | `f5c02a28` | 85.133 | 82.618 | 1.456 | 1.0241 | 1.0212–1.0304 |
| 2 | `0` | 94.994 | 92.398 | 1.342 | 1.0281 | 1.0094–1.0641 |
| 2 | `badc0de1` | 94.252 | 93.051 | 1.324 | 1.0093 | 0.9978–1.0134 |
| 2 | `5eed2026` | 94.145 | 92.557 | 1.321 | 1.0172 | 1.0075–1.0344 |
| 2 | `f5c02a28` | 94.862 | 92.340 | 1.310 | 1.0272 | 1.0222–1.0316 |

Each interval is the exact 3,125-resample bootstrap for that seed's five
paired ratios. The table uses the first clean block per seed, selected before
reading its ratios. Four one-thread preflight attempts failed CPU-pressure
qualification before launching a measured child. One two-thread block
completed 22 calls but had one contended sample and was discarded; its
replacement was clean. Every failure and discarded call is retained in the
raw archives.

The preregistered one-thread full-call median had to reach **1.05×** with a
lower interval bound above 1.02. All four medians were 1.023–1.026×, so the
candidate fails that gate despite a 1.46× row-build phase gain. It also fails
the smaller-case guard in the two-thread `badc0de1` block: n24 degree 3 had
a 0.9358× paired median versus its A/A minimum of 0.9878×. No complete-call
2× result was observed. The runtime option was removed after archiving.

The exact measured source diff is
[`measured_candidate.patch.gz`](measured_candidate.patch.gz), SHA-256
`0ca522cabea38da7486b328408c45fc55a6dd411ce4b17f2a2e6346ca2f480c8`.
The [raw one-thread archive](runs/x86_36760909104/segments-t1.tar.gz)
contains 34 files, SHA-256
`3144f35a6628c4144aca5aaffebbcfe63735bb25a5c886a28c64774c54c3cb06`;
the [raw two-thread archive](runs/x86_36760909104/segments-t2.tar.gz)
contains 27 files, SHA-256
`d0d6ca8ec14b9d0d1344fc104ac8554b355c273338f39bd26d9673de3b6e6fbd`.
Adjacent `MANIFEST-t1.json` and `MANIFEST-t2.json` list every file, byte
count and digest. This is a Boolean matrix-F5 solver-stage diagnostic, not
one-target IC online time or a paired Pollard-rho comparison.
