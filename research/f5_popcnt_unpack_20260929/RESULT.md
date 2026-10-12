# POPCNT direct F5 output unpack: rejected after qualified paired runs

The [protocol](PROTOCOL.md) was committed before candidate code or timing. The tested implementation compiled the unchanged direct output loop with runtime-gated x86-64 POPCNT and selected it with `KIC_F5_POPCNT_UNPACK=1`. Its [measured source patch](measured_candidate.patch) and [candidate workflow](WORKFLOW_CANDIDATE.yml) are archived. The opt-in runtime path and active one-off workflow were removed after the frozen gate failed. This is a matrix-F5 solver-stage result, not a one-target IC DLP or rho speedup.

[CI run 36671280788](https://github.com/aburan28/crypto/actions/runs/36671280788) passed the release shared-elimination and F5 exactness tests, then completed one- and two-thread paired jobs with Rust 1.98.1 on Linux x86-64. The one-thread host was AMD EPYC 9V74; the two-thread host was AMD EPYC 7763. Their absolute times are not compared across hosts. Both jobs used one binary, SHA-256 `055d3fd9e09e8b4b5827b23ec2d4677ef4b9bc496241ee6db09815b65cb2d343`, for explicit reference/candidate modes. The candidate F5 source SHA-256 was `1ba6d597e31a374cd0fe0591e7fa8f2ca8a8b4045aa3838dfc01873bc5c2d405`, at PR head `422327c75dbec2e30c0b36a1d1fb5f0812ea5492`; Actions checked out merge SHA `8a37bd7c68d28f621b129f78f226a11013ce84fa`.

For each seed and thread count, the **first** isolated exact block was selected before reading any timing. Every selected block had two warmups, five A/A reference pairs and five alternating A/B pairs: 22 processes, all seven F5 cases per process. All **176 selected processes** matched raw and canonical fingerprints, rank, terms, columns, builder/criterion counts and reduction word operations. The candidate route was selected only where the frozen column threshold applied. Every selected reservation left zero eligible user threads on the reserved CPUs and recorded zero contended samples. The one-thread job also preserved five preflight refusals with zero benchmark calls; the two-thread job selected all four first attempts.

The table reports the primary n24 degree-4 **complete F5 call**. Ratios are paired reference/candidate medians; intervals are exact 3,125-resample bootstrap 95% intervals from the five A/B pairs. Values above one favor POPCNT. A/A ranges measure reference/reference variation on the same seed and host.

| Rayon threads | Seed | Reference median, ms | Candidate median, ms | Full-call ratio (95% interval) | A/A range | Unpack ratio |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 1 | frozen | 86.254 | 86.221 | **1.0005×** (0.9944–1.0067×) | 0.9792–1.0098× | 1.0073× |
| 1 | holdout A | 85.782 | 85.519 | 0.9982× (0.9968–1.0113×) | 0.9913–1.0079× | 1.0078× |
| 1 | holdout B | 85.073 | 84.326 | 1.0094× (0.9949–1.0172×) | 0.9919–1.0013× | 1.0127× |
| 1 | holdout C | 85.623 | 86.434 | 0.9950× (0.8855–1.0127×) | 0.9979–1.0053× | 1.0068× |
| 2 | frozen | 95.717 | 92.910 | **1.0258×** (1.0130–1.0378×) | 0.9759–1.0089× | 1.0761× |
| 2 | holdout A | 95.702 | 93.055 | 1.0317× (1.0215–1.0414×) | 0.9781–1.0034× | 1.0786× |
| 2 | holdout B | 95.118 | 93.188 | 1.0177× (1.0070–1.0474×) | 0.9851–1.0103× | 1.0664× |
| 2 | holdout C | 95.763 | 93.200 | 1.0245× (1.0235–1.0365×) | 0.9821–0.9984× | 1.0713× |

The frozen one-thread median and lower bound miss the preregistered 1.05× and 1.02× promotion gates. The one-thread holdout-C full-call median also falls below its own A/A minimum. Two-thread unpacking improved on its different host, but the two-thread full call stayed below 1.05× and several smaller cases regressed below their own A/A minima. POPCNT supplies no accepted incremental one-thread gain and no further 2× result.

## Raw receipts and reproduction

- [One-thread manifest](runs/36671280788/MANIFEST-t1.json) and `segments-t1.tar.gz`: all nine attempts, including five preflight refusals; deterministic archive SHA-256 `f6eaf3bcb8bbdf10c0724fc39139a525379862acdbb5fef00ee893f21435c177`.
- [Two-thread manifest](runs/36671280788/MANIFEST-t2.json) and `segments-t2.tar.gz`: four first-clean attempts; deterministic archive SHA-256 `1ee7bd7ac1ec7c870055ed2f94df45b3243fda7e980bb5d5e15d12f6897cf543`.
- [Run metadata](RUN_CANDIDATE.json) records the successful CI jobs. Each manifest lists every file's byte count and SHA-256, selected attempt, binary/source hashes, CPU details, exactness and all paired ratios.

Extract an archive with `tar -xzf segments-t1.tar.gz` or `tar -xzf segments-t2.tar.gz` from the run directory. To rebuild the measured candidate, apply `measured_candidate.patch` to the protocol's pinned source, restore `WORKFLOW_CANDIDATE.yml` to `.github/workflows/f5-popcnt-unpack-segments.yml`, and replay the same-binary `=0/1` modes. Absolute times from a different host are not a denominator for the paired ratios.
