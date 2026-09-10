# Stage 8: degree-15 replication panel

The selected `[2,4]` divisor kernel was reused without target-dependent changes for five new scalar-blind targets. Every index-calculus process and every same-target automorphism-rho process ran sequentially under a fresh-process wall/CPU/RSS meter.

| Secret | IC relations / trials / conflicts | IC wall / core-s / peak MiB | rho iterations / additions | rho wall / core-s / peak MiB | IC/rho core ratio |
|--:|:--|:--|:--|:--|--:|
| 17 | 7 / 7 / 3,779 | 0.209406 / 0.198017 / 3.61 | 1 / 20 | 0.005237 / 0.003591 / 2.02 | 55.14 |
| 53 | 7 / 7 / 5,231 | 0.205382 / 0.203253 / 3.77 | 2 / 23 | 0.005332 / 0.003631 / 2.05 | 55.98 |
| 89 | 7 / 7 / 5,235 | 0.209155 / 0.206919 / 3.69 | 3 / 26 | 0.005535 / 0.003872 / 2.02 | 53.44 |
| 137 | 7 / 7 / 2,288 | 0.194317 / 0.192146 / 3.50 | 1 / 20 | 0.005229 / 0.003706 / 2.05 | 51.85 |
| 199 | 7 / 7 / 3,238 | 0.194118 / 0.191941 / 3.59 | 3 / 26 | 0.005249 / 0.003743 / 2.03 | 51.28 |

All five scalars were recomputed from the recovered values and matched their public targets. All five runs had zero direct-relation recoveries and zero invalid models. Median index-calculus cost was 0.198017 core-seconds and median rho cost was 0.003706 core-seconds, for a median per-pair ratio of 53.44. Total measured core time was 0.992276 seconds for index calculus and 0.018543 seconds for rho.

This panel establishes repeatability at one toy degree and one public factor base. It does not establish asymptotic scaling, a crossover, external reproduction, novelty, or SOTA.
