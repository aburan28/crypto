# Stage 5: GitHub-hosted external-environment reproduction

GitHub Actions run [34426104136](https://github.com/aburan28/crypto/actions/runs/34426104136) completed successfully on `ubuntu-latest` (`Linux 6.17`, Azure x86-64, Python 3.12.3). The workflow built from a fresh checkout, ran both scalar-blind index-calculus controls and both signed-Frobenius quotient-rho controls, applied a fail-closed verifier, passed the 34 focused Koblitz library tests, and uploaded the retained metrics and dependency lock as `koblitz-sota-gate-external-reproduction`.

| Secret | IC relations/trials/conflicts | IC wall / core-s / peak MiB | rho iterations/additions | rho wall / core-s / peak MiB |
|--:|:--|--:|:--|--:|
| 53 | 5 / 5 / 122 | 0.004259 / 0.003957 / 13.34 | 2 / 23 | 0.002392 / 0.002120 / 13.34 |
| 101 | 5 / 5 / 125 | 0.004335 / 0.004016 / 13.36 | 4 / 29 | 0.002454 / 0.002180 / 13.32 |

The scientific counters match the macOS producer runs exactly. The external host recovered and point-verified both scalars, observed no direct-relation shortcut or invalid SAT model, confirmed that the factor-base predicate uses no discrete-log labels or subgroup enumeration, and recovered the same targets with the quotient rho walk. Rho again consumed fewer core-seconds in both pairs.

This is an external-host correctness and accounting reproduction, not an independent design: the workflow and verifier are part of this PR. It does not reproduce WDSat, CryptoMiniSat, Magma, the n=31/41/59 matrix, or the public GGMP discovery. It also is not an external novelty review. Those distinctions keep gate 7 open in part.
