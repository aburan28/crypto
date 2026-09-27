# Matrix-F5 echelon output: first paired result

The always-echelon arm preserved the exact canonical row-space fingerprint, rank, pruning and criterion counts for every seven-case call on the frozen seed and two holdouts. Its raw row fingerprint changed as expected; on the frozen primary it changed from `ed5234ba018bc079` to `3f659516eff553b8` while the canonical row-space fingerprint remained `ed5234ba018bc079`.

CI run [36304304816](https://github.com/aburan28/crypto/actions/runs/36304304816) at PR head `e5c6f3e7` completed 66 successful calls on a pinned one-thread AMD EPYC 7763 Linux runner. The [full raw receipt](runs/36304304816/ci-result.json) has SHA-256 `4b4026da1f5c429b9b6eb124c8f56f3279acba9244150ee44b35fc618a768b26`. Ratios are paired reference/candidate medians with exact five-pair bootstrap 95% intervals.

| `f5_n24_m24_d4` workload | Elimination ratio, 95% interval | Complete F5 ratio, 95% interval |
| --- | ---: | ---: |
| Frozen | 2.330, 2.286–2.342 | 1.519, 1.505–1.525 |
| Holdout A | 2.295, 2.252–2.304 | 1.510, 1.505–1.522 |
| Holdout B | 2.375, 2.321–2.395 | 1.548, 1.526–1.565 |

This clears the primary gain gate but fails the frozen no-smaller-regression gate. The `n12_m12_d4` complete-call ratio was 0.830 on the frozen seed, 0.821 on holdout A and 0.823 on holdout B, all below A/A noise, despite its elimination ratio exceeding one. The `n16_m16_d4` and one `n24_m24_d3` cell also missed their A/A ranges. The extra unpacking of some echelon bases outweighed their faster elimination. Always-echelon output therefore stays explicit and is not selected by default. A distinct structural selection rule is specified in the protocol before its independent run. These are internal solver-stage costs, not one-target IC or DLP speedups.

## Selective one-thread result

CI run [36305037801](https://github.com/aburan28/crypto/actions/runs/36305037801) at PR head `51e1abcc` completed 88 successful calls over the frozen seed and three holdouts. All canonical row-space fingerprints, ranks, pruning and criterion counts matched. The [complete raw receipt](runs/36305037801/ci-result.json) has SHA-256 `1b4d379329b0868d5ca4f6908e51900ff0ed717ebb40b9d7ce734f6f6a28d7fa`.

| `f5_n24_m24_d4` workload | Elimination ratio, 95% interval | Complete F5 ratio, 95% interval |
| --- | ---: | ---: |
| Frozen | 2.321, 2.305–2.348 | 1.504, 1.496–1.510 |
| Holdout A | 2.246, 2.160–2.259 | 1.488, 1.456–1.496 |
| Holdout B | 2.349, 2.340–2.374 | 1.540, 1.537–1.543 |
| Holdout C | 2.318, 2.292–2.341 | 1.527, 1.518–1.535 |

The selective rule keeps the large gain and removes the material `n12_m12_d4` regression. Three unchanged reduced-path cells narrowly missed the strict A/A gate: holdout A `n16_m16_d4` 0.9907 versus A/A minimum 0.9937; holdout A `n24_m24_d3` 0.9913 versus 0.9917; and holdout C `n24_m24_d3` 0.9863 versus 0.9915. Because both arms execute the same reduced code for these cells, this is an unresolved timing control, not evidence of a changed row space. The frozen rule still treats it as a miss. One bounded independent repeat and four-thread control are specified before the next run; the default remains RREF.
