# Boolean F4 degree-bucket pair queue

The opt-in degree-bucket queue produced the same basis fingerprints, matrix shapes, pair counts, skip counts and logical work as the flat queue on every call. The local Boolean F4 test module passed, including the deterministic pair-selection ordering test. CI run [36306721905](https://github.com/aburan28/crypto/actions/runs/36306721905) at PR head `077d3b75` completed all 66 one-thread calls on the frozen seed and two holdouts. Its [complete raw receipt](runs/36306721905/ci-result.json) has SHA-256 `76c2573c5f85ac4c2575a554b67f29ef5f50256c1bdc60838c43d9934a4c14fb`.

| `n20_m30` workload | Outside-build/elimination ratio, 95% interval | Complete F4 ratio, 95% interval |
| --- | ---: | ---: |
| Frozen | 1.040, 1.034–1.056 | 1.002, 0.987–1.003 |
| Holdout A | 1.036, 1.032–1.041 | 1.008, 1.005–1.023 |
| Holdout B | 1.029, 1.026–1.437 | 0.993, 0.986–1.076 |

Ratios are paired flat/bucket medians. No smaller complete-call median fell below its own A/A lower bound, but the frozen targeted phase misses the predeclared 1.05 threshold and the frozen complete-call interval includes one. The one-thread promotion gate therefore fails. No four-thread control or default-on replay is warranted. The flat queue remains the default; the measured candidate source is preserved in commit `077d3b75` and the final PR removes its runtime and workflow changes. This is an internal solver-stage result, not a one-target IC or DLP speedup claim.
