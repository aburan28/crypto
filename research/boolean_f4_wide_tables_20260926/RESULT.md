# Boolean F4 six-column lookup tables: first paired run

The opt-in six-column/four-matrix-budget arm matched exact bases, counted logical XORs and matrix shapes on all seven cases and both holdouts. It passed the primary one-thread timing gate but has an unresolved smaller-case holdout regression, so the four-column default remains in place pending an independent repeat and four-thread control. This is an internal Boolean F4 solver-stage result, not a one-target IC or DLP speedup.

CI run [36297452737](https://github.com/aburan28/crypto/actions/runs/36297452737) at PR head `e0d1f171` passed tests for both table widths and completed 66 calls. The [complete raw receipt](runs/36297452737/ci-result.json), SHA-256 `ddc877f55e704b154199483a4e4f87d6d7276e4afcbcc59f2770225b19257435`, contains all process outcomes, paired timings, exact output signatures, operation counts and actual table storage. The Linux x86-64 runner reported AMD EPYC 7763, Rust 1.98.1 and one pinned CPU.

| `n20_m30` workload | Four-column elimination | Six-column elimination | Paired elimination ratio, 95% interval | Four-column full F4 | Six-column full F4 | Paired full ratio, interval |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Frozen seed | 302.6 ms | 279.6 ms | 1.082, 1.037–1.089 | 571.7 ms | 531.3 ms | 1.074, 1.048–1.078 |
| Holdout `1ac0ffee` | 399.8 ms | 305.0 ms | 1.303, 1.278–1.351 | 664.7 ms | 555.2 ms | 1.192, 1.177–1.216 |
| Holdout `2468ace0` | 441.6 ms | 340.2 ms | 1.339, 1.211–1.473 | 712.7 ms | 665.1 ms | 1.176, 1.003–1.304 |

Milliseconds are each arm's five A/B median; the paired ratio median can differ from the ratio of marginal medians. The frozen primary matrix peak was 44.4 MB in both arms, while actual table peak rose from 21.4 to 55.2 MB; holdout A tables rose from 36.9 to 94.5 MB. Performed word XORs on the frozen primary fell from 371.6 to 308.1 million, including table construction. These extra tables stayed within the candidate's declared four-matrix budget.

The holdout B `n20_m20` complete-call paired median was 0.881 against its A/A range 0.891–1.067, with bootstrap interval 0.851–1.007. All other smaller-case medians were within A/A noise or improved. The runner's holdout B A/A ranges were wide across multiple cells, so this one result cannot establish a repeatable regression or be ignored. The next frozen control repeats all cases at one thread and adds four threads; the feature remains opt-in until those results settle the gate.
