# Boolean F4 lookup-table width study

The opt-in six-column/four-matrix-budget arm matched exact bases, counted logical XORs and matrix shapes on all seven cases and both holdouts. It passed the primary one-thread timing gate but has an unresolved smaller-case holdout regression, so the four-column default remains in place pending an independent repeat and four-thread control. This is an internal Boolean F4 solver-stage result, not a one-target IC or DLP speedup.

CI run [36297452737](https://github.com/aburan28/crypto/actions/runs/36297452737) at PR head `e0d1f171` passed tests for both table widths and completed 66 calls. The [complete raw receipt](runs/36297452737/ci-result.json), SHA-256 `ddc877f55e704b154199483a4e4f87d6d7276e4afcbcc59f2770225b19257435`, contains all process outcomes, paired timings, exact output signatures, operation counts and actual table storage. The Linux x86-64 runner reported AMD EPYC 7763, Rust 1.98.1 and one pinned CPU.

| `n20_m30` workload | Four-column elimination | Six-column elimination | Paired elimination ratio, 95% interval | Four-column full F4 | Six-column full F4 | Paired full ratio, interval |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Frozen seed | 302.6 ms | 279.6 ms | 1.082, 1.037–1.089 | 571.7 ms | 531.3 ms | 1.074, 1.048–1.078 |
| Holdout `1ac0ffee` | 399.8 ms | 305.0 ms | 1.303, 1.278–1.351 | 664.7 ms | 555.2 ms | 1.192, 1.177–1.216 |
| Holdout `2468ace0` | 441.6 ms | 340.2 ms | 1.339, 1.211–1.473 | 712.7 ms | 665.1 ms | 1.176, 1.003–1.304 |

Milliseconds are each arm's five A/B median; the paired ratio median can differ from the ratio of marginal medians. The frozen primary matrix peak was 44.4 MB in both arms, while actual table peak rose from 21.4 to 55.2 MB; holdout A tables rose from 36.9 to 94.5 MB. Performed word XORs on the frozen primary fell from 371.6 to 308.1 million, including table construction. These extra tables stayed within the candidate's declared four-matrix budget.

The holdout B `n20_m20` complete-call paired median was 0.881 against its A/A range 0.891–1.067, with bootstrap interval 0.851–1.007. All other smaller-case medians were within A/A noise or improved. The runner's holdout B A/A ranges were wide across multiple cells, so this one result cannot establish a repeatable regression or be ignored. The next frozen control repeats all cases at one thread and adds four threads; the feature remains opt-in until those results settle the gate.

## Independent one-thread repeat and four-thread control

CI run [36298595953](https://github.com/aburan28/crypto/actions/runs/36298595953) at PR head `0d14ad99` completed 66 calls at each thread count with exact outputs and logical counters. The raw [one-thread](runs/36298595953/ci-result.json) and [four-thread](runs/36298595953/ci-result-threads4.json) receipts have SHA-256 `05d1078078c349f1c320da9fc11727cebceba0938b8377caad1f3bfa0b881a38` and `a41ac3811f2cffc2a0a1d881a167c7abaf2559b697026334be45fcd0327c761a`, respectively. This run used a Linux x86-64 AMD EPYC 9V74 runner; comparisons stay within each run and thread count.

| Thread count, `n20_m30` | Frozen elimination ratio, 95% interval | Frozen full-F4 ratio, interval | Holdout A full ratio | Holdout B full ratio |
| --- | ---: | ---: | ---: | ---: |
| One | 1.096, 1.088–1.111 | 1.078, 1.053–1.086 | 1.206 | 1.255 |
| Four | 1.095, 1.076–1.113 | 1.090, 1.081–1.110 | 1.130 | 1.114 |

The primary gain repeats. Yet `n20_m20` at one thread now has a frozen complete-call ratio of 0.982 against A/A 0.993–1.007, and holdout A is 0.994 against A/A 0.997–1.006. Holdout B is within A/A, as are all four-thread smaller-case medians. The smaller-case slowdown recurred on a different runner and workload, so six-column tables **fail the preregistered no-regression gate** and remain opt-in. The next five-column arm is a distinct experiment, not a reinterpretation of the six-column result.

## Five-column first paired run

CI run [36300043072](https://github.com/aburan28/crypto/actions/runs/36300043072) at PR head `2a0c4518` passed five-column correctness tests and completed 66 one-thread calls with exact bases, counted logical XORs and matrix shapes. The [full raw receipt](runs/36300043072/ci-result-five.json), SHA-256 `5654da89e41063189ec096ce200f7b8fdf69be9946c2f110787c0d1582d4689a`, records a Linux x86-64 AMD EPYC 7763 runner pinned to one CPU.

| `n20_m30` workload | Four-column elimination | Five-column elimination | Paired elimination ratio, 95% interval | Four-column full F4 | Five-column full F4 | Paired full ratio, interval |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Frozen seed | 302.6 ms | 277.0 ms | 1.095, 1.075–1.101 | 567.2 ms | 544.3 ms | 1.043, 1.035–1.046 |
| Holdout `1ac0ffee` | 397.3 ms | 299.7 ms | 1.326, 1.307–1.332 | 661.2 ms | 565.0 ms | 1.170, 1.163–1.177 |
| Holdout `2468ace0` | 406.0 ms | 297.3 ms | 1.367, 1.349–1.376 | 669.9 ms | 562.0 ms | 1.192, 1.182–1.196 |

Frozen primary peak table storage rose from 21.4 to 33.6 MB while the matrix peak stayed 44.4 MB; performed XORs fell from 371.6 to 324.9 million. No smaller complete-call median regressed beyond its own A/A range. The frozen A/A full-F4 range was unusually wide, 0.991–1.240, despite a paired A/B interval above one. An independent one-thread repeat and four-thread control are required before any default change.

## Five-column independent repeat and four-thread control

CI run [36301067360](https://github.com/aburan28/crypto/actions/runs/36301067360) at PR head `659a9cfc` completed 66 calls at each thread count. Every call succeeded with the same exact output signatures and logical counts. The complete raw [one-thread](runs/36301067360/ci-result-five.json) and [four-thread](runs/36301067360/ci-result-five-threads4.json) receipts have SHA-256 `cc68e3ed1601618176a857203e2f534aeeacb6c997882045e6290f238ada2404` and `3f69baa15538342b49ea239304f80825895dcb0b9d5883822172198ef9062833`.

| Threads, `n20_m30` workload | Elimination ratio, 95% interval | Complete F4 ratio, 95% interval |
| --- | ---: | ---: |
| One, frozen | 1.074, 1.068–1.095 | 1.037, 1.033–1.049 |
| One, holdout A | 1.298, 1.291–1.318 | 1.155, 1.148–1.168 |
| One, holdout B | 1.324, 1.299–1.331 | 1.166, 1.156–1.174 |
| Four, frozen | 1.060, 1.008–1.071 | 1.020, 1.007–1.053 |
| Four, holdout A | 1.180, 1.151–1.224 | 1.060, 1.050–1.103 |
| Four, holdout B | 1.163, 1.141–1.205 | 1.071, 1.056–1.096 |

The independent repeat confirms the primary and holdout gains, and no four-thread smaller case regressed beyond A/A noise. One one-thread smaller case misses the frozen strict gate: holdout B `n14_m14` has a complete-call paired median of 0.9964 against its A/A range 1.0007–1.0159. Its five-pair bootstrap interval is 0.9931–1.0161. This is small and uncertain, but the gate is not met. The four-column default stays in place while this specific cell is adjudicated; five and six columns remain opt-in.

## Bounded five-column adjudication

CI run [36302666376](https://github.com/aburan28/crypto/actions/runs/36302666376) at PR head `8fe9ac2f` completed 66 successful calls at each thread count, with exact output signatures and logical counters. The raw [one-thread](runs/36302666376/ci-result-five.json) and [four-thread](runs/36302666376/ci-result-five-threads4.json) receipts have SHA-256 `ceef082d635505b2b41c413bda3fe88ff5373a72ce95220dba05a64321491af7` and `947ef7d141d1649cfa505f7cf8c4f0b3177a233980566655895d611ca23fd2bd`.

| Threads, `n20_m30` workload | Elimination ratio, 95% interval | Complete F4 ratio, 95% interval |
| --- | ---: | ---: |
| One, frozen | 1.073, 1.039–1.105 | 1.032, 1.008–1.051 |
| One, holdout A | 1.314, 1.300–1.347 | 1.165, 1.157–1.174 |
| One, holdout B | 1.363, 1.330–1.377 | 1.184, 1.169–1.196 |
| Four, frozen | 1.061, 0.984–1.062 | 1.010, 0.995–1.026 |
| Four, holdout A | 1.176, 1.144–1.190 | 1.071, 1.045–1.085 |
| Four, holdout B | 1.177, 1.115–1.241 | 1.081, 1.024–1.098 |

The disputed one-thread holdout B `n14_m14` median was 0.9923, inside its own A/A range 0.9779–1.0042; its bootstrap interval was 0.9842–1.0031. No other smaller median fell below its A/A range at either thread count. The one-thread frozen primary clears the 1.05 elimination gate; the four-thread control does not show a material regression, though its frozen intervals include one. The bounded adjudication clears the frozen promotion rule, so the next change selects five columns by default with a four-column reference and a final default-on replay. Six columns remain opt-in after its repeated smaller-case regression.

## Final default-on replay with batched pairs

CI run [36305490281](https://github.com/aburan28/crypto/actions/runs/36305490281) at PR head `11345661` completed 66 successful calls at each thread count. Both binaries include the merged batched-pair default from PR #882. The reference explicitly selects four-column tables; the candidate uses the new five-column default. Every case matched its exact basis, matrix shape and logical work counters. The raw [one-thread](runs/36305490281/ci-result-five.json) and [four-thread](runs/36305490281/ci-result-five-threads4.json) receipts have SHA-256 `6e07266627baea542c3f0afa6b8bc81bac2249b235072a668f0ab752a0bc5844` and `8a011a6d3dbd8bc3d354b10a55f28370e2350085a9ac656c78e44eb170847988`.

| Threads | `n20_m30` workload | Elimination ratio, 95% interval | Complete F4 ratio, 95% interval |
| ---: | --- | ---: | ---: |
| 1 | Frozen | 1.066, 1.057–1.072 | 1.079, 1.070–1.137 |
| 1 | Holdout A | 1.359, 1.350–1.367 | 1.237, 1.231–1.242 |
| 1 | Holdout B | 1.354, 1.349–1.368 | 1.233, 1.229–1.247 |
| 4 | Frozen | 1.058, 1.046–1.066 | 1.068, 1.059–1.096 |
| 4 | Holdout A | 1.159, 1.132–1.201 | 1.136, 1.104–1.144 |
| 4 | Holdout B | 1.207, 1.174–1.219 | 1.139, 1.099–1.149 |

Both holdout primary gains lie outside their same-run A/A ranges. No smaller complete-call median fell below its own A/A lower bound at either thread count. The one-thread frozen elimination and complete-call gates, the four-thread control, and exact-output checks all pass. Five-column tables are the measured default within the stated table budget; four columns remain available with `f4-four-tables`, while six columns remain opt-in with `f4-wide-tables`. These are Boolean F4 solver-stage gains, not a measured one-target IC or DLP speedup.
