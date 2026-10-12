# Stage 38: algebraic relation-collection window sweep

Stage 38 integrates the windowed, walked relation collector developed in draft PR #204 with the current Stage-35 public-target and evidence boundary. The upstream draft measured pseudo-random subgroup bases; this stage tests only the algebraically discovered Frobenius-union/two-torsion-saturation base.

A full three-summand pair-table scan finds a triple once with each summand in the scanned position, while retaining only one sorted witness. A cyclic window keeps any witness whose first, second, or third index falls in the window, then verifies the unsorted witness in the group. Additional probes are walked in deterministic runs of 64, so each run pays one scalar multiplication and shared additions. Work-unit boundaries still reproduce the same relation sequence.

The integration retains Stage 35's explicit `descent_summands` report and adds the effective `collection_window`. It also carries `summands_scanned` in every relation-unit document and reports initial and final totals. Older units deserialize the new counter as zero; new measurement admission requires the exact nonzero charge `trials × window`, or `trials × |F|` for the full-scan control.

The projected signed-Frobenius orbit map now uses single-word group arithmetic when the field fits in a word. A point-for-point equivalence test compares it with the general big-integer map on several subgroup bases and the exact n=41 algebraic recipe. The first orbit-map construction previously fell between the select and collect timers; fresh workflows now charge it to selection.

The hosted sweep runs six cells on one inherited Linux CPU from one clean binary:

- full scan;
- windows 595, 297, 149, 74, and 37;
- 8,192 probes per work unit, one initial unit, and at most 64 units;
- the exact Stage-35 algebraic search, five public targets, two-summand descent, and signed-Frobenius rho in every cell.

The winner minimizes `collection + logs`, the two stages the window changes. Every cell still charges the repeated factor-base search and five-target descent, but their timing noise cannot choose this parameter. Selection occurs only after all six cells finish and does not change their public target stream or factor-base candidates. A later holdout must confirm any chosen window.

A local same-binary feasibility sweep found all six cells complete. Relative to the new full-scan control, windows reduced the final third-summand lookup charge while keeping 29 certified columns. Window 595 gave the lowest local precompute, 3.661 seconds versus 4.373 seconds full scan; window 149 gave 3.632 seconds in a consecutive run. These timings select no claim: the hosted sweep determines the candidate and a separate run must validate it.

This is a setup constant optimization. It does not replace the SAT matrix, change the algebraic factor-base requirement, supply Magma or external review, change the asymptotic exponent, or establish a Koblitz index-calculus SOTA result.

## Hosted result

GitHub run [`34698066872`](https://github.com/aburan28/crypto/actions/runs/34698066872) completed from merge commit `647481caed3c703d270800922d1e262a5b683834`. The downloaded artifact and a fresh local replay of its sealed verifier agree. Every cell reconstructed the same 4,759-point algebraic base, certified all 29 factor-base logarithms from collected relations, and solved the same five public hash-to-curve targets without constructed target scalars.

| Window | Relations | Trials | Summand lookups | Select (s) | Collect + logs (s) | Precompute (s) | Five-target IC (s) | IC/rho amortized |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full | 42 | 16,384 | 77,971,456 | 5.684 | 10.007 | 15.692 | 16.115 | 31.58x |
| 595 | 43 | 32,768 | 19,496,960 | 5.691 | 3.883 | 9.574 | 9.996 | 19.84x |
| 297 | 56 | 81,920 | 24,330,240 | 5.679 | 4.502 | 10.181 | 10.590 | 20.76x |
| **149** | **31** | **81,920** | **12,206,080** | **5.637** | **3.361** | **8.999** | **9.418** | **18.45x** |
| 74 | 40 | 172,032 | 12,730,368 | 5.621 | 3.725 | 9.346 | 9.770 | 19.26x |
| 37 | 29 | 294,912 | 10,911,744 | 5.641 | 3.886 | 9.526 | 9.902 | 19.46x |

Window 149 wins the predeclared `collection + logs` metric. Against the same-run full scan it is 2.98 times faster on those stages, performs 6.39 times fewer charged third-summand lookups, and reduces the five-target IC time by 1.71 times. Against Stage 35, its five-target IC time improves by 1.54 times and the amortized IC/rho ratio improves from 28.86 to 18.45. Its online descent remains 1.22 times faster than signed-Frobenius rho, while the complete five-target algorithm remains slower.

The six-cell scientific sweep used 69.708 seconds wall, 69.564 core-seconds, and at most 346,611,712 bytes sampled process-tree RSS. Including the fresh four-core build gives 167.434 seconds sequential wall, 396.613 core-seconds, and 1,667,104,768 bytes maximum sampled process-tree RSS. The complete workflow, including the validation job and runner queue gap, took 319 seconds wall. These sweep-selection resources are charged evidence; Stage 39 must measure the selected window once on untouched relation and target streams before it can replace the Stage 35 result.
