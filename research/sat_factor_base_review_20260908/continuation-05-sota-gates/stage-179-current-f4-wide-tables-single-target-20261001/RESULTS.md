# Stage 179: six-column BlockTables on the frozen target

Both feature configurations passed all eight Boolean-F4 tests and all six target processes returned exhaustive UNSAT with identical algebraic and row-equivalent F4 counters.

| pair | six / five wall | six / five core | six / five RSS |
|---:|---:|---:|---:|
| 1 | 1.507128 | 1.032098 | 0.862692 |
| 2 | 0.797520 | 1.049225 | 1.019843 |
| 3 | 1.108173 | 0.976286 | 0.940548 |

Median paired ratios are 1.108173 wall, 1.032098 CPU, and 0.940548 RSS. Severe host descheduling makes absolute wall snapshots unsuitable, but the candidate also performs 1.106291x as many actual word XORs and uses 1.396484x the table memory.

The two clean builds, two feature-test runs, and six paired queries charge 2023.503767 wall seconds, 3234.101524 core-seconds, and 4237049856 bytes maximum RSS.

Six columns is rejected and the five-column repository default remains selected. This does not change any SOTA gate.
