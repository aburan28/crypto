# Orbit-disjoint v2 cold panel: six eligible CPU no-gos

The one-shot [run 36803331080](https://github.com/aburan28/crypto/actions/runs/36803331080) measured the Q panel frozen by #1102 and validated by #1105, on `main` at `38c173072e16f580237447c93b085634275ba048`. All six `measure` jobs completed. A second host replayed every sealed rank and public-Q logarithm and reproduced the hosted ratios. That replay did not retime the cells.

The single-host verifier still records `timing_eligible: false` and `second_host_replay_required: true`. The preregistered rule, applied after that replay, accepts all six cells. Each has `contended_samples = 0`, a compact A/A median inside `[0.90, 1.10]` whose 95% interval contains 1, and `host_timing_candidate: true` on both the hosted receipt and the replay. Idle thread affinity stays in the isolation record and is not treated as measured contention.

**Every eligible interval lies wholly above 1.** These six W64 policies are fixed-policy native-CPU no-gos. They do not close other index-calculus routes. Common operation-equivalent `S` stays unset. There is no method speedup, no n=83 result, and no n=131 transfer.

| cell | K | compact CPU / rho [95%] | A/A median | verdict |
|:--|--:|--:|--:|:--|
| n37/L1 | 7 | 2.810 [2.718, 2.874] | 0.996 | no-go |
| n37/L1024 | 42 | 3.031 [2.998, 3.047] | 0.999 | no-go |
| n41/L1 | 85 | 13.378 [9.367, 18.491] | 0.996 | no-go |
| n41/L1024 | 255 | 1.314 [1.288, 1.359] | 1.000 | no-go |
| n53/L1 | 220 | 22.675 [17.904, 31.595] | 1.006 | no-go |
| n53/L1024 | 440 | 1.605 [1.566, 1.666] | 0.992 | no-go |

The ratio is the preregistered geometric mean of the two compact child CPU times divided by the same-block rho CPU time. Class is accounting: the compact and rho programs are the frozen sources, and this round prices that frozen policy under the v2 isolation rule. The hosted cold-run SHA-256 values and the matching replay are in [`SUMMARY.json`](SUMMARY.json). Do not dispatch this panel again and do not replace a Q.
