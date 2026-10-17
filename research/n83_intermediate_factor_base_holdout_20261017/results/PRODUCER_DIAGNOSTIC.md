# N83 intermediate factor-base producer diagnostic

K3072/root11 produced the shortest target-online geometric mean across the
20 sealed public points in cloud run `39042774613`: 6.205612 ms. The table
uses the per-target producer `target_ms` fields from all four archived
strata, paired by point index. Each rho/IC ratio uses a separately solved
same-point rho receipt. These are producer diagnostics while independent
replay of the preserved archives completes; the frozen selector has made no
choice.

| Factor-base policy | Points | Geometric target-online ms | Geometric same-point rho/IC | K3072 faster point count | Geometric paired policy/K3072 |
| --- | ---: | ---: | ---: | ---: | ---: |
| K2304/root11 | 20 | 30.605300 | 49.945708 | 18/20 | 4.931875 |
| K2543/root11 | 20 | 17.279708 | 143.728003 | 15/20 | 2.784529 |
| K2816/root11 | 20 | 18.418924 | 76.196215 | 15/20 | 2.968108 |
| K3072/root11 | 20 | 6.205612 | 308.309304 | — | 1.000000 |

All raw before/after calibration receipt pairs were within the 10% gate,
including when relative spread is divided by the smaller of the two times.
The original cloud verifier stopped at its module lookup error before
independent relation and scalar replay. The corrected verifier's local
replays and provider binding determine the admissible result. Cold setup,
index construction, rank work, resources, and per-point statuses remain in
the raw archives listed in `PROVIDER_ARCHIVES.md`.
