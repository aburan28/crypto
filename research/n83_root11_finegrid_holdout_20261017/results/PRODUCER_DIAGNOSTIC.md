# N83 root11 finegrid producer diagnostic

K3072/root11 produced the shortest target-online geometric mean across the
24 sealed public points in cloud run `39043441757`: 8.934360 ms. The table
uses the per-target producer `target_ms` fields from all four archived
strata, paired by point index. Each rho/IC ratio uses a separately solved
same-point rho receipt. These are producer diagnostics while independent
replay of the preserved archives completes.

| Factor-base policy | Points | Geometric target-online ms | Geometric same-point rho/IC | K3072 faster point count | Geometric paired policy/K3072 |
| --- | ---: | ---: | ---: | ---: | ---: |
| K2543/root11 | 24 | 22.056770 | 67.045900 | 15/24 | 2.468759 |
| K2816/root11 | 24 | 15.638650 | 144.457600 | 15/24 | 1.750395 |
| K3072/root11 | 24 | 8.934360 | 219.605800 | — | 1.000000 |

Stratum 3 has two raw before/after calibration pairs above 10% at point 18:
K3072 at 12.471% and K2816 at 13.355%, measured relative to the smaller
time. Its host gate must be evaluated by the corrected independent verifier;
the full 24-point Tier 1 selector remains censored if that stratum fails.
The original cloud verifier stopped at its module lookup error before
independent relation and scalar replay. Cold setup, index construction,
rank work, resources, and per-point statuses remain in the raw archives
listed in `PROVIDER_ARCHIVES.md`.
