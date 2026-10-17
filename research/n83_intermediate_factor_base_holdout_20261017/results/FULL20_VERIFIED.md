# N83 root11 20-point Tier 1 result

**K2304/root11 is the selected complete cold-process factor base among the
four frozen policies.** Its mean cold process was 231.858 seconds across 20
fresh public points, 44.396 seconds shorter than K2543, 99.024 seconds
shorter than K2816, and 159.663 seconds shorter than K3072 on paired points.
Each of those three two-sided exact sign-flip p values is
`0.0000019073486328125`, below the preregistered Bonferroni threshold
`0.008333333333333333`.

The primary target-online comparison is `BOUNDED_UNKNOWN`: K3072/root11 has
the shortest observed geometric online interval, 6.206 ms, but its comparison
with K2816 has two-sided exact sign-flip p `0.016202926635742188`, above the
six-comparison threshold. The fresh 32-point K2816/K3072 head-to-head in
[PR #2102](https://github.com/aburan28/crypto/pull/2102) tests that remaining
pair with a predeclared two-arm rule.

| Root11 factor base | Actual usable points | Folded columns | Geometric one-target online ms | Mean complete cold process s | Geometric same-point rho/IC online |
| --- | ---: | ---: | ---: | ---: | ---: |
| K2304 | 382,464 | 2,304 | 30.605300 | **231.858** | 49.945708 |
| K2543 | 422,138 | 2,543 | 17.279708 | 276.254 | 143.728003 |
| K2816 | 467,456 | 2,816 | 18.418924 | 330.882 | 76.196215 |
| K3072 | 509,952 | 3,072 | **6.205612** | 391.521 | 308.309304 |

All 80 IC target results and their separately seeded same-point rho references
passed independent scalar and relation replay; each rho online interval was
longer than its paired IC interval. Each IC cold process includes factor-base
and index construction, queries, decomposition and scalar checks. The primary
online interval starts with target-dependent work after that reusable setup.
The [selection JSON](panel_39042774613/selection.json) retains every exact
pairwise p value, effect interval and per-host rho/IC ratio; the
[validated rows](panel_39042774613/validated_rows.json) retain target phases,
resource receipts, base hashes and correctness checks.

Four distinct GitHub Actions runner IDs and job IDs bound to the provider
snapshot; their execution windows overlapped from `15:40:36` to `16:56:07`
UTC. Every five-point stratum passed `PASS_TIER1`,
`PASS_10_PERCENT_PER_PAIR`, `PASS_INDEPENDENT_REPLAY`, resource and event-order
gates. The original cloud verifier's module lookup error is retained in each
raw ZIP; a corrected verifier independently produced each analysis from the
archived result members. All four local replays exited 0 with empty stderr.
The cloud jobs accordingly retain their original `failure` conclusions while
the aggregate records the authenticated local replay and distinct provider
runner identities. The [final provider jobs](provider_jobs_final.json) snapshot
has SHA-256
`22e3bdb66f49405e954a9fb78d845b3f8d31549b7c111b98f4ab4d7cceeb7fc2`.

| Stratum | Points | Provider artifact | ZIP SHA-256 | Local analysis SHA-256 |
| ---: | --- | ---: | --- | --- |
| 0 | 0–4 | 12109216834 | `ad4426d5ca3b9a275a2865c92ace3d021b0e945d9d181fb3b2fe4137af404f0b` | `fd32303e6be10b20bfaf376fc3b819479b9b9e7059c83deb518ccc0d246eddb0` |
| 1 | 5–9 | 12108212991 | `4bbed43eec58085be04fbc21fa46e0a1745230419f6e09015f317b120b6dc83c` | `b566c5f186cdfdc10d00214373142de21898c35aacdb09cd56a68055995c2591` |
| 2 | 10–14 | 12107554011 | `96a5ac837e4e7491bb47af097024d61a6e31c93059834acd287bc306675e07b2` | `f628ecc990d58ca7b09d3d41b46653ec25ddb5491635b03aba11b4190af98b2f` |
| 3 | 15–19 | 12108462941 | `090651938efa02504b0f9014d841949f5b738c3ab626160d0032b6d626bf6f88` | `58c2ada0aa86d2fcc28ce89ddfa082245f49dfddbd4fd59da6ad421830443db9` |

Each stratum has versioned `stratumN_analysis.json`,
`stratumN_provider_receipt.json`, and `stratumN_local_replay_receipt.json` beside
this report. The approximately 820 MiB of raw producer ZIPs remain under
`/Volumes/SSD990/crypto/artifacts/n83_intermediate_tier1_20261017_raw`.
