# N83 intermediate factor-base holdout: admitted stratum 2

K3072/root11 had the shortest target-online geometric mean among the four
frozen factor bases on the five sealed public points assigned to stratum 2:
3.548589 ms. The geometric ratio of five separately solved same-point rho
intervals to K3072 IC intervals was 540.375508. The cloud producer archived
the original verifier lookup failure; the corrected verifier then independently
replayed the raw ZIP locally, exited 0 with empty stderr, and generated an
analysis with `HALF_ADMITTED`, `PASS_TIER1`,
`PASS_INDEPENDENT_REPLAY`, and `PASS_10_PERCENT_PER_PAIR`. This is one
admitted five-point stratum; the frozen 20-point selector awaits the other
three strata and complete provider binding.

| Factor-base policy | Usable points | Folded columns | Geometric target-online ms | Mean cold process ms | K3072 faster points | Geometric paired policy/K3072 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| K2304/root11 | 382,464 | 2,304 | 26.554219 | 187,307.782 | 5/5 | 7.483037 |
| K2543/root11 | 422,138 | 2,543 | 10.516688 | 224,232.060 | 3/5 | 2.963625 |
| K2816/root11 | 467,456 | 2,816 | 29.645094 | 268,788.587 | 4/5 | 8.354052 |
| K3072/root11 | 509,952 | 3,072 | 3.548589 | 323,583.905 | — | 1.000000 |

| Point index | K2304 online ms | K2543 online ms | K2816 online ms | K3072 online ms | Same-point rho for K3072 ms |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 10 | 49.913677 | 17.284025 | 13.016262 | 22.747453 | 1,968.763954 |
| 11 | 7.119442 | 32.593728 | 14.624759 | 4.732233 | 1,640.297401 |
| 12 | 67.984941 | 13.941329 | 158.563682 | 3.332936 | 1,539.641877 |
| 13 | 27.076360 | 1.120300 | 41.495168 | 3.950998 | 1,555.783833 |
| 14 | 20.183624 | 14.621014 | 18.280494 | 0.396959 | 3,351.748318 |

The [analysis JSON](stratum2_analysis.json) records every point's verified
relation and scalar replay, calibration, host topology, cgroups, event order,
exclusive target phases, cold stages, and resources. The [provider receipt](stratum2_provider_receipt.json)
binds artifact `12107554011`; its 214,980,278-byte ZIP has SHA-256
`96a5ac837e4e7491bb47af097024d61a6e31c93059834acd287bc306675e07b2`.
The [local replay receipt](stratum2_local_replay_receipt.json) binds the
corrected verifier SHA-256
`cf5af2d981663b871370535ea7dfb4a6af305adfd45eb253d3d5ee7d5fb06557`
and canonical analysis SHA-256
`f628ecc990d58ca7b09d3d41b46653ec25ddb5491635b03aba11b4190af98b2f`.
The full raw archive remains at
`/Volumes/SSD990/crypto/artifacts/n83_intermediate_tier1_20261017_raw`.
