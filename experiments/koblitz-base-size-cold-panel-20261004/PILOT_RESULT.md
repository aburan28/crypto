# Frozen pilot decision: smaller compact-orbit bases cost more cold

The 30 preregistered one-target pilot pairs are complete: three fresh-process
IC/rho pairs at each of five `K` values for `n=41` and `n=53`. The native
[`koblitz_base_size_pilot`](../../examples/koblitz_base_size_pilot.rs)
analyzer reads only pilot paths, verifies input identity and successful
independent replay, and recomputes the selection from raw rows in
[`PILOT_ANALYSIS.json`](PILOT_ANALYSIS.json). All rho arms and IC processes
exited 0. **Twenty-seven** pairs recovered and independently replayed the
same public-point scalar. All three n41 `K=20` IC processes reported a
failed target despite exit 0; their replayers rejected the rows, so they
are failures and contribute no verified cold time. No timeout or OOM was
observed. These are exploratory, unisolated macOS measurements.

| `n` | `K` | Actual base points | Root-index entries | Verified pairs | Median verified IC cold, ms | Median rank PDP / index / matrix, ms | Median rank attempts / failed attempts |
| ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| 41 | 20 | 1,640 | 16,340 | 0/3 | unknown | 8,991 / 3.5 / 11.0 | 43 / 23 |
| 41 | 32 | 2,624 | 41,888 | 3/3 | 3,895.0 | 3,781.5 / 9.3 / 12.1 | 33 / 1 |
| 41 | 48 | 3,936 | 94,320 | 3/3 | 907.1 | 859.2 / 21.9 / 20.7 | 48 / 0 |
| 41 | 64 | 5,248 | 167,744 | 3/3 | 561.9 | 481.3 / 41.2 / 30.8 | 64 / 0 |
| 41 | **85** | 6,970 | 295,970 | 3/3 | **471.2** | 341.4 / 76.7 / 44.9 | 85 / 0 |
| 53 | 80 | 8,480 | 338,960 | 3/3 | 17,816.6 | 17,434.4 / 110.8 / 43.6 | 80 / 0 |
| 53 | 100 | 10,600 | 529,700 | 3/3 | 15,385.0 | 15,086.2 / 176.1 / 62.6 | 100 / 0 |
| 53 | 128 | 13,568 | 867,965 | 3/3 | 9,039.6 | 8,575.4 / 289.2 / 89.3 | 128 / 0 |
| 53 | 160 | 16,960 | 1,356,314 | 3/3 | 6,902.5 | 6,288.5 / 460.0 / 132.6 | 160 / 0 |
| 53 | **220** | 23,320 | 2,564,528 | 3/3 | **4,765.2** | 3,622.4 / 869.1 / 231.2 | 220 / 0 |

The index and rank-matrix stages did shrink as the base shrank. Rank PDP
grew more than those savings: relative to the fixed baseline, the fastest
eligible **smaller** base cost 1.193× as much at n41 and 1.449× at n53 on
the pilot point. At n41 `K=20`, median producer-reported cold cost was
9,351 ms even though no target scalar was recovered; it is a failure cost,
not a verified time. This supports only a fixed-producer, fixed-point
selection decision, not a general base-size optimum or ECC2K-130 forecast.

The protocol requires a held-out test of the fastest eligible smaller `K`
at each size even when it loses to the baseline in the pilot. That rule
selects **n41 `K=64` versus `K=85`**, and **n53 `K=160` versus `K=220`**.
This selection, all 300 raw pilot files, their SHA-256 manifest, analyzer
source/binary hash and machine-readable output are committed and pushed
before any held-out IC run. The held-out public points were already frozen
separately; no held-out IC result was used to make this choice.

The next decision is whether the smaller-base penalty persists on those
held-out points with six same-point IC/rho pairs per arm. Every missing
replay, rank or target failure will remain in the result. A formal method
comparison still needs `ecbench`, a sealed workload and host isolation;
these pilot wall times cannot be promoted to a controlled speedup.
