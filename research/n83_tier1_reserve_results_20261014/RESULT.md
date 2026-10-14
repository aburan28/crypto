# N83 K2940: five fresh same-point Tier 1 IC/rho pairs

The selected 488,040-point N83 factor base recovered five separately frozen
public points with a **210.25× median paired online rho/IC ratio** on a
four-CPU cloud instance. All five pairs passed independent scalar and
four-point-relation replay and the Tier 1 noise gate. The target-block
bootstrap 95% interval for the median of five paired ratios is
**[14.60×, 295.17×]**; the individual ratios remain visible below.

The candidate is
`IC1N83Ckb1a1fb488040PDP4rootRCguidedLAgaussTDdirectISO0h751130834947`:
`K=2940` folded columns and full rank, with 488,040 distinct usable points
before sign/Frobenius folding. Its curve is
`EC1N83Ckb1a1haacb69ee9f37`, over the polynomial representation of
`GF(2^83)` with modulus integer `9671406556917033397649557`, model
`y^2+xy=x^3+x^2+1`, subgroup order `8569786107849059`, and generator
`(6614244651398184860292015, 8868260410029249391539529)`.
The [candidate manifest](../n83_factor_base_100_results_20261014/candidate_K2940_root22.json)
fixes the exact base digest, quotient rule, decomposition, source and stage
configuration.

| Fresh public point | IC/rho order | IC online ms | Rho online ms | Rho / IC | IC setup s | IC peak GiB |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 1 | IC, rho | 16.968549 | 3567.668962 | 210.251858 | 521.811 | 12.718 |
| 2 | rho, IC | 16.519472 | 1292.000018 | 78.210733 | 479.534 | 12.717 |
| 3 | IC, rho | 18.677175 | 5513.018825 | 295.174127 | 477.263 | 12.718 |
| 4 | rho, IC | 89.561193 | 1307.512005 | 14.599091 | 470.477 | 12.719 |
| 5 | IC, rho | 11.600030 | 2965.326394 | 255.630925 | 466.600 | 12.718 |

Each row has its own canonical one-target `workload_id` and `run_id` in
[workloads.json](workloads.json) and [measurements.json](measurements.json).
The five public points were frozen before this run from the SHA-256 domain
`n83-tier1-reserve5-20261014-v1`, using
`1 + little_endian(SHA256(domain || u32_le(index))) mod (r-1)` as the
fixture scalar rule; only points were passed to the producers. Their point
file has SHA-256
`644c56a7bae2c92085704bb6f58ab70b38b278710de90152ece7631058c41dc9`.
Every pair used one fresh IC process with an empty index and log cache and
one strong-rho process on that same point. The IC clock began with
target-dependent query work after reusable setup and ended after recovery
verification; rho used its matching target-dependent interval. The
exclusive IC target phases and process costs are in the run rows and raw
receipts. Setup ranged from 466.600 to 521.811 seconds and is reported
separately from the online interval.

The Tier 1 runner exposed an Intel Xeon Platinum 8370C, four selected CPUs,
and 16,372,432 KiB host-reported memory. Affinity succeeded. All ten child
processes recorded zero cgroup throttle and OOM events; the measured steal
time on selected CPUs was zero. The largest before/after calibration median
divided by the smallest was 1.0126, inside the frozen 1.15 gate. The
five-pair artifact was produced by [GitHub Actions run 38590234551](https://github.com/aburan28/crypto/actions/runs/38590234551)
at source commit `552e379748cd0befe082ca66b2560e7316769b61`. GitHub artifact
`11913299459` has SHA-256
`00901da76e71ee5e3a70f645b4c8dbb54befe68b75afcb3b342feec489a60764`;
its `results/analysis.json` has SHA-256
`ece8b00a4b2e9a332c5039aa9bad3b40c0e23f3c161173f2b59b11d2c996e84b`.
An independent local replay of the extracted artifact returned
`PASS_TIER1_PAIRED_SPEEDUP` for all five pairs and reproduced that exact
analysis digest. It rechecked the complete base, group relation, target
scalar, exclusive phase sum, same-point pairing, and noise counters.

The 1,000-point K refinement and four-byte fixed30 layout comparison use
different secondary workloads to select the next factor base under the same
14-GiB observed RSS guard. A changed selection proceeds to its separately
frozen fresh-point Tier 1 pair run.
