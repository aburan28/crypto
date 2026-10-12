# Fixed n41/n53 one-target diagnostic: cold setup dominates

The frozen native compact S3-root-index producer recovered one previously
unseen public point at each of two Koblitz sizes in six fresh-process IC/rho
pairs. Every pair reached full rank, recovered the same scalar in both arms,
and passed independent general-curve replay. The target was a point generated
by public hash-to-curve and cofactor projection, not a scalar supplied to IC.
These are **exploratory, unisolated macOS timings** from example binaries,
not a formal `ecbench` cross-method claim or a controlled speedup. All paired
ratios below are descriptive of these two fixed points and this host. The
machine-readable result keeps `candidate_id` and `workload_id` null until a
complete IC1 manifest and sealed `ecbench` workload exist; the exact curve
IDs, factor-base digests, seeds and Q coordinates are retained here.

| Curve | Frozen Q | Actual base points / folded columns | Root-index entries | Recovered scalar |
| --- | --- | ---: | ---: | ---: |
| `icv1-f2m41-tm2308219-7f48b14a` (`EC1N41Ce0he09550ab560a`) | `[847899556326,1841367230650]` | 6,970 / 85 | 295,970 | 546043989315 |
| `icv1-f2m53-tm56619371-dac20a85` (`EC1N53Ce0hb097de99be9a`) | `[4024915783397303,7453281239598237]` | 23,320 / 220 | 2,564,528 | 17977960908532 |

Each cell used the exact same Q and base digest in all six pairs, with rho
first on odd repeats and IC first on even repeats. Rank policy was guided:
decompose `[a]G - R_j` for the first pivotless column `j`. At n41 all six
runs took 85 attempts for 85 verified rank rows; at n53 all six took 220
attempts for 220 rows. There were no recorded rank failures or dependent
rows. This guided rank construction is **not** a measurement of natural
ordinary-query PDP yield. Both scalars passed `[k]G = Q` and the independent
replayer checked every base-point log and every rank row.

| n | Pair | IC online ms | Rho online ms | Rho/IC online | IC cold ms | Rho cold ms | Rho/IC cold |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 41 | 1 | 3.399 | 28.462 | 8.373 | 642.389 | 28.822 | 0.0449 |
| 41 | 2 | 3.502 | 18.048 | 5.154 | 460.587 | 18.293 | 0.0397 |
| 41 | 3 | 3.361 | 18.252 | 5.431 | 456.710 | 18.494 | 0.0405 |
| 41 | 4 | 4.410 | 22.320 | 5.062 | 525.696 | 22.616 | 0.0430 |
| 41 | 5 | 3.425 | 27.544 | 8.042 | 530.801 | 27.790 | 0.0524 |
| 41 | 6 | 3.460 | 18.228 | 5.269 | 461.240 | 18.477 | 0.0401 |
| 53 | 1 | 5.167 | 110.326 | 21.352 | 4864.455 | 110.679 | 0.0228 |
| 53 | 2 | 5.047 | 100.598 | 19.933 | 4888.873 | 100.965 | 0.0207 |
| 53 | 3 | 5.041 | 118.481 | 23.506 | 4753.418 | 118.915 | 0.0250 |
| 53 | 4 | 5.269 | 106.957 | 20.301 | 4629.621 | 107.308 | 0.0232 |
| 53 | 5 | 5.593 | 98.697 | 17.646 | 5033.657 | 99.035 | 0.0197 |
| 53 | 6 | 5.553 | 100.432 | 18.088 | 5022.901 | 100.782 | 0.0201 |

The medians and observed ranges quantify run-to-run timing spread on this
one point; they do not estimate variation over fresh targets or a confidence
interval for an isolated host.

| n | IC online median (range), ms | Rho online median (range), ms | Paired Rho/IC online median (range) | IC cold median (range), ms | Rho cold median (range), ms | Paired Rho/IC cold median (range) |
| ---: | --- | --- | --- | --- | --- | --- |
| 41 | 3.442 (3.361–4.410) | 20.286 (18.048–28.462) | 5.350 (5.062–8.373) | 493.468 (456.710–642.389) | 20.555 (18.293–28.822) | 0.04176 (0.03972–0.05236) |
| 53 | 5.218 (5.041–5.593) | 103.777 (98.697–118.481) | 20.117 (17.646–23.506) | 4876.664 (4629.621–5033.657) | 104.137 (99.035–118.915) | 0.02170 (0.01967–0.02502) |

For the dashboard's IC/rho direction, the medians of the **paired** online
ratios are 0.186966 at n41 and 0.049714 at n53; the cold ratios are
23.969714 and 46.186164. These are medians of six individual IC/rho
ratios, not reciprocals of the rho/IC medians in the table.

The IC online interval begins after reusable base, root index, rank logs and
linear algebra, and contains target query, PDP, relation check, descent and
recovery check. Its five exclusive phases sum exactly in each producer row
and are checked again by the replayer. Rho online is walk plus collision
resolution and scalar-validation time after jump setup. The supplementary
in-process cold IC interval is `setup_complete_ns/10^6 + online_ms`; the
rho cold interval is `setup_ms + walk_ms + validation_ms`. The public-point
fixture generation and process launch are outside both charged intervals.
The producer's 15 exclusive cold phases sum to its cold interval, and the
replayer rejects an altered phase. Raw `/usr/bin/time` wall/CPU measurements
are archived separately. The 16 GiB RSS gate was observed with `getrusage`,
not enforced by the kernel; maxima were 45.5 MB IC / 3.5 MB rho at n41 and
337.1 MB IC / 7.6 MB rho at n53.

## Decision from the charged stages

At n41 the median rank PDP cost was **354.595 ms**, root-index build
**74.549 ms**, and rank-matrix work **51.105 ms**, against **20.555 ms**
median rho cold. At n53 those costs were **3687.618**, **912.892** and
**241.364 ms**, against **104.137 ms** rho cold. Even an impossible
zero-cost rank PDP leaves a paired observed cold remainder of at least
**123.302 ms** at n41 and **1126.221 ms** at n53, above the largest
observed rho cold times of 28.822 and 118.915 ms. Zeroing both rank PDP
and root-index build still leaves at least 52.759 and 256.516 ms,
respectively. These are **counterfactual fixed-pipeline stage bounds**, not
measured alternative algorithms or asymptotic lower bounds. The root-index
and matrix stages individually exceeded the largest rho cold observation
at each size. Optimizing only the S3 root kernel cannot make these fixed
K=85/220 cold pipelines competitive on this panel.

The next bounded experiment should vary the actual usable base size and
index architecture together on the same held-out point, preserving every
failed target, rank attempt, and full cold charge. A smaller base can reduce
both root-index construction and rank-matrix cost, but may raise target PDP
failure; that tradeoff must be measured rather than inferred from the online
numbers here. The subsequent cross-method claim needs the root-index method
inside native `ecbench`, a sealed session and operation-unit calibration,
an independent host audit, and the repository's host-isolation receipt.
This panel says nothing directly about ECC2K-130 or an isogeny descendant.

## Reproduction and failure record

The [protocol](PROTOCOL.md) and frozen [build](BUILD.md) were pushed before
either public point was generated. Run the six pairs with the
[runner](run_pair.sh), one `n`/repeat invocation at a time, and analyze them
with `cargo run --offline --locked --release --example
koblitz_cold_online_panel`; its complete machine-readable output is
[analysis.json](analysis.json). [Analyzer hashes](ANALYSIS_SHA256SUMS) pin its
source and binary. The [raw hash manifest](RAW_SHA256SUMS) covers
both point files, the preflight failure, every base/rank trace, producer row,
stderr, and independent replay receipt. The literal Q values were searched
against committed prior research/experiment files; n41 matched only its
own frozen preflight/target artifacts and n53 had no prior match.

The first n41 rho solve produced and verified the frozen Q, but the outer
`/usr/bin/time -l` wrapper exited one because macOS denied
`sysctl kern.clockrate`. Its raw files remain under
`runs/n41/preflight_time_l_failure/`; the six accepted pairs reused the
same Q after the runner switched to `/usr/bin/time` without `-l`. This was
an environment failure, not a discarded unfavorable cryptanalytic outcome.
No other pair failed, timed out, exceeded the RSS gate, or failed replay.
Three n41 negative controls were rejected by the independent replayer:
changing the recovered scalar (`target relation or recovered scalar failed
replay`), the first rank-row coefficient (`rank row 0 differs`), and one
rank-PDP cold phase by 1 ms (`IC exclusive cold phases do not match charged
intervals`). `rustfmt --check` on the touched Rust examples, `sh -n` on
the runner, and `git diff --check` passed.
