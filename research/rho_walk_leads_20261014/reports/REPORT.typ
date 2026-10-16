#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#set par(justify: true)
#show heading.where(level: 1): set text(size: 14pt)

= Pollard rho walk leads: four proposals measured

Study `research/rho_walk_leads_20261014`. Main run `run-20261014T1511Z` (commit `e531dbc19`, 94 jobs), split-sign controls `run-20261014T1513Z-split` and batch sweep `run-20261014T1513Z-batch` (commit `be3610fa1`). Modal CPU containers (16 vCPU requested, 32 threads visible, CPU model not reported), no isolation receipt; step counts are exact, ns timings descriptive. Curves $y^2 = x^3+a x+b$ over $FF_p$, $p equiv 3 space (mod 4)$, prime order $N$ by baby-step giant-step from seed 20261014; main cells $N = 53594884787 approx 2^35.6$. Every arm of a trial shares secret, table and start; every collision is solved and checked against the secret. Speedup = reference mean steps / arm mean steps, paired by trial, 95% bootstrap (4,000 resamples).

== Proposals and predictions (derived before the runs)

+ *Best of $x(P+R)$, $x(P-R)$.* Choosing the best-scoring of $k$ images gives image $Y$ expected in-degree $k u^(k-1)$ ($u$: score rank); collision time scales as $1 slash sqrt(E[lambda^2])$, so steps speedup $sqrt(k^2/(2k-1))$ = 1.155 ($k=2$), 1.512 ($k=4$). Cost: $+1M+1S$ per step (shared denominator).
+ *Lazy reduction before the DP check.* Lazy DP test: no change (merged walks share their history, so share representations). Lazy table index: walks meeting through different predecessors disagree half the time, about $1 slash sqrt(2)$.
+ *Correlated walks sharing precomputed sums.* A shared (time-dependent) index stops the walk being a function of the point; a pair table $R_i+R_j$ is one step per addition.
+ *Representative rule changing per round.* Per point: just another canonical rule. Per epoch: walks meeting in different epochs diverge.

#pagebreak()
== Results

#figure(image("speedups.svg", height: 22.5cm), caption: [Steps speedup per arm against each panel's reference (`runs/run-20261014T1511Z/summary.jsonl`, split panel from `runs/run-20261014T1513Z-split/summary.jsonl`). Blue: measured with 95% bootstrap interval; red tick: derived prediction.])

#table(columns: 4, stroke: 0.4pt, inset: 4pt,
  [*lead*], [*arm*], [*measured (95% CI)*], [*predicted*],
  [1], [best of P±R, DP method], [1.129 (1.081–1.184)], [1.155],
  [1], [best of P±R_i, P±R_j, DP method], [1.463 (1.401–1.527)], [1.512],
  [1], [sign from hash of P, DP method], [0.984 (0.941–1.028)], [1],
  [1], [random function, best of 2 / 4], [1.151 / 1.511], [1.155 / 1.512],
  [1], [best of P±R, one walk, doubling escape], [1.319 (1.277–1.362)], [1.155],
  [1], [same, shift escape], [1.170 (1.134–1.205)], [1.155],
  [1], [same doubling setup, signs split over two entries], [1.159 (1.124–1.196)], [1.155],
  [2], [DP test on lazy x], [0.998 (0.995–1.001)], [1],
  [2], [table index on lazy x], [0.733 (0.703–0.766)], [0.707 (+ DP overhead)],
  [3], [index shared across walks], [0.073 (0.070–0.076), 176/1000 censored], [loss],
  [3], [pair table, one addition], [1.027 (0.982–1.073)], [1],
  [4], [negation map, fixed rule], [1.368 (1.311–1.428)], [1.414 (1.387 with DP overhead)],
  [4], [negation map, rule per point], [1.330 (1.275–1.388)], [1.414],
  [4], [negation map, rule per epoch], [0.828 (0.792–0.869)], [loss],
)

Fruitless cycles: any walk stepping by $-R$ as well as $+R$ falls into 2-cycles $P -> P-R_i -> P$. Without an escape (`e1noesc`, $N approx 2^32$, 200 trials) the ±R arms stop on a useless repeat after a mean of 821–1,974 steps and solve 0–1 of 200; the plain walk solves 200/200 after 69,210. With $r = 1024$ the measured 2-cycle rates are 488 per $10^6$ steps for random signs (= $1 slash (2r)$), 651 for best of P±R and 1,117 for best of four; the per-point and per-epoch negation rules leave the rate at 487.

The single-walk doubling-escape excess for best of P±R (1.32) disappears under the shift escape (1.17), under the DP method (1.13), and when the signs use different table entries (1.16): it belongs to that measurement setup, not to the walk.

#figure(image("cost_vs_batch.svg", width: 100%), caption: [Batched-affine step cost, branch-free Montgomery arithmetic (selftest passed), median of 21 interleaved rotating repetitions (`runs/run-20261014T1513Z-batch/summary.jsonl`). Left: cost of best of P±R over the plain step against walks per simultaneous inversion; the dotted line is the 1.155 break-even. Right: absolute ns per step.])

Cost of best of P±R relative to a plain step, P-256: 1.215 at 1024 walks per inversion (net $1.155 slash 1.215 = 0.95$), 1.157 at 128, 1.119 at 64, 1.082 at 32, 1.050 at 16, 1.029 at 8 (net 1.12). One-limb field: 1.356 at 1024, 1.084 at 8. Best of four costs 2.43 (P-256) and 2.59 (one limb) at 1024, net about 0.6. Absolute cost favours large batches (P-256 plain 179 ns at 1024 against 938 ns for best of P±R at 16), and the Fermat inversion used here overstates the inversion share that makes small batches favourable.

== Observations

- Lead 1 gives a real per-step gain at the predicted size, but the extra $1M+1S$ outweighs it wherever large inversion batches are possible. It is a candidate only where the batch per inversion is forced small (GPU kernels such as `gpu/nist`, W = 16–32), and must then be measured in that kernel with its own inversion. Fruitless-cycle handling is required.
- Lead 2 is harmless for the DP test and harmful for the table index; it saves at most one conditional subtraction.
- Lead 3 has no mechanism that keeps the walk a function of the point; the measured variant loses more than 13×.
- Lead 4 per point changes nothing (same speedup, same 2-cycle rate); per epoch is worse than not using the negation map.
