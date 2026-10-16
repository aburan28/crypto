# Pollard rho walk leads (2026-10-14)

Four proposals for cheaper Pollard rho walks on prime-field elliptic curves,
each turned into a measurement with a derived prediction, controls and
paired trials, run on Modal CPU containers. Report: [`reports/REPORT.pdf`](reports/REPORT.pdf).

| lead | proposal | what was measured | result |
|:--|:--|:--|:--|
| 1 | step to the better of x(P+R), x(P−R) (shared inversion) | steps to a useful collision, 3 setups; ns per batched step, 6 batch sizes, 2 fields | **real per-step speedup, at the predicted size** (1.13× [1.08, 1.18] under distinguished points; 1.155× predicted). Its cost (1.22× per step on P-256 with 1024 walks per inversion) makes it a net loss with large batches and a small net gain only with small batches (cost 1.03–1.12× at 8–64 walks per inversion on P-256). Any ±R walk needs fruitless-cycle handling, which was not in the proposal |
| 2 | lazy reduction before the distinguished-point (DP) check | steps to collision with the DP test or the table index on a redundant x | DP test on lazy x: no change (0.998× [0.995, 1.001]); table index on lazy x: 0.733× [0.703, 0.766], about the predicted 1/√2. Fine where it is harmless, and it saves at most one conditional subtraction per step |
| 3 | correlated walks sharing precomputed sums | an index sequence shared by all walks; a precomputed pair table R_i+R_j | shared index: 0.073× [0.070, 0.076], 176/1000 trials censored, so this is a lower bound on the loss; pair table: 1.027× [0.982, 1.073], one addition still advances one step |
| 4 | change the class-representative rule each round | negation-map walks, fixed rule / rule per point / rule per epoch | per point: 1.33× [1.28, 1.39], same as fixed 1.37× [1.31, 1.43] and the same fruitless 2-cycle rate (487 per 10^6 steps, 1/(2r) = 488); per epoch: 0.83× [0.79, 0.87], **worse than no negation map** |

Speedups are reference mean steps / arm mean steps, paired by trial seed, with
95% bootstrap intervals over trials (4,000 resamples). Each trial seeds
everything except the curve: the secret, the table and the starting points.
Every collision is turned into a discrete log and checked against the
secret.

## Lead 1 in detail

*Prediction (derived before the run).* Taking the best-scoring of k
candidate images by a score on the image gives point Y an expected in-degree
λ(Y) = k·u^(k−1), where u is Y's score rank. The collision time of a random
mapping scales as 1/√E[λ²], so the steps speedup should be
√(k²/(2k−1)): **1.155 for k = 2, 1.512 for k = 4**. This is the whole effect.
Steering towards distinguished points adds nothing, because a DP is only a
marker.

*Measured.*

- The positive control (abstract random function) hits the prediction:
  1.151 [1.114, 1.188] and 1.511 [1.461, 1.562].
- The two-table control (best of P+R_i, P+R_j) gives 1.152 [1.117, 1.188].
- Under the distinguished-point method, which is how rho is actually run,
  best of P±R gives 1.129 [1.081, 1.184] and best of four gives
  1.463 [1.401, 1.527].
- Choosing the sign from a hash of P (no selection) gives no change:
  0.984 [0.941, 1.028].

*Fruitless cycles.* Stepping by −R as well as +R creates 2-cycles
P → P−R_i → P, the negation-map problem. Without an escape, every ±R arm
repeats on a useless collision after a mean of 821–1,974 steps (the plain
walk's useful one takes 69,210 at N ≈ 2^32): 0 or 1 solved out of 200,
against 200/200 for the plain walk (`e1noesc`). With r = 1024 the
cycles occur at 651 per 10^6 steps for best of P±R and 1,117 for best of
four. Every ±R arm here escapes 2-cycles by doubling and restarts any walk
that runs 20/θ steps without a distinguished point.

*An artifact, run down.* In the single-walk, store-everything setup with a
doubling escape, best of P±R measured 1.32 [1.28, 1.36], above the
prediction. Three checks show this comes from that setup and not from the
walk:

- the shift escape gives 1.170 [1.134, 1.205];
- the distinguished-point method gives 1.129;
- splitting the signs over two table entries gives 1.159 [1.124, 1.196].

Only same-entry ± candidates combined with the doubling escape exceed the
prediction. The deployable figure is the distinguished-point one.

*Cost.* Best of P±R adds about 1M+1S per step, because the two candidates
share the denominator x_R − x_P. Measured with branch-free Montgomery
arithmetic, the step costs this much more than a plain step:

| field | walks per inversion | cost vs plain step | net vs 1.155 |
|:--|:--|:--|:--|
| P-256, 4 limbs | 1024 | 1.215 | 0.95 (a loss) |
| P-256 | 128 | 1.157 | 1.00 |
| P-256 | 64 | 1.119 | 1.03 |
| P-256 | 32 | 1.082 | 1.07 |
| P-256 | 16 | 1.050 | 1.10 |
| P-256 | 8 | 1.029 | 1.12 |
| 1 limb | 1024 | 1.356 | 0.85 |
| 1 limb | 8 | 1.084 | 1.06 |

Best of four costs 2.4–2.6× a plain step, and its net is about 0.6.

Absolute cost still favours large batches. On P-256 a plain step costs 179
ns at 1024 walks and best of P±R costs 938 ns at 16. The inversion here is
plain Fermat, which overstates the inversion share that small batches pay.

*Open lead.* GPU rho kernels batch few walks per thread, for example
`gpu/nist` at W = 16–32. A kernel test there, with the kernel's own
inversion, decides whether best of P±R is worth its 1M+1S.

## Layout

- `harness/`: the Rust harness (no dependencies). Subcommands `e1`, `rf`,
  `e2`, `e4`, `e1dp`, `bench`, `selftest` and `summary`. Every number in
  this study comes from it.
- `modal_run.py`: the launcher. It builds `harness/` in a
  `rust:1.88-slim-bookworm` image, shards trials over 16-vCPU containers and
  writes stdout unchanged beside `launch.json`.
- `plot.py`: lays out the `summary.jsonl` files as SVG/PNG.
- `runs/run-20261014T1511Z/`: the main run (94 jobs, commit `e531dbc19`).
- `runs/run-20261014T1513Z-split/`: the split-sign controls (commit
  `be3610fa1`).
- `runs/run-20261014T1513Z-batch/`: the batch-size cost sweep (commit
  `be3610fa1`).
- `runs/run-20261014T1506Z/`: the first run (commit `d2d04dc27`), raw
  records kept, **not used**. Its benchmark used data-dependent branches,
  and walks stuck in fruitless cycles timed faster. Its `e1noesc` records
  reused the trial ids of `e1` at 32 bits. It has no summary, and the run
  above replaces it.

Curves: y² = x³ + ax + b over F_p, p ≡ 3 (mod 4), prime order N found by
baby-step giant-step. They are generated from seed 20261014: bits 36 gives
p = 53595032867 and N = 53594884787 (about 2^35.6). Parameters are in each
JSONL's `curve` line.

## Reproduce

```sh
cd research/rho_walk_leads_20261014
cargo build --release --manifest-path harness/Cargo.toml
harness/target/release/rho_walk_leads selftest
pip install modal python-socks          # python-socks: Modal through an HTTP proxy
MODAL_TOKEN_ID=… MODAL_TOKEN_SECRET=… modal run modal_run.py --plan full --out runs/<tag>
cat runs/<tag>/*.jsonl | harness/target/release/rho_walk_leads summary > runs/<tag>/summary.jsonl
# committed raw records are xz-compressed: xzcat runs/<tag>/*.jsonl.xz | … summary
# (the batch summary is one `summary --in` per bench-m<walks> file, prefixed with batch_walks)
python plot.py runs/<main> runs/<split> runs/<batch> reports
```

Host: Modal CPU containers (16 vCPU requested, 32 threads visible, CPU model
not reported). No isolation receipt. The ns timings are descriptive: they are
medians of 21 interleaved, rotating repetitions on a shared host.
