#set page(paper: "a4", margin: 1.6cm)
#set text(size: 9.2pt)
#set table(stroke: 0.4pt)
#show table: set text(size: 8pt)

= Run report: Joux–Vitse $n = 5$ ladder and complete DLPs, 2026-10-15T05:18Z

Descriptive record of the runs (observations only; interpretation is in `../../README.md` and `../../REPORT.pdf`).

*Spec.* `PROTOCOL.md` (frozen, commit `ffccaf115`), orchestration `run.sh {ladder,dlp11,dlp31}`, binary `examples/joux_vitse_ladder.rs` at commit `1e8592da0` (rustfmt only after the freeze). Host: Claude Code cloud container, Intel Xeon \@ 2.10 GHz, 4 vCPU, 15 GB, rustc 1.97.0; isolation none; the three orchestrations ran concurrently on the 4 vCPUs.

*Timeline.* First launch 04:33Z, killed by a container restart before 05:18Z (partial records in `data_killed_20261015T0518Z/`, not used). Relaunched 05:18Z. DLP `p = 11` seeds 1 and 2 completed 07:33Z and 11:02Z; seed 3 stopped by hand at 12:14Z with 0 relations (its curve has a one-point factor base; `data/dlp/p11_s3.STOPPED.md`). DLP `p = 31` completed 12:09Z (inside the 10-hour limit). Ladder completed 15:07Z.

*Requested and completed.* Ladder: 12 primes × 2 seeds × (32 F4 trials + `500 p` meet-in-the-middle trials), all 24 files complete. DLP: 4 arms requested, 3 completed and verified, 1 stopped. Added after the freeze (addendum): one re-run of the committed `jv_quintic` plain test at `p = 271`, completed.

*Measured (unit: `F_p` multiplications; counts, not wall time).*

#table(columns: 5, align: (left, right, right, right, right),
  [quantity], [value], [n], [spread], [record],
  [trial cost, all rungs], [$5.9$–$9.2 dot 10^9$ (rung means)], [768 trials], [min $4.6 dot 10^9$, max $3.3 dot 10^10$], [`data/ladder/*.jsonl`],
  [trials with a decomposition], [491 of 1,422,768], [—], [per-rung Wilson 95 % in `summary.json`], [`data/ladder/`],
  [F4 vs oracle disagreements], [0], [768], [—], [`f4_trial` rows, `agree`],
  [relations failing to re-add], [0], [all], [—], [`relations_verified`],
  [DLP `p = 11` s1: IC online / rho], [$1.74 dot 10^13$ / $3.84 dot 10^5$], [1 target], [—], [`data/dlp/p11_s1.jsonl`],
  [DLP `p = 11` s2: IC online / rho], [$2.67 dot 10^13$ / $1.52 dot 10^5$], [1 target], [—], [`data/dlp/p11_s2.jsonl`],
  [DLP `p = 31` s1: IC online / rho], [$5.31 dot 10^13$ / $7.78 dot 10^5$], [1 target], [—], [`data/dlp/p31_s1.jsonl`],
  [`jv_quintic` plain test, today's engine], [$1.104 dot 10^10$], [3 random + 2 planted], [—], [`data/jv_quintic_rerun/`],
)

#figure(image("trial_cost.svg", width: 80%), caption: [Cost of one trial against $p$ (mean, min–max; 64 trials per rung). Plotted from `summary.json` by `../../figures/plot.py`.])
#figure(image("trials_per_relation.svg", width: 80%), caption: [Trials per decomposition against $p$ (Wilson 95 %), from the same records.])
