#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9pt)
= Run report: run-20261015T0728Z (D4 cyclotomic falsifier: exact yield + ecbench cost)

*Commands.* `research/cyclotomic_falsifier_20261015/run_yield.sh` (three invocations of `examples/cyclotomic_falsifier_yield.rs`, `--nulls 1999 --seed 20261015 --family 15`), then `run_cost.sh` (`ecbench run --wait --cpus none --allow-busy`, `verify --replay-all`, `table`, `compare --save`) on `specs/cost-a1-n23.json` and `specs/cost-a0-n23.json`. Protocol frozen at commit `9972ea0e8`; amendments `2984c8ca5`, `3512a53db`. Yield start 2026-10-15T07:28Z, cost start 07:36Z (UTC). Binaries: yield sha256 `be6c9ea3…311c`, ecbench `4d8bd56d…27f6`.

*Host and isolation.* 4 vCPU shared container, five lanes running; one thread, `nice 10`, 2.5 GB address-space cap; ecbench L0 on every record. Wall times descriptive only.

*Workloads.* Yield: $K_1\/F_(2^23)$ (r = 4196903, h = 2), $K_0\/F_(2^23)$ (r = 2095853, h = 4), $K_1\/F_(2^11)$ (r = 991); exact over all r targets; 1999 null bases per rung. Cost: 16 public targets per curve (seeds 202610151001, 202610151003), 2 rounds + 1 warm-up.

*Completion.* Yield: 3/3 rungs written (349 s, 125 s, < 1 s); all self-checks passed. Cost: 240/240 records verified on each session; replay-all 160/160 identical on each; audits `ok` (receipts `62cb9d9c…cbf0`, `b7590fbe…0e4f`). The earlier session `sessions/void-a1-n23-full-unpriced` errored on every record (missing `full_prices.json` row) and was stopped (amendment A2).

*Measured values (yield, n = 23).* All 15 cyclotomic tests: effect vs null median between −0.49% and +0.65%, two-sided rank p between 0.363 and 1 (critical 0.00333). Positive control: p = 0.001 on all 12 n = 23 tests (−35% to −64% coverage, +246% to +293% energy). Full table: study README §3; raw: `yield/yield-*.json`.

#figure(image("yield-vs-null.svg", width: 100%), caption: [Exact coverage and energy relative to the null median, n = 23 cells; thin line = min–max of 1999 null bases, bar = 2.5–97.5%; seed 20261015.])

*Measured values (cost, mean S, 95% interval as `ecbench table` recorded; IC values are lower bounds).*

#table(columns: 6, inset: 4pt,
[*session*], [*curve*], [*arm*], [*mean S*], [*95%*], [*÷ rho-frob*],
[`ECBS1h6e7055fb8555`], [K1 n23], [rho-frob], [0.439], [[0.375, 0.512]], [1],
[], [], [ic-cyc-w1], [4.562], [[3.021, 6.754]], [10.40],
[], [], [ic-cos-c1], [5.162], [[3.825, 6.637]], [11.76],
[], [], [ic-cyc-w2], [13.154], [[13.146, 13.161]], [29.98],
[], [], [ic-cos-c30], [12.872], [[12.864, 12.880]], [29.34],
[`ECBS1h6c2295d00f35`], [K0 n23], [rho-frob], [0.385], [[0.323, 0.456]], [1],
[], [], [ic-cyc-w2], [7.042], [[7.032, 7.054]], [18.29],
[], [], [ic-cos-c17], [6.276], [[6.264, 6.288]], [16.30],
[], [], [ic-cyc-w3-c3], [1.028], [[0.913, 1.147]], [2.67],
[], [], [ic-cos-c3], [1.052], [[0.937, 1.183]], [2.73],
)

#figure(image("cost-s-by-arm.svg", width: 100%), caption: [Mean S per arm, sessions as in the table; bars = recorded 95% interval.])

*Raw records.* `session-archives/<curve>.tar.gz` (the full session directories; `records.jsonl` inside), `audits/audit-<curve>.json`, `tables/`. Interpretation is in the study README and REPORT.pdf.
