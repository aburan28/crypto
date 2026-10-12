#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
= Run report: ic_semiregular_law_20261005, round 1

*Observations only.* Interpretation is in `RESULTS.md` of the same study.

== What ran
- *Command:* `lanes.sh` (study directory), which runs `examples/macaulay_dense.rs` once per (cell, draw, arm, degree D), scanning D = 1, 2, … to the registered cap d_reg + 2 and stopping at the first refuting D.
- *Registration:* commit `1c1092005` (PREREGISTRATION.md, predictions.jsonl, dump.sha256), before any registered draw ran.
- *Systems:* exported by `examples/rr_degree_ladder.rs --dump-dir`, seed 20261005, four rootless draws per cell; hashes in `runs/registered/dump.sha256`.
- *Hosts:* this container (2 lanes × 2 threads, 4 cores, 16 GB), 2026-10-05 16:58Z to 2026-10-11 23:27Z, with one restart; three remote 4-core hosts (lanes c, d, e). Isolation: none beyond one process per lane; degrees are deterministic, wall times advisory.
- *Budget per process:* 86,400 CPU-s, 13.5 GB address space; the engine refuses a basis above 12 GB.

== Runs requested and completed
- 68 units registered (17 cell-arm pairs × 4 draws); 68 resolved, none censored, none killed or refused.
- Remote lanes a, b and f stopped at the account's weekly usage limit before measuring past D = 8. Their six units were re-run from D = 1 locally. Lane b's partial log is kept in `runs/registered-remote-b` and not read.
- Excluded by registration: cell K1n61l4 (seen in a smoke test); `x4` on K1n23l6 (not run, as registered).

== Measured refutation degrees (one reading per draw)
`triv` marks a draw whose system contains the constant 1 (refutes at D = 1). Window = registered `[d_reg, d_reg + 1]`; fit = the fixed-ℓ value.

#table(columns: 6, align: left, [*cell*], [*arm*], [*readings*], [*window*], [*fit*], [*in window*],
[K0n61l4], [rr], [5 5 5 5], [(5, 6)], [6], [yes],
[K0n61l4], [x4], [8 8 8 8], [(8, 9)], [9], [yes],
[K0n61l6], [rr], [6 6 6 6], [(6, 7)], [8], [yes],
[K0n61l6], [x4], [10 10 10 10], [(9, 10)], [11], [yes],
[K1n23l5], [rr], [7 7 triv 7], [(7, 8)], [8], [yes],
[K1n23l5], [x4], [10 10 10 10], [(9, 10)], [10], [yes],
[K1n23l6], [rr], [8 8 8 8], [(7, 8)], [8], [yes],
[K1n29l5], [rr], [6 6 6 7], [(6, 7)], [8], [yes],
[K1n29l5], [x4], [10 10 10 10], [(9, 10)], [10], [yes],
[K1n45l5], [rr], [6 6 6 6], [(6, 7)], [8], [yes],
[K1n45l5], [x4], [9 9 9 9], [(9, 10)], [10], [yes],
[K1n59l6], [rr], [6 6 6 6], [(7, 8)], [8], [below],
[K1n59l6], [x4], [10 10 10 10], [(10, 11)], [11], [yes],
[K1n61l5], [rr], [5 5 5 5], [(6, 7)], [8], [below],
[K1n61l5], [x4], [9 9 9 9], [(9, 10)], [10], [yes],
[K1n61l6], [rr], [6 6 6 6], [(6, 7)], [8], [yes],
[K1n61l6], [x4], [10 10 10 10], [(9, 10)], [11], [yes])

Resolved draws: 68. Inside the window: 60. Below: 8. Above: 0.

#figure(image("degrees.svg", width: 100%), caption: [Recorded refutation degree for every resolved draw (dots), one column per cell, with each cell's registered window (shaded band) and the fixed-ℓ fit (dashed). Run ic_semiregular_law_20261005 round 1; 4 draws per cell; values from `runs/registered*/K*.jsonl`. Trace-constant draws and the excluded cell are omitted.])

== Raw records
- `runs/registered/K*.jsonl`, `*.stderr`, `progress.txt` (local lanes)
- `runs/registered-remote-c`, `-d`, `-e` (remote lanes, with `HOST.txt`)
- `runs/registered/readout.md` (printed by `analyze.py`)
- `runs/registered/topforms.jsonl` (post-hoc diagnostic; not part of the registered readout)
