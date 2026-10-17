#set page(paper: "a4", margin: 1.7cm)
#set text(size: 9.3pt)
#set table(stroke: 0.4pt)
#show table: set text(size: 8.3pt)
= Run report: Part A ladder and ablation (E14), 2026-10-15T05:30Z

Command `RAYON_NUM_THREADS=3 research/gaudry_quartic_d4_char2_20261015/run_ladder.sh` (registered with `PROTOCOL.md`, commit `666438375`). Code `5a5828371`; `examples/gaudry_quartic_bin` sha256 `25a20a67…`. Host: Xeon \@ 2.10 GHz, 4 vCPU, shared with session B1. 22 configurations requested, 22 completed (05:30–05:45); counts are deterministic.

Arms: plain (degree 14), weighted (2,2,2,1; degree 26), Frobenius (2,2,1,1; degree 11), each with/without sparsest pivots and replay, at $m = 8$ (10 planted targets, seed 1); Frobenius + sparsest + replay at $m = 5, 6, 7, 8, 10, 12$ (planted 15–20 and uniform 50 targets, MITM check), and $m = 15$ cost only.

Measured (mean $C_4$, $FF_q$ multiplications): plain+s+r $9.54 dot 10^7$; weighted+s+r $3.60 dot 10^6$; Frobenius+s+r $1.98 dot 10^5$. Derived paired ratios (`raw/summary.json`, 95 % bootstrap, 9–10 pairs): plain/weighted 26.5 (26.4–26.6), weighted/Frobenius 18.1 (18.08–18.18), plain/Frobenius 481 (479–483); sparsest 1.33; replay 1.0004. Ladder means $1.13 dot 10^5$ ($m = 5$) to $2.50 dot 10^5$ ($m = 15$).

Failures recorded: 512 closed solves, 0 disagreements with MITM. Unsolved: 5 planted at $m = 5$, 3 planted at $m = 6$ (every decomposition of each contains a pair $\{t, t+1\}$); 0 uniform. Raw records: `raw/*.jsonl`, `raw/*.log`.

#figure(image("c4_ablation.svg", width: 90%), caption: [$C_4$ by option at $m = 8$, 10 planted targets each.])
#figure(image("c4_ladder.svg", width: 70%), caption: [$C_4$ against $m$ (mean, 95 % bootstrap); 15–50 targets per point.])
