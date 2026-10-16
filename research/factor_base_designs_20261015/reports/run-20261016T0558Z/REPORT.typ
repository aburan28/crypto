#set document(title: "Factor-base designs: G1, m = 4", author: "fb_designs_slope")
#set page(paper: "a4", margin: (x: 1.8cm, y: 1.8cm), numbering: "1")
#set text(size: 9.5pt)
#set par(justify: true)
#show raw: set text(size: 8pt)

= Factor-base designs: G1, the m = 4 series on p#sub[1] (run-20261016T0558Z)

*Kind:* solver-axis stage diagnostic (AGENTS.md §2). Observations only: nothing here prices a pipeline, reports S or shows a speedup. *Protocol:* PREREGISTRATION.md Amendment 4 (registered in commit #raw("113f31659") before any G1 row), with Amendments 1–3 for everything it does not change.

== Run

- *Command:* #raw("target/release/examples/fb_designs_slope run --arms A0,A11 --m 4 --t-min 4 --budget 3600 --grace 300 --vmem-kb 15000000 --distinct-targets --file 'G1-{arm}-m{m}' --out-dir research/factor_base_designs_20261015/runs")
- *Commit measured:* #raw("113f31659af224c2e8eb46ebf86e5a18fb0993a3"), recorded in every row.
- *Host:* "vm", Intel Xeon @ 2.10 GHz, 4 cores, 15 GB RAM, no swap, a Claude Code cloud container. *Isolation:* one system at a time.
- *Window:* 2026-10-16 05:58Z to about 06:39Z (UTC).
- *Limits:* F4 budget 3600 s, hard kill 3900 s, #raw("ulimit -v 15000000") per child. The session also runs inside a memory cgroup that the harness does not set. Its OOM killer ended two systems at about 14 GB resident (below).
- *Seeds:* as in Amendment 2, with G2 on (#raw("--distinct-targets")). G2 redrew nothing (#raw("target_redraws") = 0 on every row). #raw("compare-rows") confirms that the t = 4 systems are the F4 systems that stopped at 300 s (#raw("runs/G1-vs-F4-crosscheck.txt")).
- *Records:* #raw("runs/G1-{A0,A11}-m4.jsonl") and #raw(".log"), #raw("runs/analysis-amendment4.txt") (native #raw("analyze --amendment4")) and #raw("runs/G1-oom-kernel-log.txt").

== Graph

#figure(image("m4_d_cols_vs_b.svg", width: 100%),
  caption: [m = 4 on p#sub[1]: D (left) and widest F4 matrix (#raw("max_cols"), log#sub[2] scale, right) against b, for A0 and A11. t = 2, 3 come from F4 (300 s budget) and t ≥ 4 from G1 (3600 s, 15 GB). Both targets are plotted (n = 2 per cell; they coincide). Crosses mark memory aborts, which have no D. Rows from #raw("runs/F4-*-m4.jsonl") and #raw("runs/G1-*-m4.jsonl"), plotted by #raw("plot.py").])

== Every system in the m = 4 series

#table(columns: 9, align: center, stroke: 0.4pt + gray,
  [*arm*], [*t*], [*source*], [*vars*], [*eqs*], [*D*], [*max_cols*], [*time (s)*], [*status*],
  [A0], [2], [F4], [10], [11], [6 / 6], [1911 / 1911], [0.3 / 0.3], [ok, ok],
  [A0], [3], [F4], [14], [15], [7 / 7], [21619 / 21619], [96.4 / 93.8], [ok, ok],
  [A0], [4], [G1], [18], [19], [– / –], [– / –], [404.0 / 406.6], [aborted, aborted],
  [A11], [2], [F4], [21], [22], [3 / 3], [34 / 34], [0.0 / 0.0], [ok, ok],
  [A11], [3], [F4], [28], [29], [5 / 5], [2159 / 2159], [0.3 / 0.3], [ok, ok],
  [A11], [4], [G1], [35], [36], [6 / 6], [88740 / 88740], [751.9 / 738.6], [ok, ok],
  [A11], [5], [G1], [42], [43], [– / –], [– / –], [73.5 / 65.9], [aborted, aborted],
)

- *A0, t = 4:* both systems aborted on an allocation failure under #raw("ulimit -v 15000000") (SIGABRT), at about 405 s. F4 had stopped the same systems at 300 s, at partial D = 7.
- *A11, t = 5:* both systems were killed with SIGKILL after 66–74 s. The kernel log shows the memory-cgroup OOM killer (#raw("CONSTRAINT_MEMCG")), at 13.97 GB resident, below the 15 GB address-space limit. The same OOM event also killed an idle shell of the session that was waiting on the run. That shell was not part of the measurement.

== Verdicts (Amendment 4, G1; §1 thresholds)

#table(columns: 5, align: (left, center, center, center, left), stroke: 0.4pt + gray,
  [*arm*], [*completed cells (t)*], [*σ*], [*Δlog#sub[2] max_cols/Δb*], [*verdict*],
  [A0], [2 (2, 3)], [1.0000], [3.500], [inconclusive (fewer than 3 cells; t = 4 aborted for memory)],
  [A11], [3 (2–4)], [1.5000], [5.675], [no middle ground at this arity (σ ≥ 0.5)],
)

- *Prediction "A0's σ ≥ 0.5 at m = 4":* not decided (2 completed cells). Amendment 4 makes no further budget or memory change without a new amendment.
- *A11:* no prediction was registered. Its m = 4 series is D = 3, 5, 6 over t = 2–4.
- *Width (secondary, descriptive):* A11's #raw("max_cols") over A0's is 0.018 at t = 2 and 0.100 at t = 3. There is no shared completed cell at t = 4.

== Observations

- At m = 4, A11 starts below A0 on both D and width: D 3 vs 6 at t = 2, and 5 vs 7 at t = 3. But its D rises by 3 over two bits, and its width grows about 5.7 bits per bit of factor base, against A0's 3.5 over its two cells.
- A11's partial D at F4's 300 s stop (6) already equalled its final D at t = 4. The 3600 s budget changed its status, not its degree.
- The binding limit at m = 4 is memory for both arms: A0 at t = 4, A11 at t = 5. On this host the effective ceiling is the session's memory cgroup (about 14 GB resident), not the #raw("ulimit") the amendment raised.

== Deviations

None in the protocol. The amendment's memory limit took effect as written, but a memory cgroup this lane did not set ended two systems below it. Amendment 4 already counts an OOM kill as an abort, and those rows record #raw("child exited: signal: 9") as their reason. The kernel log is committed beside them.
