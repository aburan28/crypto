#set page(paper: "a4", margin: 2cm)
#set text(size: 10pt)
#set par(justify: true)

= Run report: XL vs MutantXL closure on committed descent draws

*Study* `research/dreg_closure_mutantxl_20261016/` (pre-registered in
`PREREGISTRATION.md`, commit `4966c1190`). *Run* 2026-10-16T02:35Z to
02:39Z UTC (scored run); soundness control 2026-10-16 about 02:42 to 03:21 UTC. *Class* stage
diagnostic: no `S`, no rho ratio, no speed claim.

== What ran

- Command: `research/dreg_closure_mutantxl_20261016/run.sh BINARY`, which
  runs `examples/dreg_closure --max-cols 70000` once per input file. The
  input is the 41 committed run files of studies 2026-09-23 to 2026-10-06.
  Binary SHA-256 `f1734a17…a0011` (`runs/binary.sha256`).
- Systems: the chained $S_3$ decomposition systems (b = 1) of the $E_0$
  descent studies. For m = 3 there are 14 cells ($N$ = 11–25); for m = 4
  there are 9 cells ($N$ = 12–32). There are four unsatisfiable draws per
  cell, 92 in total, each rebuilt from its committed subspace basis and
  $x(R)$.
- Arms, at each degree $D$ from 3 to the committed XL degree, with
  $C(N, <= D) <= 70000$:
  - plain XL;
  - the MutantXL fall closure (`ReducedBasis::from_system_closed_full`),
    saturated.
- Host: the 4-core, 16 GB container, which is not isolated (`AGENTS.md`
  §10). Wall times are practicality notes only.

== Completion

- 41 of 41 files exit 0, and no file hit the 3 h limit (`runs/queue.log`).
- 92 of 92 draws were measured.
- The XL-identity assertion held on every degree built: the plain basis
  refutes exactly at the committed degree.
- Every closure saturated.
- Peak memory 1.0 GB; slowest single degree 39 s.

== Measured values

#figure(image("closure_vs_xl.svg", width: 100%),
  caption: [Recorded values from `runs/*.jsonl` (92 draws, run
  2026-10-16T02:35Z). A: committed XL refutation degree (triangles; the
  tri-star marker is a $>=$ bound) against the measured closure refutation
  degree (dots), four draws per cell. B: the recorded columns and counted
  64-bit word ops of the closure build at its refuting degree, against the
  unknowns $N$. Source: `plot.py` (layout only).])

- The closure refutation degree is 3 or 4 on every draw.
- The committed XL degree is 5, 6, 7 or $>= 7$.
- Per-cell values, gaps and the registered verdicts are in `score-output.txt`
  (native `dreg_score closure`):
  - Q1 below in 20 of 20 testable cells, and in the 3 exploratory cells;
  - Q2: 26 draws with gap 1 and 58 with gap $>= 2$;
  - Q3: m = 3 medians 4, 4, 4, 4 over $ell$ = 2–5.

== Control

Soundness, added after the scored run and disclosed as such: every committed
satisfiable draw in the same files, at degrees 3–5 within the cap.
0 of 345 satisfiable draws refuted by XL or the closure, all built to degree 5 (one above every closure refutation). Run in two invocations: the first was killed by the session tool's background time limit after 315 complete draws; the remaining 30 ran to exit 0 (`satisfiable-control-2`, peak 12.3 GB). The kill is a tool limit, not evidence. Records: `runs/satisfiable-control{,-2}.{jsonl,log}`; binary
SHA-256 `dcffc6b6…cc094`.

== Raw records

`runs/*.jsonl` (one line per draw, with every degree's XL and closure ranks,
columns, products, word ops and seconds), `runs/*.log`, `runs/queue.log`.
The interpretation is in `../../RESULTS.md`.
