#set page(paper: "a4", margin: 2cm)
#set text(size: 10pt)
#set par(justify: true)

= Run report: MutantXL closure at degree 4, N = 26–42

*Study* `research/dreg_closure_scale_20261016/`, pre-registered in
`PREREGISTRATION.md` (commit `f10c86d38`). *Run* 2026-10-16, 17:05 to
18:09 UTC. *Class:* stage diagnostic, with no `S`, no rho ratio and no
speed claim.

== What ran

- *Draws.* `run.sh LADDER CLOSURE gen` generated the draws with seed
  20261016 (exact solution count; XL built to degree 3 only), and they
  were committed before any measurement.
- *Measurement.* `run.sh … measure` ran one process per unsatisfiable
  draw, as
  `prlimit --as=15e9 timeout 7200 dreg_closure --to 4 --max-cols 130000 --draw K`.
  `run.sh … control` ran the soundness control.
- *Binaries.* SHA-256 values are in `runs/binary.sha256`.
- *Cells.*
  - Committed m = 3 ℓ = 5 draws: (11,5), (13,5) and (15,5), N = 26–30, XL ≥ 7.
  - Fresh m = 3 draws: (16,6), (18,6), (19,7) and (21,7).
  - Fresh m = 4 draws: (12,3), (14,3) and (13,4).
  - Four unsatisfiable draws per cell, 40 in total.
- *Host.* The 4-core, 16 GB container, which is not isolated (`AGENTS.md`
  §10). Wall times are practicality notes only.

== Completion

- *Processes.* 41 of 41 (40 draws plus 1 control) exited 0. None was
  killed by the memory or time limit.
- *Degree 4.* It was built in full on every draw, and every closure
  saturated.
- *Peak RSS.* 7.1 GB at (13,4) m = 4. The slowest draw took 283 s.

== Measured values

#figure(image("closure_scale.svg", width: 100%),
  caption: [Recorded values from `runs/*.d*.jsonl` (40 draws, run
  2026-10-16T17:05Z). A: the closure refutation degree per draw against
  the unknowns $N$. Triangles mark draws whose degree-4 closure saturated
  without refuting. B: the counted closure word operations against the
  columns $C(N, <= 4)$ at the refuting degree. Source: `plot.py` (layout
  only).])

- *Closure degree.* 4 on all 32 draws of eight cells.
- *Above 4.* All 8 draws of m = 3 ℓ = 7 (N = 40, 42).
- *XL.* It refuted at degree 4 on none of the 40 draws.
- *Work ratio.* log₂(word ops) ÷ log₂(columns) was 2.18–2.33 on the
  refuting draws.
- *Per-draw values.* In `score-output.txt` (native `dreg_score scale`).

== Control

Soundness, as registered: the first two satisfiable draws of each fresh
cell with $N <= 36$. Only one exists, `m4-12-3` draw 0, which has 24
solutions. Built to degree 4, it was refuted by neither XL nor the closure.
Records: `runs/control-m4-12-3.d0.{jsonl,log}`.

== Raw records

- `runs/*.d*.jsonl`: one line per draw, with every degree's XL and closure
  ranks, columns, products, word ops and seconds.
- `runs/*.log`: peak VmHWM per process.
- `runs/queue.log`; `draws/*.jsonl` (generator output).
- The interpretation is in `../../RESULTS.md`.
