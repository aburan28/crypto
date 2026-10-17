#set page(paper: "a4", margin: 2cm)
#set text(size: 10pt)
#set par(justify: true)

= Run report: what sets the MutantXL closure degree (ℓ, N or m·ℓ)

- *Study:* `research/dreg_closure_driver_20261016/`, pre-registered in
  `PREREGISTRATION.md` (commit `b85dd3226`).
- *Run:* 2026-10-16, 20:49 to 23:06 UTC.
- *Class:* stage diagnostic. No `S`, no rho ratio, no speed claim.

== What ran

- *Draws.* `dreg_ladder` with seed 20261017 (exact solution count; XL
  built to degree 3 only).
  - Each cell's draws were committed before it was measured.
  - The first generation pass of (24,6) was cut off by my wrapper's
    1,700 s limit and rerun. The three draws both passes produced are
    identical apart from the timing field.
- *Measurement.* One process per draw:
  `prlimit --as=15e9 timeout 7200 dreg_closure --to 4 --max-cols 155000 --draw K`.
- *Cells.* m = 3 (22,6) and (24,6) at N = 40 and 42, and m = 4 (11,5) and
  (12,5) at N = 42 and 44. Four unsatisfiable draws per cell, 16 in all.
- *Host.* The 4-core, 16 GB container, not isolated (`AGENTS.md` §10).

== Completion

- *Processes.* 20 of 20 (16 draws and 4 controls) exited 0.
- *Closure.* Every closure saturated.
- *Peak RSS.* 7.8 GB. The slowest draw took 730 s.

== Measured values

#figure(image("closure_driver.svg", width: 100%),
  caption: [Recorded outcomes of this study (black outlines) and the scale
  study, by (ℓ, N). Filled markers refute at degree 4 on every draw.
  Hollow markers saturate at degree 4 without refuting on every draw.
  Source: `plot.py` (layout only).])

- *Closure.* It refutes at degree 4 on 16 of 16 draws.
- *XL.* It refutes at degree 4 on none.
- *Work ratio.* log₂ ops ÷ log₂ cols is 2.18–2.33.
- *Controls.* Four satisfiable m = 4 draws, none refuted.
- *Per-draw values.* In `score-output.txt`.

== Raw records

- `runs/*.d*.jsonl` and `.log`.
- `runs/queue.log`.
- `draws/*.jsonl`.
- Interpretation: `../../RESULTS.md`.
