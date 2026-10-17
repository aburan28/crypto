#set page(paper: "a4", margin: 2cm)
#set text(size: 10pt)
#set par(justify: true)

= Run report: one fixed variable on the m = 3, ℓ = 7 draws

- *Study:* `research/dreg_closure_specialise_20261017/`, pre-registered in
  `PREREGISTRATION.md` (commit `9e0accfc4`).
- *Run:* 2026-10-17, 02:11 to 05:38 UTC.
- *Class:* stage diagnostic. No `S`, no rho ratio, no speed claim.

== What ran

- *Script:* `run.sh` runs `dreg_closure --to 4 --max-cols 130000 --draw K --fix …`, one process per (draw, arm), under `prlimit --as=15e9 timeout 7200`.
- *Binary:* SHA-256 `be7109a4…1006`, recorded in `runs/binary.sha256`.
- *Draws:* the 8 unsatisfiable m = 3, ℓ = 7 draws of the closure-at-scale study, (19,7) and (21,7) at N = 40 and 42.
- *Arms:*
  - A0 fixes a summand coordinate to 0 (dimensions 6, 7, 7).
  - A1 fixes the same coordinate to 1 (6 as a coset, 7, 7).
  - B fixes one coordinate per summand (6, 6, 6); this is the positive control.
  - C fixes one intermediate coordinate; this is the null control.
- *Pre-registration check:* fixing a non-occurring variable reproduced a recorded outcome and both recorded ranks exactly.
- *Host:* the 4-core, 16 GB container, not isolated.

== Completion

- 32 of 32 processes exited 0.
- The queue was resumed twice, once after a background-job limit stopped it and once after a container restart. Each time, the run in flight left an empty output and was rerun with the same verified binary (see `runs/queue.log`).
- Peak RSS was 6.5 GB.

== Measured values

#figure(image("specialise.svg", width: 100%),
  caption: [Recorded outcomes per draw and arm. Filled: the degree-4
  closure refutes. Hollow: degree 4 saturates without refuting. Without
  any fixed variable, none of these draws refutes at 4. Source: `plot.py`
  (layout only).])

- *Refuting at degree 4:* A0 on 6 of 8, A1 on 6 of 8, B on 8 of 8 and C on 6 of 8.
- *Registered controls:* B is valid. C is not, since it needed ≥ 7 of 8 failures.
- *Verdict:* uninformative for H_sum vs H_max.
- *Per-run values:* in `score-output.txt`.
