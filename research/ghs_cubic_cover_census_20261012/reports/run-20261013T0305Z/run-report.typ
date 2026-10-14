#set page(paper: "a4", margin: 1.8cm, numbering: "1")
#set text(size: 9.5pt)
#show table: set text(size: 8pt)
= Run report: cover census, run 20261013T0305Z

Descriptive record of one run. Interpretation is in the study `README.md` / `REPORT.pdf`.

== What ran
- Command: `rustc -O cover_census.rs && ./cover_census 5 && ./cover_census 7 && ./cover_census 11 120`, 8 threads.
- Spec: `PROTOCOL.md`, frozen before the run. Fixed seed 20261012 for the q = 11 sample.
- Isolation: none; shared host under heavy load. No timing is a result. All recorded values are exact counts and deterministic.
- Raw records: `results/run-20261013T0305Z.jsonl` (`rec` = field, curve, summary, totals).

#raw(read("../../results/run-20261013T0305Z.env.txt"), block: true)

== Runs requested and completed
#table(columns: 4, stroke: 0.4pt,
  [*Arm*], [*Requested*], [*Completed*], [*Notes*],
  [q = 5], [all classes, all structures], [252 curves x 16275 structures], [],
  [q = 7], [all classes, all structures], [690 curves x 120099 structures], [],
  [q = 11], [120 sampled curves and their 2-isogeny neighbours, all structures], [120 curves x 1786323 structures], [first attempt failed on a shell quoting error and produced no output; rerun with the same binary],
)

== Recorded values
Unit: number of isomorphism classes of curves (q = 5, 7) or sampled curves (q = 11). Counts are exact; no uncertainty is attached. `mismatches` (genus formula against the sum over quadratic subcovers) is 0 in all three arms.

#figure(image("genus3_fraction.svg", width: 100%), caption: [Run 20261013T0305Z. Fraction of curves whose minimal cover genus is 3, by \#E(K) mod 4; n under each group is the number of curves. Derived as a ratio of two recorded counts. q = 11 is a sample of 120.])

#figure(image("../branch_types.svg", width: 100%), caption: [Schematic of the branch configurations behind the measured (m, r) types.])

== Missing data
- Isogeny-class minimum at q = 11 (sample only).
- Hyperellipticity of the m = 2 genus-3 covers: not decided by branch data.
