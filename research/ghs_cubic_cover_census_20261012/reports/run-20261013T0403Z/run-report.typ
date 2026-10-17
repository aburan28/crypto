#set page(paper: "a4", margin: 1.8cm, numbering: "1")
#set text(size: 9.5pt)
#show table: set text(size: 8pt)
= Run report: cover census, run 20261013T0403Z

Descriptive record of one run. Interpretation is in the study `README.md` / `REPORT.pdf`. This run re-executed the same deterministic tool after it gained aggregate `mod4` records and a neighbor-model count; its per-curve records are identical to run `20261013T0305Z` (checked by diff of the `curve` records).

== What ran
- Command: `rustc --edition 2021 -O cover_census.rs && ./cover_census 5 && ./cover_census 7 && ./cover_census 11 120`, 8 threads.
- Spec: `PROTOCOL.md` (sha256 in the environment block below), protocol-versioned at commit, not independently preregistered. Fixed seed 20261012 for the q = 11 sample.
- Isolation: none; shared host under heavy load. No timing is a result. All recorded values are exact counts and deterministic.
- Raw records: `results/run-20261013T0403Z.jsonl` (`rec` = field, curve, summary, mod4, totals).

#raw(read("../../results/run-20261013T0403Z.env.txt"), block: true)

== Runs requested and completed
#table(columns: 4, stroke: 0.4pt,
  [*Arm*], [*Requested*], [*Completed*], [*Notes*],
  [q = 5], [all classes, all structures], [252 curves x 16275 structures], [],
  [q = 7], [all classes, all structures], [690 curves x 120099 structures], [],
  [q = 11], [120 sampled curves and their 2-isogeny neighbours, all structures], [120 primary + 127 neighbor models = 247 models, each over 1786323 structures], [the earlier run's first attempt failed on a shell quoting error with no output; its command and stderr were not kept (see environment block)],
)

== Recorded values
Unit: number of isomorphism classes of curves (q = 5, 7) or sampled curves (q = 11). Counts are exact; no uncertainty is attached. `mismatches` (genus formula against the sum over quadratic subcovers) is 0 in all three arms. The plotted fractions and the n under each group are the tool's `mod4` records, not computed by the plotting script.

#figure(image("genus3_fraction.svg", width: 100%), caption: [Run 20261013T0403Z. Fraction of curves whose minimal cover genus is 3, by \#E(K) mod 4, as recorded in the `mod4` records (`genus3_fraction`); n under each group is the recorded curve count. q = 11 is a sample of 120.])

#figure(image("../branch_types.svg", width: 100%), caption: [Schematic of the branch configurations behind the measured (m, r) types.])

== Missing data
- Isogeny-class minimum at q = 11 (sample only).
- Hyperellipticity of the m = 2 genus-3 covers: not decided by branch data.
