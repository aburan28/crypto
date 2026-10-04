# F6-IC geometric closure: prepared panel and full small-curve control

F6-IC is a new **IC-specific hybrid solver design**: it lets exact usable
factor-base geometry refute or close branches inside an inherited Boolean F4
search. It is a project name, not Faugère's F6 or a new complexity theorem.
Each version has a distinct source commit, candidate manifest, previously
unseen public target, workload ID and immutable run receipts. All target
queries, failed PDP attempts, point arithmetic and scalar replay stay inside
the five exclusive online phases. The state supplies 62 actual usable base
points and 29 folded relation columns; those numbers are not interchangeable.

| Version and frozen result | F4 / F6 reductions across the target's attempts | F4 / F6 median online ms on that version's own target | Geometric mechanism | Decision |
| --- | ---: | ---: | --- | --- |
| [v0](pilot/V0_RESULT.md) | 16 / 16 | 3.842 / 4.114 | Two fixed x-codes, no closure reached | Negative |
| [v1](pilot_v1/V1_RESULT.md) | 379 / 163 | 47.813 / 46.246 | One fixed x-code, general point arithmetic | Reduced search; high point cost |
| [v2](pilot_v2/V2_RESULT.md) | 911 / 383 | 118.850 / 82.966 | Single-word exact point arithmetic | Better complete call; below 2x |
| [v3](pilot_v3/V3_RESULT.md) | 696 / 295 | 88.125 / 58.892 | Batched exact point arithmetic | Below 2x; stop branch |

Every row used a **different fixed public target** and a matched inherited
F4 arm, so absolute milliseconds between rows must not be compared as a
version-to-version speedup. Each row has three process repetitions; the raw
receipts retain all statuses and attempts. All 24 processes completed,
recovered the matching scalar within each workload, and passed replay.

A [matched F4/F5/F6-IC diagnostic](three_way/RESULT.md) used one new target
for all three arms. Its exploratory median complete online times were
32.134 ms for F4, 6295.054 ms for F5, and 20.906 ms for F6-IC. F6 beat
both on that target, yet missed the preregistered 2× F4 gate in every
repetition. Its corrected stage-source manifests have new candidate IDs
even though the F4/F6 executable is the same as v3.

The [eight-target panel](eight_target/RESULT.md) froze eight more public
points. All 40 process runs verified. Across targets, the exploratory
complete online F4/F6 ratio ranged **1.291–1.583×** (median **1.502×**),
and F5/F6 ranged **287.921–536.449×** (median **375.898×**). F6 won on
every target, but none met the 2× F4 gate. These are observed sample ranges,
not asymptotic or controlled CPU bounds.

The [full cold n9 control](small_cold/RESULT.md) built its own 14-point
usable base, collected and checked eight ordinary relations, solved both
folded matrix columns, and recovered a new public logarithm with F4, F5
and F6. Independent native replay checked all relations, column logs and
the target. On that one tiny curve the exploratory inside-worker cold
intervals were **4.707, 61.937 and 3.150 ms**, respectively; the online
intervals were **0.372, 3.118 and 0.272 ms**. No prepared log table was
imported in this control.

The CPU host was not isolated. All wall times are exploratory; no controlled
speedup is established. These IC-variant comparisons did not run rho. The primary
IC-versus-rho online ratio, operation-normalized `S`, general natural relation yield,
and any effect on the n53 `PDP4root` path remain **unknown**. The existing
strong-rho crossover cannot be inferred from a faster prepared n17 PDP.

The branch reached a clear decision under its preregistered gate: geometric
closure consistently reduced the prepared call, but its own group work kept
the complete F4/F6 ratio short of 2x. The cold n9 control shows the algorithm
can run through relation collection and final LA, without implying that
this tiny instance predicts larger curves. Next work should expose F6 in
the native `ecbench` pipeline, then use isolated-host receipts and a
same-target strong-rho arm before any attack-level speed claim.
