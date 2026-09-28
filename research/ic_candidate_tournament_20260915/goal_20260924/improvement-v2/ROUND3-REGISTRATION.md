# Third and final bounded IC attempt: premeasurement registration

The exact [eleven-arm panel](round3.json) is frozen by its byte SHA-256,
`30b04a370e77e2a2bb7c01b7e96fd1ad8b1785104240113ac71f85eac28beca8`.
The [version-two protocol](PROTOCOL.md) remains the decision rule. This file is
a registration, not an improvement result. Seed `2026092553`, 72 fresh
confirmation points (12 in each of six cells), three process repetitions per
point, 3,480 maximum native/profile pairs, 8 GiB, 180 seconds per child, one
CPU and one worker thread are fixed. No prior measured execution is repeated.

The previously qualified `pairinv` cold incumbent, `prepared_both` online IC,
and separate cold/online rho references remain unchanged. The candidate source
is the independently checked round-two Rust worker, with source-manifest SHA-256
`8582e4ab4b63e98696a0ff00ee296e2902923a2c39c325ab0f3e3950ffbb2b28`.
`round3-v1` is a new **configuration registry** for that exact source, not a
claim of new executable code. Its ten new complete-pipeline combinations were
chosen using only the round-one and round-two **development** stages. None
duplicates an earlier registered configuration. Confirmation and replay data
were inspected only for the prior decision and are excluded from tuning.

The development evidence showed competing costs: half pair tables improved
single-target online time but increased complete cold instructions; five- and
six-orbit bases with word rows saved some cold cost but lost online time. This
panel combines three- through eight-orbit Frobenius bases, heuristic/half/cover
pair-table coverage, full/word/bounded scalar rows and one two-row collection
batch. The deliberately small three-orbit and broad eight-orbit arms preserve
exploration around that tradeoff. Every combination builds its own factor base,
collects verified ordinary-query relations, solves its own final scalar-field
matrix, descends the supplied public point and replays the recovered scalar.
No candidate inherits another candidate's cheapest measured phase.

The tournament's family key now recognizes exact factor-base recipes, orbit
target, pair-table coverage, collection batch and adapter. That prevents a
different mechanism being silently counted as the same family when the
development portfolio reserves
its diversity slots. The frozen evaluator retains this selector; prior rounds
are reverified with their original frozen evaluators and decisions.

The [runner](../../run_improvement_v3.py) restores both sealed prior rounds,
extends their exact public-point exclusions (including failed preparation),
and pins the three supplemental readiness/qualification/control fixture
corpora. It independently replays both prior decisions before generating a
new point. The accepted reference archive supplies the exact frozen workers.
The [workflow](../../../../.github/workflows/ic-improvement-round3.yml) runs
archived fixed-vector candidate controls during review and starts the measured
round once when this new panel first lands on `main`. The workflow retains
partial command logs and every failure; an interrupted or incomplete schedule
is not turned into a win or automatically retried.

Only independently verified complete one-target results can enter promotion.
The headline compares online IC and rho on the same supplied point. The
supplementary cold instruction and native-time gates both require at least a
20% improvement against the qualified incumbent, confirmation and replay,
the original familywise rule and the maximum per-cell regression rule. A
negative result retains the incumbent. The fixed small-field panel, even with
a win, cannot establish the globally fastest index-calculus algorithm or an
ECC2K-130 crossover. F4, F5, SAT and large sparse LA are already represented
in the broader generic autolab but are not promoted by this specialized
pair-table panel; their own complete, source-bound qualification remains
future work.
