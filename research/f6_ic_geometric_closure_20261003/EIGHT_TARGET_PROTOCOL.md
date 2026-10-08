# F6-IC variability panel: preregistration

The prior matched comparison used one public target and found a descriptive
1.54x F6-IC/F4 complete online advantage. Repetition on that point did not
measure target variation. This panel freezes eight further public targets,
with fixture seeds `20261004031` through `20261004038` inclusive, before
constructing any point or seeing any outcome. Use algorithm seed
`20261004039` for every arm. The question is the empirical range of paired
F4/F6 and F5/F6 complete online costs on those eight targets, including
non-completions. The target is a fresh hash-to-curve point; no scalar is
constructed or supplied to a solver.

Reuse the exact source, binary and candidate manifests in `three_way/`:
one certified n17 prepared log state, 62 actual usable base points, 29
folded columns, three summands, degree 3, 8192 nodes, one worker, one query
per trial, at most 32 trials. Each candidate uses its production search and
split rule. Freeze one workload ID and three byte-hashed inputs per point.
Commit all fixtures, workloads, inputs, source/binary hashes and this
schedule before the first timed process.

For odd target indices, execute F4, F6, F5, F6, F4; for even indices, F6,
F4, F5, F4, F6. Thus F4 and F6 each have two fresh processes per target,
and F5 has one. Use `RAYON_NUM_THREADS=1`, disable artifact cache, and cap
each process at 180 seconds. Keep every stdout, stderr, exit status and UTC
start/end. Every verified online interval includes all failed target PDP
attempts and must equal the sum of the five exclusive phases. Keep timeouts,
exhaustions, OOMs and unverified scalars as outcomes, never as wins.

For each target with all three solvers verified, use the median of its two
F4/F6 process timings and its one F5 process timing. Report all eight target
rows; across complete targets, report the minimum, median and maximum paired
F4/F6 and F5/F6 ratios and the count on which F6 is faster. The two F4
replicates provide an A/A timing-noise diagnostic; report their spread.
Also report reductions, splits, F6 additions/lookups/batch groups and PDP
share. If any F5 arm times out, report its observed wall lower bound and
leave its verified speed ratio unknown. The engineering gate remains 2x
F4/F6 complete online on at least six of eight targets; a failure does not
authorize cherry-picking a target or changing the panel.

This Mac is unisolated, so wall ratios are exploratory. These are prepared
one-target IC variants, not a full cold IC/rho comparison. The full small-
curve test is a separate protocol and must charge factor-base construction,
ordinary relation attempts, matrix and final linear algebra, target descent
and replay, with a same-target rho reference through `ecbench` before any
IC/rho claim. The best n53 `PDP4root` path remains outside this panel.
