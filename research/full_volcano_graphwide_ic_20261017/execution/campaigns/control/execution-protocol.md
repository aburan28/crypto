# Execution addendum: controls and bounded graph-wide survey

This addendum preserves the original full protocol. It does not waive its
registration, optimized-build, wall-clock qualification or launch gates.
The assumption that an existing IC recipe had already been verified on these
large-cofactor models was unsupported. The original raw-abscissa recipe is
retained as a control; the new, separately identified projected-abscissa base
must establish correctness before a comparison can be promoted.

The hypothesis tested here is that a cofactor-projected, negation-paired base
supports checked subgroup relations on every stored vertex. This is not a
faster-than-rho hypothesis. Each successful answer must satisfy [x]G=Q in the
runner and auditor. Any wrong answer or failed audit halts the driver.

Modes are separate immutable campaigns:

* `control`: p=113, a=19, b=2, N=112, r=7; two public targets, one round,
  no warmup. Projected size 2, scan cap 64, 1,000 relation trials.
* `pilot`: the first vertex of each complete graph, two public targets,
  one round, no warmup. Projected sizes 16 and 64, scan cap 20,000,
  4,096 relation trials. Each child has a 3-second watchdog.
* `survey`: all 24,823 vertices in round-robin graph order, with precisely
  the pilot's parameters. It measures bounded attempts, including exhaustion
  and watchdog censoring; it does not replace the full million-trial matrix.
* `full`: the original grid 16/32/64/128, 8 public targets, 5 rounds and
  1 warmup, 1,000,000 trials. Planning is supported. Launch is deliberately
  refused until the original documented gates have been satisfied.

Arms: matched rho.negation, BSGS, an A/A rho repeat, the unchanged raw base
with MITM, and projected bases with MITM and subtraction. All direct arms
use checked full-rank incremental linear algebra; each base column's log is
point-checked. No fallback solver is relabeled index calculus. Cold setup is
rebuilt and charged for every one-target run. Generic sqrt(r) costs and the
automorphism-aware collision boundary are recorded by ecbench.

Target seed is 20261017, algorithm-order seed 1, alternate ordering, all
method parameters are resolved into the sealed ecbench plan. L2 is requested
with whole-core auto placement and busy-host downgrade recorded. The present
unoptimized binary is suitable for deterministic counts and correctness,
not optimized runtime claims. Two targets/one round do not qualify wall
comparisons. Retain raw wall observations without calling them admitted.

The native driver freezes an immutable queue, binary hashes and exact input
file hashes. It uses SQLite checkpoints and an advisory single-writer lock.
Each attempt receives a new session directory; no sealed file is overwritten.
Run and replay are sequential. Every deterministic measured record is
replayed. Failed attempts remain visible and are retried into new directories
on resume. Finite-field IC is explicitly unsupported for these actual
embedding degrees; the exact per-vertex applicability records are preserved.
Exact representation registration remains pending where indicated by inputs;
these are diagnostic records, not qualified scientific promotions.

Disk and elapsed-time limits are checked between vertices. A vertex may
overshoot the elapsed-time limit by one bounded session and replay. The disk
guard reserves headroom for one attempt; it is not a filesystem quota. A
stop file or signal requests a graceful pause. No background uploader is
assumed: the process pauses with data intact for explicit Drive archival.
Timeouts have unknown relation counts when the child emitted no final report;
they are not zero-yield observations. Exhausted reports retain counts.

Report verified/exhausted/error/timeout and not-run denominators separately.
Min/max/arithmetic mean use only present observations with their sample
counts; expose missing/censored values. Emit per-graph, stratum and arm views,
phase costs and all raw unpriced counters, with common 0/0.1/1 GAE bands.
There are no speedup or scaling claims from this diagnostic.
