# Fresh natural CryptoMiniSat preparation at a one-million-conflict cap

This is a new target-free, source-bound registration. It follows the
[disclosed two-instance pilot](../native-sat-budget-pilot-v1/RESULT.md), which
showed that this cap can resolve one previously inconclusive negative and one
previously inconclusive feasible source instance. Those selected inputs do
**not** estimate the natural relation rate. The earlier 512-query CMS and F5
registrations, the disclosed target controls and all three historical
confirmation sets remain consumed and closed.

## Frozen question

Can external one-thread CryptoMiniSat, with the exact accepted Weil-descent
exporter and a 1,000,000-conflict limit, produce 29 independent verified
relation rows and independently replayable factor-base logs on a **new**
natural n17 panel? The input law is `probe_scalar(2026100517, trial, 65587)`
for trials 0 through 511, ordered as generated. The public query point is
derived by the native worker; no target scalar is used as a solver answer.
The exact [config](config.json) fixes 512 attempts, export nonce 2026100521,
one CMS thread, a 30-second exporter watchdog, a 120-second solver watchdog
per attempt and a 12-hour outer watchdog. The same n17 curve, geometric
three-summand equation, six-dimensional factor-base construction, subgroup
filter and sign/Frobenius folding are checked from the frozen source. The
expected inventory is 63 geometric points, 62 usable subgroup points and 29
effective columns. These are three different quantities.

The worker must finish all 512 starts/completions even if rank reaches 29
early. It stops at the fixed count or the outer watchdog, with no refill,
resume or retry. Preserve every witness, source UNSAT, conflict-budget
inconclusive, watchdog timeout, model rejection, OOM and failure; do not
reinterpret an inconclusive result as a negative. A partial prefix is a
retained partial result. The independent geometric pair-complement oracle
checks claims and labels feasible misses but never contributes a relation
row. No planted input, F5 output or old CMS output fills the matrix.

## Admission and interpretation

The original frozen audit must pass source and binary inventory, exact
external exporter/CMS roles and child drains, all retained inputs and
stdout/stderr, point law, model-to-full-point lift, relation verification,
rank, all column logs and exclusive phase closure. Promotion to a later
complete SAT one-target experiment requires the full fixed panel and rank
29/29 with every recovered log independently replayed. Otherwise the result
is a trustworthy failed or incomplete preparation, not a SAT one-target
solver. Report witness, accepted-row and novel-rank rates per all 512 queries
with uncertainty, plus outcomes and costs for failed attempts and memory
where available. A zero-useful-row cost stays undefined, never zero.

This fresh seed is different from the 2026100311 seed in the old paired F5/CMS
panel. The old and new preparation costs are not same-point paired costs.
Any changed yield may reflect input variation as well as budget; the
disclosed two-point pilot tests only those exact source instances. Host CPU
wall-time ratios are exploratory because [this host](host-context.json) has
no auditable exclusive-CPU isolation receipt. No rho or IC speedup follows
from target-free preparation.

## One-use custody sequence

1. Commit and push this protocol, config and all worker/controller source.
   Freeze from that clean source commit with `ordinary-control-freeze`,
   `validation_only=false`, checked Homebrew Rust tools and the accepted
   external native archive. Freeze may build and run `build-identity`; it
   must make zero scientific worker or SAT calls.
2. Independently inventory the unconsumed original capsule, publish the
   **full** scientific archive with `ordinary-control-publish-registration`,
   and data-only replay it separately. Commit and push the external seal,
   publication and replay before dispatch.
3. From the original frozen controller, run `ordinary-control-execute` once
   with that exact publication and external seal, under the Mac's shared
   `busy` serializer. The command consumes a create-only claim before worker
   dispatch. Audit with the same original frozen checker, then archive the
   original execution and audit bytes without editing them.

Only after a passing full-rank SAT preparation may a **new** one-target
registration bind these logs to a public point and compare the complete SAT
pipeline with F5, the incumbent and strong same-point rho. That later
comparison needs a point-exposure census, canonical EC1/IC1/workload/run
identities, exclusive online phase clocks, independent scalar replay and
host-level isolation before any CPU speedup claim. A local stage winner does
not automatically become the whole-pipeline winner.
