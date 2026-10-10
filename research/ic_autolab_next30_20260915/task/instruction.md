# Improve complete index-calculus solve cost

Work in `/app`. The incumbent is the independently verified fast point-lifting winner:
pair_table, sparse linear algebra, batch_trials=1, collection_window=8,
max_trials=4096, summands=3. The new goal requires another 30% reduction in
complete instruction cost (candidate/incumbent <=0.70), fresh confirmation and
replay, all correctness checks and unchanged support. Do not claim an arithmetic
complexity improvement or a rho crossover from these small fixtures.

## Iteration

1. Read `/app/AGENTS.md`, candidate.json, and the relevant Rust implementation.
2. Write a concrete bottleneck hypothesis and expected phase change to NOTES.md.
3. Edit candidate.json or allowed Rust sources. Then run `python3 /opt/harness/probe.py`.
   A probe uses 8 public development targets, 3 repetitions, incumbent and candidate
   plus matched rho. It prints the location and summary; raw evidence is retained.
4. Inspect cost and correctness. Revise or revert, and repeat within the budget.
   Read recent `/app/probes/*/result.json` before choosing your next change.
5. The harness retains the lowest verified development-cost checkpoint in `/app/best`.
   Leave a concise research note. Final holdouts run after editing has ended.

There are at most 12 development probes and 60 minutes of total agent time.
Each model session has at most 15 minutes. Keep shell commands under 12 minutes.
Do not launch nested agents or background processes. Do not install dependencies.

## Measured bottlenecks and hypotheses

Prior development profiles show general-purpose scalar multiplication dominating
cofactor-class validation, relation checking and column-log certification.
Candidate A dispatches exactly representable points to the existing single-word
FastCurve implementation, retaining the general fallback and every required check.
The matched development result is exploratory; the candidate is locked before
fresh confirmation. This launch uses Harbor's no-op agent and performs no model calls.

## Fixed experiment contract

Allowed code edits: src/binary_ecc.rs, src/cryptanalysis/koblitz_index_calculus.rs,
src/cryptanalysis/koblitz_sparse_la.rs, src/cryptanalysis/koblitz_factor_base_search.rs.
Other source files, Cargo files, build flags, worker instrumentation, oracle and
measurement code are frozen. Configuration may change batch_trials, collection_window,
solver (pair_table or enumerate), linear_algebra (dense or sparse), and sparse options.
max_trials=4096 and summands=3 remain fixed. Factor-base seed, size and support are fixed.
Do not add answer tables, special cases for known targets/seeds/sizes, external data,
profiling detection, counter control, skipped phases, caches across runs, GPU code,
or a substitute generic DLP solver. Preserve the IC relation/rank/descent pipeline.

Every score charges startup through termination, including setup, failed attempts,
all required checking, log recovery, reporting and cleanup. The unit is Valgrind
3.22 amd64 Ir (user-space guest instructions), not curve operations. Native timings
are diagnostic only. Unverified solves, changed support, missing phases, timeouts
and missing costs cannot win. Development results are exploratory. Final verification
requires fresh targets, 60 confirmation targets in five curve cells, matched repeats,
a 95% paired interval below 1, no cell over 10% worse, and a fresh-process replay.
