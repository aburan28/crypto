# Supplementary retained five-summand gate

Frozen before launch on 2026-10-10, following the successful factored-S4
construction and v2 object replay. This uses the existing owner-authorized
two-hour-total local allowance, with no new allocation or cloud worker.

Hypothesis: the clean source-matched native chained-S3 worker can construct
the all-affine m=5 K64 retained-domain model under one Docker CPU, 4 GiB,
zero swap, no network, and a 60-second outer wall cap. Fixed limits are
150,000 variables and 3,000,000 finite-domain clauses. These limits admit
the documented algebraic source bound while remaining externally capped.

Only if capacity passes, run a separate zero-trial preflight under 60 seconds.
Only if both pass and 620 seconds remain before the total deadline, run one
public target-zero trial with one model, 10,000 conflicts, a 590-second worker
wall cap and a 600-second outer cap. Use the exact K64 hash-policy seed
2026100801 point set, source commit 85e14930ff0aab0d50e8fb12ee7ca382e5e2f8e2,
the frozen original public corpus and known-answer sidecar, and the same
immutable locally present image as the v2 guard. Known answers are used only
after candidate output. No source changes or cap increases during a run.

Preserve commands, hashes, container state, stdout, stderr, events, worker and
outer receipts for every outcome. Check original group verification and
relation/rank counters before admitting a worker summary. A variable/domain,
model, conflict, memory or wall cap is UNKNOWN, never UNSAT. A checked target
would remain a stage-only result: column-log verification, fully charged
construction, independent-host replay and matched total-runtime comparison
are still required. Total runtime and selected winner remain null.

This is an additional fixed source gate, not execution of every frozen tuple
and not a change to the published v1/v2 ordinal mapping. L0 shared-host times
remain budget and construction diagnostics; no comparative timing claim.
