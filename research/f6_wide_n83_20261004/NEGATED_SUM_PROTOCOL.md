# K_0 n=83 stored negative pair sums

Registered before fresh baseline runs on 2026-10-04. The current wide pair
index stores positive pair sums and creates a negated vector in every query
chunk. The candidate stores each negative pair sum once during construction,
retains the positive sum only as the exact lookup key, and stores pair
indexes in a parallel vector. The target residual is still `T - (P_i+P_j)`;
the same representative pair and exact replay rule apply.

Run the release index probe on the pinned K_0 dimension-8 and dimension-10
bases and public T001, with one warmup and three build/query observations per
size. Run the full dimension-12 probe once for all 8,219,485 pairs, one
exact query, and peak process RSS. Record baseline runs before editing,
then run the identical probes with the candidate on the same host and
compiler. Preserve every duration, outcome, pair count, source hash, exit
status, and failure. The one-query full-base comparison is exploratory and
may be sensitive to host variation.

Retain the candidate only if exact outcomes and pair counts match; its
full-base query time is at most 0.90 times baseline; each small-base query
median, full-base build time, and peak RSS are at most 1.10 times baseline.
Otherwise revert the implementation and retain the patch and raw receipts.
These are unisolated stage diagnostics, not controlled CPU speedups or
complete F6/F4/F5 or one-target IC comparisons.
