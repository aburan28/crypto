# Frozen protocol: Boolean row representation audit

2026-09-29 UTC. Follow-up to crypto PR #912. Scope: standalone synthetic
Boolean algebra with 1 <= n <= 10 polynomial variables. No curve data,
cryptanalytic integration, scalable solver, or mathematical novelty claim.

Hypothesis: packed integer GF(2) rows reduce representation overhead versus
the recovered sparse arrays at identical algebraic work. This does not change
the full-ideal dimension ceiling of 2^n. Separately inspect input interaction
width and actual fill-in; neither proves a faster solver.

Four arms: unchanged sparse exhaustive/frontier, and packed
exhaustive/frontier. Preserve graded-colex column order, input order,
degree/serial queue order, products, budgets, and certificate semantics.
No multiplication cache, row deduplication, degree truncation, or adaptive
representation selection. Packed integer XOR replaces sparse symmetric
difference; products still enumerate terms and use the same rank/unrank.
This isolates the row backend, not a state-of-the-art comparison.

Freeze inputs.json before implementation. Eight original cases (n=6,8)
are recovered exactly from the preceding results.json. Four new n=8 controls:
chain-linear, cycle-quadratic, random sparse quadratic (seed 20260929), and
dense quadratic (seed 20260930). No seeds or cases selected after timing.

Correctness gates: unchanged independent certificate verifier (raw-mask
columns), equality to unchanged sparse exhaustive output, complete truth-table
solution sets including zero coordinates, rank=2^n-number_of_solutions.
Same-schedule arms must have identical certificate bytes and deterministic
counters: submitted, products_formed, rank, stored_terms, xor_steps,
peak_row_terms, peak_pending_terms and max_scheduled_degree. Include mutation,
resource-budget and seeded random controls. Any failure stops claims.

Budgets: 100000 submissions, 100000 terms/row, 1000000 stored terms and
1000000 queued terms; 30 seconds per worker. Each worker runs single-threaded
on the first allowed logical CPU. Record visible NUMA topology and attempt
memory binding if supported; report actual failures. Pinning is not exclusive
CPU reservation. Record CPU model/features, OS, Python, visible memory,
cgroup quota, load/process summaries and affinity. DDR type may be unknown.

Before A/B: five A/A pairs of the unchanged sparse frontier arm, fresh worker
for each sample. Then seven paired sparse/packed rounds for each schedule,
alternating AB/BA. Each worker does three fresh computations, charging
construction, normalization, conversion, products, scheduling, elimination
and certificate creation. Separately time independent verification and report
compute+verify totals. Process startup/import, JSON and fingerprints excluded
from these stage times. No warmup is discarded. All three results must match.

Record minimum, median and all raw samples. Median paired ratios and a
deterministic percentile bootstrap interval (10000 resamples, seed 912) use
the seven worker means; these exploratory intervals are not corrected for
multiple comparisons. Report A/A spread and keep wall-time claims local to
this virtualized host. Mathematical row-operation counts are deterministic
diagnostics, not calibrated hardware work or ECDLP operations.

One separate traced worker per arm/case records peak Python allocations;
untraced workers record process ru_maxrss in KiB (includes interpreter/imports
and verifier; may be too coarse to resolve differences). Never infer a memory
gain from payload bytes alone. Record sparse index payload and packed minimal
bit payload separately from Python heap and RSS.

Structural diagnostics: input degree histogram, equation-level primal graph,
deterministic min-fill order/width and fill edges (a heuristic upper bound on
treewidth), maximum input-row terms versus observed intermediate row terms.
These diagnostics do not implement chordal elimination.

Success: exactness and identical work on all inputs, plus a packed/sparse
compute-time paired interval below 1 on all four unseen controls for a named
schedule and a median decrease greater than 10%. Otherwise record mixed or
negative results. Stop after this suite; do not tune based on observed times.
No ECC improvement, m=83 confidence, asymptotic or novelty claim is admissible.
End-to-end cryptographic costs and speedup remain null.

Source and input hashes are committed before measurement. Preserve the prior
artifact and every run directory; never overwrite delivered measurements.
