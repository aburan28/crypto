# Optimized continuation, frozen before execution

The existing dev-build survey and its 99-vertex checkpoint are preserved.
It has audited every one of the 51 conductor strata but not every vertex.
This continuation uses a separately named campaign and a release-profile
binary; no binary or sealed session in the earlier campaign is replaced.

Keep all survey mathematics and budgets unchanged: actual largest prime
subgroups, two public one-target workloads, seed 20261017, one alternating
round, zero warmups, projected bases 16/64 with MITM/subtraction, raw-base
MITM control, matched rho/BSGS and rho A/A, 4,096 relation trials and a
3-second child watchdog. Freeze the new binary, exact inputs, source and
driver identity into the new campaign. Compiler profile is Cargo release
(opt-level 3, no added target-cpu flags), debug symbols disabled, incremental
off, two build jobs. Preserve build failures and resource observations.

Hypothesis: completed deterministic operation reports are independent of
compiler optimization. Before launch, pass the small control and the same
three root-vertex pilot with the release binary. Audit all sessions and replay
both successful and exhausted deterministic reports. Match records across
profiles on workload id, arm, round and algorithm seed. For pairs with final
deterministic reports, require identical outcomes, recovered scalars, total
GAE, phase group/native counts, field operations, counters and factor-base
identities. Any disagreement halts launch and is retained. A censored record
has no final report to compare; expose that missing pair, never call it equal.

Do not infer an optimized runtime speedup from these debug/release observations.
The primary result is counted work and correctness. L2, qualified A/A and
adequate rounds remain necessary for runtime comparisons; absent admission,
retain elapsed values as descriptive only. A changed watchdog outcome can
change coverage even when completed operation counts agree. Keep profiles
separate in all summaries, reports, means and sample denominators.

After the matched controls pass, run the whole frozen 24,823-vertex survey
queue with all 51 stratum representatives first. Use a 24-hour invocation
limit and 2 GiB spool ceiling with the existing 64 MiB headroom guard. Pause
cleanly on a guard, correctness/audit failure, STOP marker or signal; retain
all attempts and the completion ledger. No unattended Drive upload or external
worker is implied. Provisioning and original full-matrix gates remain open.

Success is honest per-vertex coverage with checked successes or explicit
failures and reproducible records. No faster-than-rho, scaling, finite-field
IC, complete original A--E, full eight-target/five-round/million-trial matrix,
independent implementation replay or exact representation-registration claim
is made by this optimized bounded survey.
