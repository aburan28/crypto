# Target-blind shared-rank folded-table gate (preregistered)

This protocol was committed as `5d13c68ea7d9ac89df91d44e2f418a32ccdcc7f7`
before either n37 producer run. The preserved outcome and the later
base-wall-timer correction are in [RESULT.md](RESULT.md).

This is a bounded native Rust correctness and cost-accounting gate for the
shared setup required by a future many-target n37 index-calculus comparison.
It is **not** a new IC/rho performance comparison. The prior
[full-rank n37 one-target result](../full_rank_compact_orbit_ic_20261003/RESULT.md)
reached rank 43/43 on 16 measured Q/round pairs, but recomputed the rank
system for each target. The [checked-column change](../full_rank_column_check_20261003/DECISION.md)
points to a shared, target-blind rank database as the next necessary step.
The older [cold W64 batch panel](../disjoint_cold_v2_outcome_20261001/RESULT.md)
already priced a different complete 1,024-target policy against batched rho;
it is neither this folded table nor this rank producer.

## Frozen hypothesis and input

H1: the exact 42-column `compact-orbit-scan` source base and counted
`mitm-frobenius-counted:m=3` folded oracle can collect rank 42/42 from
target-blind `[a]G` probes within 1,000,000 trials. Every accepted witness
must sum to its probe as a full point, every matrix row must match its
published witness, and all 42 solved column logs must satisfy
`[coefficient × column_log]G = [cofactor]P` for a column representative.

The curve is ICV1 `icv1-f2m37-tm534059-32aad96b`, `a=0`, `n=37`, prime
subgroup order `r=230603167`. The base is
`compact-orbit-scan:columns=42,raw_x_cap=1000000`; the oracle is
`mitm-frobenius-counted:m=3`. The rank matrix is dense incremental Gauss
with 42 unknowns, not the 43-unknown target-coupled matrix. The rank seed
is decimal `202610031137`, fed to `StdRng::seed_from_u64(seed XOR
0x534841524544524b)`; each trial draws `a` uniformly from `1..r` and
queries only `[a]G`. The cap is 1,000,000 trials. No public Q, fixture
scalar or target-dependent state is an input to this gate. The Cargo.lock
and source commit used will be recorded in the result before acceptance.

The base's ordered points, column indices and coefficients are hashed as
SHA-256 of `ic-shared-rank-base-v1\0`, then the 64-bit little-endian column
count, then for each ordered point: 64-bit little-endian x and y, one byte
for infinity, 64-bit little-endian column and coefficient. The exact
digest is a measured identity output, not a value chosen after rank is
known. The existing archived-source parity test must pass before this
new digest can be treated as the original 42-column policy.

## Evidence and stop rule

The Rust producer emits every accepted relation's trial index, scalar,
full-point coordinates, ordered witness indices, exact dense row, right
hand side and rank transition, plus all failed-probe counts and per-phase
operation/native-work ledgers. A separate Rust verifier must reconstruct
the base and every `[a]G` with the general binary-curve group law, check
the witness and row, solve the matrix independently, and check every base
log by full-point multiplication. It must not call the producer's folded
table, fast group law or elimination. Both raw output and replay receipt,
their SHA-256 hashes, source/input hashes, complete failures and the
decision will be committed in this PR. Preserve any failed attempt.

Pass only if rank is 42/42, all 42 column logs are independently verified,
every accepted relation replays, and no inconsistent or invalid row is
reported. An exhausted trial cap, any mismatch, or a failed independent
replay is a **gate failure**, not a partial success. Record all of these
outcomes. Report the complete target-blind setup cost (base, folded table,
failed/successful probes, relation verification, elimination and base-log
checks) in counted group-addition equivalents alongside unpriced native
work and wall diagnostics. Do not divide it by single-target rho or call
it a batch speedup: no target descent, actual batched rho, memory cap or
fully priced native work is measured in this gate.

If H1 passes, the next separately preregistered campaign freezes new
orbit-disjoint point-only Q, a residual policy, batch lengths, a native
strong signed-Frobenius batch-rho reference and isolated cold-process
timing. It charges this rank setup once and every descent, checks every
target log and reports one matched whole-batch cost table. If H1 fails,
stop this folded-table batch route and diagnose base support or relation
rank before spending on target blocks. Neither outcome alone establishes
an ECC2K-130 crossover or n131 transfer.
