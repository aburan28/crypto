# One-target static-SAT IC pipeline on n17a1

The [fresh ordinary-query sample](../static-cms-s4-natural/RESULTS.md)
found six independently verified relations in 32 queries with the
source-receipted static wide-S4/CryptoMiniSat stage. This new
[panel](panel.json) freezes a **complete** three-summand IC attempt before
execution. It uses a distinct relation-query seed, one previously unseen
public target point constructed by a declared hash-to-curve law without a
known scalar, a separate target-descent seed, and fixed attempt/watchdog
limits. The target point is input; its construction is outside the online
clock. No sealed confirmation set is reused or opened.

The public target is derived by SHA-256 of the ASCII domain
`ic-static-sat-target-v1` followed by a zero byte, the little-endian
eight-byte seed `2026092938`, and a little-endian eight-byte counter starting
at zero. Take the first eight digest bytes little-endian modulo `2^17` as
`x`, select rational lift number `digest[8] mod lift_count` in the source
lift order, then multiply by the curve cofactor. Skip nonexistent lifts and
the identity. The first valid counter is zero, yielding `[114119,85674]`.
The generator scalar is never used in the IC attempt. The 256 planned
ordinary scalars are distinct and disjoint from the prior F5 and natural
SAT ordinary samples; this is checked again at admission.

First verify the exact n17a1 field, subgroup, generator, Frobenius eigenvalue,
63 geometric factor-base points, 62 nonidentity subgroup-usable points and
29 sign/Frobenius-folded columns against the source-bound parent evidence.
Bind the exact Rust exporter source/dependencies and the Phase-B static
CryptoMiniSat source/dependency build receipt. Preflight the copied executable
from its final path. Every ordinary relation query uses the independently
replayed trial-keyed nonzero scalar law at seed `2026092936`. Export its
public `[a]G` point as the wide-S4 XOR-DIMACS system, run one static solver
child, verify any complete model against **all** source clauses/XOR rows,
lift the three coordinates into the exact base, and group-readd the result.
Retain all 256-or-fewer attempts, including failed exports, solver errors,
UNSAT/Unknown, timeouts, nonlifting models, duplicates and dependent rows.
Do not convert a capped or failed attempt to a negative result.

For each witnessed relation, independently construct the row of
`[h]P` coefficients under the signed Frobenius orbit quotient and right-hand
side `h·a (mod r)`. Use exact dense modular Gaussian elimination over the
verified prime subgroup order, recording every matrix row, rank transition,
dependency and duplicate. Stop collection only at full rank **and** after
every solved orbit-column logarithm independently scalar-replays against
its point. The row implementation is already controlled on the 29 actual
source-bound F5 relations and reproduces all 29 certified logarithms; a
forged group row is rejected. A timeout or insufficient rank remains a
retained incomplete run.

Once the reusable base/index/log preparation is ready, start the one-target
online interval. Target query `i` draws `(a,b)` from the separately frozen
seed `2026092937`, constructs `R=[a]G+[b]Q`, and gives only `R` to the SAT
exporter. Charge **all** target-dependent attempts, including failures, to
this one target. For a checked decomposition, derive
`d = (Σ orbit_coeff·column_log - h·a)/(h·b) mod r`, then independently verify
`[d]G=Q`. Stop at that verified scalar or 64 attempts. The online interval
ends after scalar replay. Record exclusive target-query, target-PDP,
target-relation-check, target-descent and recovery-check wall phases whose
sum equals the measured online interval. Keep setup/build/precomputation,
process launch of the controller and fixture construction outside that
online interval; preserve separate process and cold-preparation costs.

The complete runner and every executed helper must be committed and hashed
before launch. The exact immutable `IC1` candidate/curve/workload manifests,
method, source manifest and one run ID are sealed with the panel before the
first job. The source manifest includes the local transitive runtime Python
imports plus the Rust source/dependency and static solver build receipts.
The run retains the candidate/workload IDs even if rank or target descent
fails. A complete result needs an independent archive replay of source
systems, group relations, matrix/rank/logs, target-query law, phase ledger and
final scalar. Failed, timed-out, OOM or unverified runs remain rows without
imputed online cost or speedup.

The subsequent comparison must run source-bound F5, incumbent IC and
single-target rho on **the same Q** in a new paired protocol and one frozen
resource envelope. The headline is verified one-target online wall for each,
with rho/IC only when both solve. Stage times or this one unpaired IC run
cannot establish a speedup or global optimum.
