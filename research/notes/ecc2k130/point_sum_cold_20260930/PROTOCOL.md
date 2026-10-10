# Fresh-Q cold CPU gate for the compact full-point query

Status: preregistered before generation of this panel's Q or any cold timing.
The candidate is the optional `KIC_QUERY_BACKEND=point_sum` merged in
[PR #1070](https://github.com/aburan28/crypto/pull/1070), after the
protocol-only [PR #1068](https://github.com/aburan28/crypto/pull/1068)
and independently replayed public-Q development pilot. No result from that
pilot is eligible here. Build both compact arms from the *same binary* at
merged source commit `d8d67d3b683334db105195b0fc74d8e660dc4b5e`.
Its compact example SHA-256 is
`7b2393731cd7da047b9404585c85dc2f384ee857e43778661e806649c6cf552d`;
the matched rho v3 source SHA-256 is
`98b6e8a27d821ebfc7ae716410184dffe89aadf27dfe36bcac8ce834c39cf04c`.
The pinned runner helper SHA-256 is
`bb57f12a4f57b6921674af026a5993c661751bdcbc47d0984b6f32af64b014b8`.
The outcome must also record the Cargo.lock, compiler, binary and every
input hash. Any source change, including an apparent correctness fix,
requires a new protocol and new Q.

## Hypothesis, cells and boundary

Storing full rational pair points saves query arithmetic but adds two
group additions per indexed state, y storage and conversion per candidate.
The hypothesis is that this tradeoff lowers complete cold CPU in at least
one rank-heavy n41/n53 cell without losing rank or a target scalar.
The six fixed cells reuse the previously selected useful-base sizes and
rank seed 7. Each uses the same K, public Q, W64 query window, base
constructor and prefilter policy in control and candidate. The control
policy is the best retained current-source choice per cell: S3/off for
n37 and the n41/n53 single-target cells, S3/blocked for n41/n53 batches.
Rho is the 32-walk normal-basis signed-Frobenius v3, DP bits 4, on those
same Q. No known-answer scalar enters either solver.

| Cell | K | Q per process | S3 prefilter | Cold blocks |
|:--|--:|--:|:--|--:|
| n37/L1 | 7 | 1 | off | 20 |
| n37/L1024 | 42 | 1,024 | off | 5 |
| n41/L1 | 85 | 1 | off | 5 |
| n41/L1024 | 255 | 1,024 | blocked | 5 |
| n53/L1 | 220 | 1 | off | 5 |
| n53/L1024 | 440 | 1,024 | blocked | 5 |

For each n/L, generate Q and a separate verifier-only label file once
with the pinned rho v3 `KIC_RHO_GENERATE_ONLY=1`, corpus
`point-sum-cold-n{n}-L{L}-20260930-v1`, and seed `2026093007`.
The command is `rho n 0 signed_frobenius L 2026093007`. Independently
check every `[d]G=Q`, subgroup order and field modulus; ensure each Q is
disjoint from all available prior published point-only corpora at its n.
Commit the exact point/label files and hashes in `FROZEN.json` before the
first timed arm. A collision or generator mismatch invalidates the
corpus; record it and freeze a new protocol rather than selecting around
the result. Labels remain inaccessible to the producer and rho commands.

## Cold execution and independent replay

Each block runs fresh child processes, rotating
`control_a, point_sum, rho, control_b` by `block_index mod 4`.
No process or prebuilt index is reused. A child charges startup, curve and
normal-basis setup, point-defined base, root index plus candidate y points,
full-rank relation collection, linear algebra, every target descent,
group verification and output. Use one pinned Linux x86-64 core, one
worker, a reserved-core isolation monitor, `wait4` child user+system CPU,
wall, ru_maxrss and `/proc` peak RSS. Cap each arm at 900 seconds and
5 GiB address space/observed RSS, with a 90-minute cell-job cap. Stop a
cell on its first failed child, retaining every completed child and the
failure. Never drop a block or rerun a selected cell in this protocol.

The independent verifier must read only the frozen point file, separate
labels and raw traces. It recomputes every compact base point, selected
four-point relation, augmented rank row, representative log, target log
and `[d]G=Q`; it also checks every rho scalar and Q. Compare base hashes,
state/root table counts, rank attempts, useful K and target hit/miss
counts between compact arms. Root-order witness or probe differences are
allowed only when both four-point witnesses independently replay and no
previously solved Q becomes a miss. Preserve raw stdout, stderr, base,
rank, targets, host, isolation, commands, environment, resource record,
source/input/binary hashes and second-host replay receipts in the outcome
PR. The verifier must reject a changed raw byte or mismatched Q.

## Decision and stop conditions

For each block use `sqrt(CPU_control_a * CPU_control_b)` as the control.
Report paired candidate/control and candidate/rho ratios, median and
two-sided 95% log-t intervals, in one six-row **complete-process CPU**
table. The control-B/control-A median must be in [0.9,1.1], its paired
interval must include 1, and the reserved-core monitor must report zero
contention. Otherwise the whole cell's timing is ineligible. Repeat Q
within blocks quantifies execution noise, not target-distribution
uncertainty. The full-rank attempt floor K+L, extra indexed-point memory,
field-product counts and stage timers are diagnostics only.

A candidate/control interval wholly below 1 is a cell-specific
engineering gain. A candidate/rho interval wholly below 1 is a tentative
native CPU crossover and triggers a separately frozen disjoint-Q and
second-host confirmation, plus the still-missing calibrated all-phase
operation unit `S=operations/sqrt(r)` and n83/n131 transfer analysis.
Any missing scalar, replay failure, timeout, OOM, source mismatch or
eligibility failure censors that cell; do not impute a speedup. If all
eligible candidate/rho intervals stay above 1, report the quantitative
no-go for **this query representation at these six cells**, including
the per-cell required CPU reduction `1-1/R`. This alone does not rule out
different factor bases, higher-arity PDP or degree-263 descendant-native
policies. Update the canonical scoreboard only with fully eligible
end-to-end evidence; leave method-level speedup and n131 claim unset
until their separate gates pass.
