# n37 shared-rank point-only target recovery gate (preregistered)

Status: protocol only. No target input has been generated or queried under
this protocol. This follows the verified target-blind 42-column setup in
[PR #1302](https://github.com/aburan28/crypto/pull/1302). It is a bounded
correctness and cost-accounting gate before a separately frozen, isolated
IC/strong-batch-rho timing comparison. It makes no speedup claim.

## Fixed method and hypothesis

Use curve `icv1-f2m37-tm534059-32aad96b`, `r=230603167`, the exact
`compact-orbit-scan:columns=42,raw_x_cap=1000000` base and counted
`mitm-frobenius-counted:m=3` oracle from #1302. Build the 42-column rank
system from seed `202610031137`, cap 1,000,000 target-blind trials, and
retain its prepared folded table. Require its ordered base digest
`8460ac4c28515db701c3897a03b4ce0f28abf7cd56fd759ad095f98436a76dcf`
and all 42 independently checked column logs. The new producer receives
only a point file; it must never read the scalar fixture.

H1: all 16 newly frozen Q are recovered and independently checked within
64 three-summand descent attempts each. Attempt zero queries Q directly.
If it misses, attempts 1 through 63 query `Q + [a]G`, where each `a` is
drawn from `1..r` using `StdRng::seed_from_u64(202610031649 XOR target_index)`.
For a witness `Q+[a]G = sum(P_i)`, recover
`d = h^{-1} sum(coef_i * column_log_i) - a (mod r)` and check `[d]G=Q`
as a full point. Check the witness sum as a full point before accepting it.
Preserve misses and all target-dependent group/native work in the report.

## Input freeze and independent controls

Before querying, a native Rust preparer snapshots every tracked
`research/notes/ecc2k130/**/*.points.jsonl` file present in the
preregistration parent commit, including its path, SHA-256 and row count.
Only files for n37 contribute orbit exclusions. Add the 3,108 base points
and all 55 accepted rank-probe points to the exclusion set. The signed
Frobenius orbit key is the minimum of `x,x^2,...,x^(2^36)` in the pinned
field (sign leaves x unchanged). Candidate `j=0,1,...` is SHA-256 of
`n37-shared-rank-target-v1|j`, interpreted big-endian and reduced modulo
`r` to obtain `d`; reject zero, excluded or repeated orbits. Accept the
first 16 candidates in order. Write a point-only Q file and a separate
scalar fixture containing every accepted `j,d,Q`, rejection counts, exact
inventory and its hashes. The fixture is for replay only. Freeze and commit
both files and an independent general-curve-law input receipt before any
target producer run. Never replace a target after its outcome is known.

The producer emits the rank setup, each Q and every attempted residual
scalar, hit/miss, witness indices, recovered scalar, verification result,
exclusive phase ledgers, and whole target-dependent interval. A separate
native Rust replay uses the frozen point-defined base and general binary
curve law, not the fast producer oracle, to verify each witness, folded
coefficient, residual equation, recovered scalar, and `[d]G=Q`; it also
compares recovered logs with the scalar fixture. Mutating one witness or
one recovered scalar must fail replay. Both successful and failed target
rows remain in the committed raw output.

Pass only if all 16 inputs are orbit-disjoint under the frozen inventory,
the rank setup agrees with #1302 on all non-timing fields, every Q has a
verified scalar within 64 attempts, all counts and phase sums replay,
and the two negative controls are rejected. A miss, timeout, invalid row,
or replay mismatch fails this gate. The per-target online interval includes
all attempts and scalar replay; target-independent rank/base/table setup
is separately charged once. Report group-addition equivalents and native
unpriced work, with wall time only as an unisolated diagnostic. No rho
denominator, timing crossover or n131 extrapolation is admissible here.

If H1 passes, freeze a separate 1,024-point block with the same policy and
run the complete cold shared-rank batch against the same-Q 32-lane strong
signed-Frobenius batch rho on an isolated host, including independent replay,
memory, native conversion costs, and A/A drift controls. If H1 fails,
diagnose the exact failed Q and avoid spending on batch timing.

## Post-gate process clarification

The H1 inputs, stop rule and interpretation above were frozen before the
run and are unchanged. For the next *comparative* measurement, the primary
one-target contract takes precedence: first freeze one new Q and pair this
candidate's verified online interval with a same-point strong rho solve.
The 1,024-point experiment above follows only after that primary result;
its distinct question is shared-setup amortization. This clarification
does not reclassify the 16-point correctness cohort as a batch speed test.
