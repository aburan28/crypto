# F6-IC: geometric closure inside a Boolean Macaulay search

Status: registered design and bounded experiment **before any F6-IC execution**.
"F6-IC" names the proposed solver in this repository, not an established
successor to Faugere's F5 and not a claim of a new complexity bound.

## Problem and boundary

The IC oracle needs one verified relation `Q = P_1 + ... + P_m` with each
`P_i` in the exact usable factor base. Its present Gröbner solver reduces the
Semaev Boolean system until it finds a complete assignment, then checks the
factor-base points and group identity. A complete Boolean root that has no
factor-base lift is useless for IC. The hypothesis is that exact geometric
checks inside the splitting tree remove algebraic work without losing any
usable relation.

The comparison boundary is the same-point, one-worker, signed-Frobenius rho
online solve, with first target-dependent work through independent scalar
replay charged. The primary IC metric is the same one-target online interval,
the sum of target query, PDP including every failed attempt, relation check,
descent, and recovery check. Reusable base/log/index construction is excluded
and reported separately. Candidate, workload, and run IDs follow the
cryptanalysis repository's `AGENTS.md` convention; `fb` is actual usable
points before folding. A missing phase, censored solve, or unverified scalar
leaves the IC/rho speedup unknown.

The 2026-10-03 CPU isolation gate applies. This Mac can provide correctness,
operation counts, and exploratory wall-time rows, but cannot promote a CPU
wall-time ratio. A controlled 2x online claim requires the isolated benchmark
service or equivalent complete host-level isolation receipt. Historical n17
F5/F4 host ratios are exploratory under that rule.

## Algorithm, version 0

Run inherited F4 on the target-specialised Semaev system. At every node,
before another matrix reduction, inspect the current partial assignment.
Only variables with known Boolean values count as fixed; a variable merely
held as a placeholder for an affine definition does not.

1. **Exact support gate.** For each summand whose entire x-coordinate code is
   fixed, map its bits through the frozen subspace basis. If no exact usable
   factor-base point has that x-coordinate, close the branch. Use the
   enumerated point set, including subgroup/cofactor policy, rather than the
   nominal subspace.
2. **Residual closure.** If at least `m-1` summands have fixed x-coordinate
   codes, enumerate the exact factor-base point choices for those codes.
   For each choice compute `R = Q - sum(P_i)` and look up `R` in the exact
   factor base. Check any already fixed bits of the remaining summand's
   x-coordinate. Return a relation only after checking the full group sum.
   If all choices fail, close this branch; no full Boolean root is needed.
3. **Algebraic continuation.** Otherwise run the unchanged inherited F4
   reduction, propagation and splitting. A spent node budget remains
   `unknown`, never a proof of no decomposition.

This is a witness-directed solver: a geometric refutation means no *usable
IC relation* extends the partial assignment, even if an x-only Semaev root
does. The support gate is sound because every usable summand belongs to the
enumerated base. Residual closure is complete for a branch because every
remaining relation must have one of the enumerated fixed-coordinate point
choices and the unique residual point. The returned relation is independently
group-verified. F4's original-equation verification remains in place for
roots reached by algebraic continuation. If variable renaming or affine
elimination makes the fixed-coordinate test uncertain, skip the geometric
test and continue with inherited F4.

No group-result cache or previous target answer may be read by this solver.
Its group lookups and additions are target-dependent PDP work. The exact
factor-base x-to-points index is reusable preparation; if built lazily during
the online interval its build cost stays charged there. Preserve a portable
CPU path and the original inherited F4 engine as an explicit reference.

## Frozen pilot and decision rules

The first pilot uses the existing n17a1 prepared, certified-log state from
`research/ic_candidate_tournament_20260915/goal_20260924/prepared-f5-runtime-v1/`
with 63 geometric points, 62 usable projected points and 29 folded columns,
`m=3`, degree 3, node budget 8192, one Rayon worker, one public target,
`batch_trials=1`, and 32 target-query attempts. Its **new** public fixture
seed is `20261003019` and its algorithm seed is `20261003034`; no known scalar
is supplied to either IC or rho. The source, binary, exact point, input JSON,
base digest, candidate records, and workload record must be committed before
the first timing run. Use exactly three fresh processes per arm in the fixed
order: inherited F4 / F6-IC / rho; F6-IC / inherited F4 / rho; inherited F4 /
F6-IC / rho. Apply the same 600-second cap to every process. Retain every
failure, timeout, OOM and zero-yield attempt. No target reselection or run
replacement follows an unfavorable result.

Before that pilot, run exact small-system controls with all solutions
enumerated: compare accepted factor-base relations from inherited F4 and
F6-IC; include non-lifting x-roots, missing x-support, affine definitions,
permuted chain variables, and both values of each split. A novel F6-IC
candidate ID is issued only after actual worker wiring, source hash, frozen
base and complete stage choices are known. `PDP3f6` is the proposed compact
stage code; this protocol alone remains a design with `candidate_id: null`.

Record reductions, splits, F4 word operations, geometric gate calls and
refutations, residual point choices/additions/lookups, Macaulay build and
reduction time, every target-query PDP status, verified relation yield, novel
rank, final LA, target descent, and replay. The pilot passes its correctness
gate only if every reported relation is independently verified, both IC arms
recover the same scalar on the same point, all five exclusive phases close,
and every pilot process exits zero. A favorable operation count or exploratory
wall time is a stage diagnostic. Promotion requires a fresh target panel and
the isolated-CPU timing gate; an unverified, incomplete, or noisy row never
counts toward a speedup.

Stop and preserve a negative receipt if the exact support index does not
match the factor-base geometry, if a branch can be refuted without a valid
group argument, or if the bounded pilot exceeds its caps. Do not change the
seed, target, candidate or accounting after seeing results.
