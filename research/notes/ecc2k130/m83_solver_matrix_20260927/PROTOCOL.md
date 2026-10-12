# Frozen m=83 solver comparison, 2026-09-27

## Hypothesis and scope

At the same Koblitz subgroup, factor base, target, and resource limit, compare
the correctness and cost of a fixed-phase Semaev-S4 decomposition oracle using
SAT, a batched Boolean Macaulay/F4-style reduction, a signature-based Boolean
F5B prototype, and exhaustive Boolean search (FES). Separately measure the
existing ternary-inline and implicit phase-choice SAT encodings. This is a
**solver-stage diagnostic**, not an index-calculus speedup or completed ECDLP.
The null is that no algebraic solver has a verified cost advantage over the
exhaustive baseline. Failed setup, unsupported engines, OOM, and timeout are
outcomes, not grounds to discard inputs.

## Frozen inputs

- Field: `GF(2^83)` in the polynomial basis modulo
  `z^83 + z^45 + z^2 + z + 1`; curve `E0: y^2+xy=x^3+1`, order
  `4 * 2417851639230796216685689`; clear the cofactor and verify every
  generator and point is in the prime-order subgroup.
- Normal-coordinate and factor-base code: copy the archived `degree-scale`
  probe and its `heldout-phase` vendor modules verbatim into this PR, record
  SHA-256 for every source. `CyclotomicFactorSpace(83, 2)` gives the same
  two-bit quotient payload in every arm; use seeds `260938` (development) and
  `260939` (holdout), with one setup per seed shared across engines. Freeze the
  chosen normal generator, subgroup generator, representatives, and Frobenius
  eigenvalue in each output manifest. No input is selected by solver runtime.
- Fixed-phase comparison: use phase triple `(0,1,2)` for the three points,
  Semaev `S4` with Boolean descent, and the six payload bits as the only
  solver variables. For each seed choose the first two distinct proper planted
  triples in lexicographic iteration over usable base-point lifts, requiring
  nonzero target and intermediate x. Then draw four independent natural
  subgroup targets using Python `random.Random(seed*1000+832)`, nonzero scalars
  from `[1,r)`, in draw order. The random seed and generated scalars are
  recorded, not supplied as witnesses to a solver. Run all six inputs for
  every available solver, including no-decomposition cases. If fewer than two
  planted triples exist, record `setup_failed` for the affected input slots.
- Phase-choice comparison: same first planted target, seeds, and phases
  `(0,1,2)`, compare `ternary_inline` and `implicit` encodings from the
  archived scale probe, on equal 3-phase workloads. The 95-variable inline
  and 347-variable implicit encodings are *not* equivalent to the six-variable
  fixed-phase workload; report them in separate tables. F5B's `2^nvars`
  layout precludes both, so record `unsupported_memory` before allocation.
- The source hashes and exact generated targets are captured in immutable
  JSONL receipts. Reproducibility does not rely on an unpublished scratch path.

## Engines, budgets, and checks

For each fixed-phase query, build one canonical Boolean ANF, hash its sorted
equations and use identical equations for the four engines. SAT translates
ANF to Z3 Boolean And/Xor, F4-style performs degree-batched Macaulay row
reduction in the Boolean quotient, F5B is the existing squarefree signature
prototype (pin/copy its source), FES exhaustively checks all 64 assignments.
The Boolean relation is necessary, but not sufficient, for a group relation:
lift every returned model in the prime-order subgroup, check that its three
points add to the target, and block algebraic false positives. Independently
enumerate the 64 assignments and compare the verified solution set, including
empty cases. Log the candidate count, false-positive count, equation hash,
operation counters when available, rank contribution, timing for preparation,
equation build, solver setup, solve, lifting and verification, plus peak RSS.
Do not count planted inputs as natural relation yield or rank.

Each subprocess has a 1 GiB virtual-address limit and 70 s wall limit.
Equation building and solver setup each have a 25 s deadline, and each solver
check has a 5 s deadline per query; an engine that cannot enforce a solve
deadline must run in a killed child. Set `PYTHONHASHSEED=0`, SAT seed 13,
one thread. Preserve timeout, OOM, unavailable dependencies, and exceptions
in the receipts. Record exact interpreter, platform, backend versions, source
hashes, input digest, and elapsed and CPU time. For phase-choice comparisons
the same limits apply to one planted target per arm and seed. No results from
the earlier `degree_53_83_phase_scale` report are silently combined with new
measurements; cite those earlier stage results separately.

## Boundaries, success, and stop conditions

The exhaustive fixed-phase reference visits at most `2^(3*2)=64` Boolean
assignments; its measured field/curve-operation counts and verified answers
are the matching local boundary. For each solver report its cost / FES cost
in *the same charged unit* where measurable. Field multiplication, curve
addition, SAT decisions, F4 rows, and F5 pairs are different units; leave
cross-backend operation ratios null without an observed conversion. Report
paired wall and CPU times only as practicality data, never as a DLP speedup.
The full ECDLP boundary is a matched Pollard-rho run including the applicable
automorphism quotient, with `S=total_operations/sqrt(r)`. No solved ECDLP,
rank-complete relation matrix, matched rho reference, or meaningful `S` is
expected at this tiny factor-base size; record each as null. If the fixed
phase holds yield zero natural verified relations, report that exactly.

Success for this bounded diagnostic means all planned inputs produce a
verified result or an explicit bounded outcome, all satisfiable claims match
the exhaustive oracle, every saved candidate passes group verification, and
all row/phase/cost accounting is reproducible. An IC improvement toward
ECC2K-130 remains **unestablished** unless a separate m83 full-pipeline
matched-target verified-DLP/rho study clears the `AGENTS.md` gate. Stop an
individual case at its deadline; do not increase it after looking at results.

## Execution plan

Commit this protocol before executing the new cases. Commit executable source,
frozen config and hashes, raw receipts including failures, independent
verification, a one-table-per-workload comparison and decision to this PR.
Update the canonical scoreboard only for measured full-pipeline performance
changes; do not add a stage-only speedup row.
