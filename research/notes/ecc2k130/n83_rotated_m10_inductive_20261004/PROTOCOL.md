# Frozen follow-up: inductive-validity n83 rotated-m10 gate

## Motivation and isolation from the parent outcome

The exact-public-target complete-per-edge gate in PR #1343 ended
`BUDGET_OR_SOLVER_STOP` after its frozen 120-second CaDiCaL window.  It built
998,338 Boolean nodes and 3,702,311 clauses, then reached 255,395 conflicts
without SAT or UNSAT.  That outcome remains fixed and is not retuned.

This follow-up tests one algebraic reduction: check the ten factor points on
the curve once, then enforce the complete group law on nine additions without
rechecking the curve equation for each intermediate point.  For valid inputs,
the complete addition relation defines the true group sum, which is again a
valid point; induction therefore makes the intermediate validity equations
redundant.

This is a paired circuit/solver follow-up on the same public Q.  It may show a
representation reduction or a verified PDP relation.  It cannot confirm a
general yield rate or a faster-than-rho ECDLP.

## Proof and regression gates

Before any capped n83 solver run:

1. Exhaustively enumerate every valid point pair, candidate result and slope
   on the type-II-ONB n=3 E0 curve.  With validity checks disabled inside the
   addition relation, an assignment may be accepted if and only if the
   candidate is the numerical group sum, and every valid pair must have at
   least one accepted slope.
2. Build both n83 planted variants.  The inductive circuit's known full model
   must replay, and its node and AND counts must be strictly below the
   complete-per-edge circuit.
3. Rebuild the old complete public circuit and require its DAG prefix BLAKE3
   `689697d223df16ad1cb65edfa6f44f4ecbb3f0032ca55f8b31587da8696525cf`
   and ordinary-CNF BLAKE3
   `bd4c5148fa97d7af749c5736bd389b77fe938a2b978aed901d893b8bd463b892`.
   Drift is a hard stop.

The n=3 exhaustive test and the n83 known-model unit test were executed while
developing the implementation.  No inductive n83 CNF was exported, no exact
inductive size was inspected or recorded, and no inductive solver was run
before this protocol.

## Frozen instance, exporter and caps

Use the same curve, type-II normal basis, slots, generator, exact public Q,
ordinary Tseitin exporter, committed CaDiCaL 3.0.1 binary and one solver
thread as the parent protocol.  Keep its 8,000,000-node, 2,000,000,000-byte,
900-second cold-process, 12-GiB and 120-second solver caps.

Run in order:

1. Emit and replay the deterministic planted inductive certificate.
2. Cold-build the inductive public-Q circuit and run CaDiCaL once for 120 real
   seconds with default configuration.
3. Preserve all counts, digests, stdout/stderr, terminal status, conflicts,
   propagations, decisions, solver time and peak RSS.  If SAT, rebuild and
   independently replay the full model, all factors and the group sum.

Do not change validity policy, branch formulas, slot order, target, solver,
preprocessing or timeout after observing the result.

## Decision rule

- A representation win requires strictly fewer nodes, AND gates, clauses and
  bytes than the frozen complete public circuit.
- `SAT_VERIFIED_RELATION` requires full Boolean and numerical replay.  It
  admits a separately frozen disjoint-target yield panel; it is not a DLP or
  rho crossover.
- UNKNOWN/time limit, cap, crash or replay failure is a bounded negative for
  this inductive/Tseitin/CaDiCaL route.
- UNSAT rejects the exact tuple domain for Q, not high-arity descent in
  general.

No algorithmic speed claim is permitted without full relation collection,
rank, factor logarithms, held-out scalar recovery and fully charged cold cost
below matched strong rho.

