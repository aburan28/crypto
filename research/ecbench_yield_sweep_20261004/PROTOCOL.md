# Yield sweep: what every base, oracle and solver yields on the same targets

Registered 2026-10-04, before either session ran. The two specs beside this
file are hashed by `ecbench` at run time.

## Question

For one target on one curve, what fraction of the points an index-calculus
pipeline tries does each factor base decompose, how much does each
decomposition oracle and solver pay per decomposition, and what shape is
the algebraic system it solves? This is the yield ledger AGENTS.md §7c's
browser reads and `docs/ecbench/schema.sql`'s `ic_yield` view returns:
one row per (curve, target, factor base, oracle, solver), with the
record's new `solver` block (Macaulay and SAT statistics) where an
algebraic oracle ran. It is a **stage diagnostic**: it ranks bases and
oracles on identical targets to inform a candidate combination, and no
row in it is a speed. The whole-pipeline `S` is still carried, with the
matched rho reference, so the table stays in the repository's unit.

## Boundary

- **Floor:** the generic floor `√(π/2A)` in `S`, recorded per run.
- **Reference:** `rho.signed_frobenius_strong` on the Koblitz curves,
  `rho.negation` on the prime curves.
- **Unit:** `S`, cold, every phase charged. An algebraic solver's work is
  counted in its own unit (conflicts, monomial operations, word XORs) and
  listed as unpriced, so every descent arm's `S` is a lower bound; the
  record says so.

## Workloads

- Koblitz: `a=0 n=13`, `a=1 n=17`, `a=0 n=19`, `a=0 n=23`; prime:
  `prime_search` 16, 18, 20 bits, seed 59297. Four planted targets per
  curve, `target_seed 20261004`; one round, no warm-up; `max_trials` at
  its default; timeout 240 s per run. Exhausted and timed-out runs are
  kept as their own status.
- Koblitz arms: the polynomial-basis subspace base at dimensions 6 and 8,
  each under `subtract`, `mitm:m=2`, and `descent-algebraic:m=2` solved
  by `sat-cdcl`, `buchberger-f2`, `f4-f2` and `crossbred-f2`. Every arm
  on a base decomposes the same targets over the same points (same
  `FB1` id), so yields are directly comparable within a base.
- Prime arms: the abscissa base at 16, 32 and 64 points under `subtract`
  and `mitm:negation_folded=1`. No algebraic oracle runs on a prime curve
  in `ecbench`.
- The orbit bases (`koblitz-orbit`) are not in this session: on two of the
  four curves they do not exist (`research/ecbench_all_candidates_20261003`,
  PROTOCOL-2), and their yields are already in the ledger from that
  session.

## Predictions

- **P1 (exact oracles agree on yield):** on one base and one target,
  `subtract`, `mitm:m=2` and every `descent-algebraic:m=2` arm find the
  same number of relations from the same trials, because each decides
  decomposability exactly: the yield column is identical across those
  arms on every (curve, target, base) cell where all verify. A solver
  that returns a budget-exceeded verdict breaks this and is reported as
  such.
- **P2 (yield falls with the field, not the target):** at fixed base
  dimension, mean yield falls as `n` grows (the base is a fixed fraction
  `2^d / 2^n` of the abscissae), and within one curve the per-target
  yields of one arm have a coefficient of variation below 0.5.
- **P3 (solver cost per decomposition):** per solver call, SAT conflicts
  and Buchberger monomial operations both rise with `n` at fixed
  dimension; `f4-f2`'s solving degree does not exceed the semi-regular
  degree the record carries for the system's shape.
- **P4 (prime base size):** on each prime curve the yield rises with the
  base size roughly in proportion (16 → 32 → 64 points), and lookups per
  relation fall for the folded table.
- **P5 (correctness and replay):** every verified execution recovers the
  planted scalar; every replay in the audit is identical, the algebraic
  arms included (their solver budget is zero, so they are deterministic
  by the harness's rule).

## Success, stop, inadmissible

- The session succeeds when P5 holds and every arm has a status on every
  workload. Predictions are scored as written whether they hold or not.
- Inadmissible after seeing results: changing curves, seeds, targets,
  bases, oracles, solver parameters, the unit or the reference; dropping
  an arm; editing a session.
- Wall time is reported at whatever level the host earns and is
  descriptive; counts and yields are the result.

## Classification

Fixed in advance: no method is changed, so nothing here is an *advance*
or *engineering*. The session adds measurements of existing methods and
is classed *accounting* on the scoreboard if it changes any figure there,
which it is not expected to.
