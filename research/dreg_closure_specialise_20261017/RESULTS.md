# Results: uninformative for H_sum vs H_max; one fixed variable usually brings the ℓ = 7 draws back to degree 4

Scored natively by `dreg_score specialise runs/m3-*.jsonl` (`score-output.txt`)
against `PREREGISTRATION.md` (commit `9e0accfc4`).

- **Records.** Every value comes from a committed record in `runs/`.
- **Run.** 2026-10-17, 02:11–05:38 UTC: 32 processes, all exited 0.
- **Interruptions, disclosed in `runs/queue.log`.** The queue was resumed
  twice: once after the session's 2-hour background-job limit stopped it,
  and once after a container restart. Each time, the run in flight left an
  empty output and was rerun from scratch with the same binary
  (`runs/binary.sha256`, verified on restart).
- **Peak RSS:** 6.5 GB.
- **Run report** with graph and PDF: `reports/run-20261017T0211Z/`.

## The measurement

Without specialisation, none of these 8 draws refutes at degree 4 (closure-at-scale study).

| cell | draw | A0 (6,7,7) | A1 (6 coset,7,7) | B (6,6,6) | C (null) |
|---|--:|---|---|---|---|
| m3 (19,7) | 0 | refutes | does not | refutes | does not |
| m3 (19,7) | 3 | refutes | refutes | refutes | refutes |
| m3 (19,7) | 5 | does not | does not | refutes | does not |
| m3 (19,7) | 7 | refutes | refutes | refutes | refutes |
| m3 (21,7) | 0 | refutes | refutes | refutes | refutes |
| m3 (21,7) | 2 | refutes | refutes | refutes | refutes |
| m3 (21,7) | 3 | does not | refutes | refutes | refutes |
| m3 (21,7) | 6 | refutes | refutes | refutes | refutes |

- **Positive control B:** refutes on 8 of 8, which is valid.
- **Null control C:** refutes on 6 of 8. It is **not valid**, because ≥ 7 of 8 had
  to fail.
- **Registered verdict: uninformative for H_sum vs H_max.**
  - Arm A refutes on 12 of 16 (A0 6 of 8, A1 6 of 8).
  - Because the null control refutes just as often, that cannot be read as
    support for H_sum.

## What the runs do show (descriptive, not registered)

- **The ℓ = 7 degree-4 failure is marginal.**
  - Fixing **any single variable** lets degree 4 refute on 6 of 8 draws.
    That holds for a summand coordinate at 0 (A0), the same coordinate at 1
    (A1), and an intermediate-point coordinate (C) alike.
  - So the step from ℓ = 6 to ℓ = 7 does not put these systems far past
    degree 4. They sit at its edge.
  - Neither "max ℓ" nor "Σℓ" alone predicts which single fixes work. The
    intermediate-coordinate fix, which changes neither, works as often as
    the others.
- **Draws differ.**
  - (19,7) draw 5 resists every single fix and refutes only under B.
  - The other seven refute under at least two of the three single fixes.
- **Both branches of one guess.**
  - Both values of variable 6 refute (A0 and A1) on 5 of 8 draws. On those
    draws, a one-variable guess (2 branches) refutes the system at
    degree 4.
  - The degree-4 matrix of each branch has `C(N−1, ≤4)` columns, about
    92k–112k here.
  - Without the guess, the system needs degree 5 or more, with
    `C(N, ≤5)` ≈ 760k–1M columns.
  - This is a **hybrid (guess-and-closure) observation**. It is measured
    only at ℓ = 7, N = 40 and 42, with a single guessed variable.

## What this means for ECC2K-130 (derived)

- **The surviving explanations.** They remain "ℓ = 7" and "m·ℓ ≥ 21", from
  the closure-driver study. This study does not choose between them.
- **A third option enters the cost model: guess-then-degree-4.**
  - With k guessed variables, the oracle per target is
    `2^k · C(N−k, ≤4)^w` instead of `C(N, ≤5)^w`.
  - At ℓ = 7 a single guess sufficed for 5 of 8 draws.
  - How many guesses degree 4 needs as ℓ grows to the priced 33–45 is
    unmeasured. If it grows like the excess over ℓ = 6 per summand, it is
    about `m·(ℓ−6)` guesses: ≈ 108 at m = 4, ℓ = 33. That is far worse
    than degree 5.
  - So this does not move any row. It names the measurement that would:
    the number of guesses degree 4 needs at ℓ = 8, 9.
- **No row beats rho.** The free-oracle floor still binds.

## Successor

- **Guesses needed at ℓ = 8.** How many fixed variables does degree 4 need
  at m = 3, ℓ = 8?
  - It needs N ≥ 43 with unsatisfiable draws; for example (22,8) at
    N = 46, S = −2, which is 179k columns.
  - That sits at this container's memory limit.
- **The two surviving explanations.** Telling them apart still needs
  m = 4, ℓ = 6 on a larger host.
