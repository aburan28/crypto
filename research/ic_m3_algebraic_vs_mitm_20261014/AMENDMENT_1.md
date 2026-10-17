# Amendment 1: the `m = 4` grid with the wide (u128) chain solver

Registered 2026-10-15. Committed before any run of the grid below.

## Why

In the original `m = 4` grid (PROTOCOL.md), F4, F6-IC and SAT were
**unsupported** at n = 23, 29 and 31. The chain system there has
`4d + 2n > 64` Boolean unknowns, and the narrow engines keep monomials in a
`u64` (`koblitz_groebner::MAX_VARS = 64`;
`problems/binary-boolean-system-64-variable-cap.md`).

The owner asked for the cap to be widened to 128 and the `m = 4` grid
rerun. The owner chose this route over rewriting the 64-bit engines (whose
counters sit behind frozen benchmarks): use the repository's existing
`u128` solver, `wide_groebner::wide_groebner_decompose`
(`MAX_WIDE_VARS = 128`). It works on the same chained-S₃ system. At each
search node it:
- branches on whole summands, over the coordinate vectors that are base
  points;
- linearises the substituted equations to a fixpoint (a Macaulay reduction
  at their own degree);
- closes the branch on a reduced `1`;
- lifts leaves in the group.

That solver is not inherited F4: its search is summand enumeration pruned
by linear algebra. So it is a **new arm**, `wide-chain`, run on **all
seven** `m = 4` cells. Its slope is then fitted on one solver throughout.

## What is frozen

Everything in PROTOCOL.md §2–§6 applies unchanged, except:

- **Arms:** `wide-chain` only. The other `m = 4` arms keep their original
  records in `runs/m4/`.
- **Cells:** `1:11:3`, `0:13:4`, `1:17:5`, `1:19:5`, `1:23:6`, `1:29:8`,
  `0:31:8`. The chain has `4d + 2n` = 34, 42, 54, 58, 70, 90 and 94 unknowns,
  all ≤ 128.
- **Targets:** seed `20261014`, a 1,024-target pool, and 6 decomposable plus
  6 non-decomposable targets. This is the same rule as before, so the pool
  and the selected targets are the same scalars as in `runs/m4/`.
- **Floor:** the floor is re-timed in this run, in `runs/m4-wide/truth-*.json`.
- **Budgets and early stop:** 120 s and 4 GiB per target; early stop after
  4 leading censored targets.
- **Output:** `runs/m4-wide/`.
- **Analysis:** `ic_algebraic_floor analyze --dir runs/m4-wide --out
  amendment-1/ANALYSIS.json`. Same estimator, bootstrap (2,000, seed
  20261014) and decision rule (§5).
- **Counted unit:** search nodes, recorded per arm and never compared
  across arms.
- **Prediction (before data):** `wide-chain` grows faster than the floor.
  Smoke solves on one non-pool target (`k = 7`) took 2, 38, 66, 530, 47,849
  and 57,587 ms at n = 11, 13, 17, 23, 29 and 31. The tail is close to the
  `|F|²` leaves of summand enumeration.
- **Follow-up:** if the verdict is "grows more slowly than the floor",
  PROTOCOL.md §7 applies.

The smoke solves above are not claim data.
