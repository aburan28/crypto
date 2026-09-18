# Crossbred α AutoLab

Work in the repository root. The incumbent is frozen: chained Crossbred
(X1) and symmetrised Crossbred (X3) each failed to fit `α` because they
could not supply four usable `m = 3` rungs. Your job is to run the
harness beats, then only if a beat opens a new rung or a growing FFD,
record it in the route-target note and the scoreboard.

Do not claim a rho crossover. Do not GPU the search while `filters = 0`.
Do not change T4, `Q_enum`, or the paired-oracle divisor with hindsight.

## Iteration

1. Read `AGENTS.md` and `RESEARCH_ECC2K130_ROUTE_TARGETS.md` X1–X5.
2. `python3 research/ecc2k130_crossbred_autolab_20260918/crossbred_autolab.py plan`
3. Launch one beat. Inspect `runs/<id>/artifacts/claim.json`.
4. A replay that matches the freeze is accounting, not progress.
5. `fit.alpha` may not invent rungs and may not mix X1 with X3.
   `x5.ffd_chained_m4` is the next measurement: FFD growing with `n`
   at chained `m = 4` would be a result (H1). The chained *symmetrised*
   `S₃` at `m = 4` is still missing; do not relabel chained-`x` FFD as
   that arm.
6. `promote` copies a passing `runs/` claim into `evidence/`. Do not
   cite a gitignored run.

## Unit

X1/X3: 64-bit word XORs. `Q_enum = C(|F|, m−1)` at one word-op per pair.
X5: first-fall degree, 16 draws for a promotion (4 is smoke).
Field-product `S` stays on the product-law table and is not computed here.
