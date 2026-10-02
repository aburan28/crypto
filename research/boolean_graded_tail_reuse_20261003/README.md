# Graded Boolean Macaulay batch reuse

`PROTOCOL.md` and `protocol.json` freeze an exact algebraic experiment before
implementation or timing. A fixed quadratic core makes the cubic Macaulay
block identical across changing affine tails. The proposed candidate reuses
the high-block elimination schedule, then recomputes the changed low-degree
tail and materializes every output. It must beat fresh matched packed controls
in complete cold batches to qualify.

Status: **protocol only**. No producer, discovery, holdout, timing, solver
cost, relation yield or rho comparison exists. The new holdout seeds are
unused. Inputs are generated public Boolean systems, with no curve or key
interface. The previous full-trace screen in draft PR #1202 is context,
not an accepted performance baseline.
