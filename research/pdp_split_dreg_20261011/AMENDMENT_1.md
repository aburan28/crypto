# Amendment 1: solution read-out defect; grid rerun as `runs-v2/`

Written 2026-10-11 after the first measured grid (`runs/`, binary in
`runs/BINARY.sha256`, commit in `runs/COMMIT`) and **before** its rerun.

## What failed

Correctness gate 2 failed on cell `3:0:13:5`, plan `split-2+1`, `g = 0`.
Three of four targets (`x_R` = 4035, 7695, 4691) returned no decompositions.
Exhaustive enumeration found 1, 3 and 1 decompositions there, and the same
targets at `g = 2` and `g = 4` found all of them. No other run in the grid
failed the gate.

## Cause

F4 was correct. Its basis does not contain `1`, and every ordering of the
true decomposition (with its `u`) satisfies every basis element. The fault
was the read-out. `pq_f4_f2::solutions_from_reduced_basis` assumes that
every linear element's tail holds free variables only. Here 24 linear
elements had 99 tail terms on other elements' leading variables, so the
function's completion step produced wrong points. It then filtered all of
them out and returned `Some([])` instead of `None`.

The diagnostic is the scratch reproduction in the PR description, and the
regression test is
`pdp_split::tests::non_interreduced_basis_still_yields_points`.

## Change

`pdp_split::solutions_from_basis` replaces that call. It row-reduces the
linear elements first, enumerates the remaining free variables, and keeps a
point only if every basis element vanishes on it. That is exact for any
generating set. F4, the systems, the grid, the seeds, the budgets and the
pricing are unchanged.

## Consequence for the protocol

- Solving degrees, word XORs and F4 wall time do not depend on the read-out.
  The rerun must reproduce the first grid's solving degrees and XOR counts
  exactly. `pdp_split_analyze --replica runs` checks this.
- Verified sets and gate 2 are taken from `runs-v2/`.
- `runs/` stays as recorded, including its failed cell.

The rerun is `run_grid.sh` with `runs` replaced by `runs-v2`.
