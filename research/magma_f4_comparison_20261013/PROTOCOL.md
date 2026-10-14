# Protocol: Magma's F4 on the systems behind our Gröbner comparisons

Registered 2026-10-13, before any Magma output exists.

## Question

Our claims against F4 (`research/pdp_root_finder_20261012/` §5: the
root-finders against `f4_fp_tower` on PKM's Kummer cells; the inherited-F4
oracle on binary decomposition systems; `C₄`, the `k = 4` Macaulay solve)
are all against this repository's own Gröbner engines.  msolve was
12–24× faster than our `f4_fp` on square systems
(`research/gb_cas_comparison_20261011/`).  Magma's F4 is the reference the
literature uses.  How long does Magma take on the *same* systems, and does
it reach the same verdicts?

## Systems (`scripts/`, exported from the committed rows, not regenerated)

- 26 PKM Kummer `tower` systems at `p₁ = 2013265921`: every finished
  `f4_fp_tower` row of the cells in the §5 table (`m = 2`, `t = 6…11`;
  `m = 3`, `t = 3…6`; `m = 4`, `t = 2…4`), two random targets a cell, by
  `pkm_tower_pilot --export-magma-rows` from the row files in
  `pkm_row_files.txt`.  Every one is a unit ideal in our runs (no
  decomposition).
- 5 `k = 4` systems: generic planted `d = 3…6` and the `S₅` system of
  `--s5 --p 269 --seed 1` system 0 (`gaudry_quartic_replay --export-magma`),
  with our rational solutions written into each script.
- 8 Boolean root systems of the `n = 17` decomposition search
  (`research/gb_cas_comparison_20261011/`), 35 variables, unit ideal by
  msolve and PolyBoRi.

Each script times `GroebnerBasis` (F4 is Magma's default over a prime field
and on Boolean rings) with `Cputime`, prints a `RESULT-GB` line, and for the
`k = 4` systems then times `Variety` and compares Magma's solution set with
ours (`RESULT-VAR`).

## Where it runs

The Magma online calculator (V2.29-10; 60 s a submission, 50,000 bytes of
input), submitted by hand in a browser: the calculator refuses scripted
clients, and none is used.  `make_bundles.py` packs the scripts into 14
submissions; heavy cells get their own.  A script cut off at 60 s is
recorded as "> 60 s minus the time already used in its bundle", not as a
result.

## Units and caveats fixed in advance

Magma `Cputime` on the calculator's (unknown) hardware against our recorded
times on a 4 vCPU Xeon @ 2.10 GHz.  Hardware differs: ratios within 3× are
called "comparable", not wins or losses.  Our PKM times are `f4_fp_tower`
runs from the committed rows (a single thread); our `k = 4` times are this
session's `solve_system4` (four threads).  Our Boolean times are a
decision search, Magma's a full basis: different computations, reported
side by side as in the 2026-10-11 study.

## Predictions

- M1. Every verdict agrees: unit ideal on all 26 PKM and 8 Boolean
  systems; the same rational solution set on the 5 `k = 4` systems.  A
  disagreement stops the comparison until explained.
- M2. Magma is at least 10× faster than `f4_fp_tower` on every PKM cell it
  finishes (msolve's margin over our `f4_fp` was 12–24×).
- M3. Magma does not finish the `m = 3, t = 6` and `m = 4, t = 4` cells in
  60 s (ours: 3–3.5 h and 55 min).
- M4. Magma finishes the `S₅` solve (basis and variety) within 60 s.
- M5. The root-finders' cost per target (`research/pdp_root_finder_20261012`:
  197 multiplications, microseconds) stays at least 10⁴× below Magma's time
  on every PKM `m = 3` cell Magma finishes.

## Decision

The paper-facing comparison quotes Magma's times wherever it finished, ours
beside them, and never our engine alone as "F4".  If M2 holds, every earlier
"against F4" ratio is restated against Magma where Magma finished.
