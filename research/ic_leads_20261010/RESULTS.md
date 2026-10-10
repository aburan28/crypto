# Lead ladder: first baseline run, 2026-10-10

Protocol: `docs/ic/LEAD_LADDER_METHODOLOGY_20261010.md`. Raw reports and exact
commands: `runs/20261010T19*/`. Tables as the summariser prints them:
`RESULTS_TABLES.md`. Binary `target/release/ic` at commit `7ab28103d`
(sha256 in each run directory). Host: Apple M4 Pro, load average 78–150
during the runs, so **no wall-clock figure here is evidence**; every ratio
below is read from counted group operations or native counters.

## What was run

| rung | curve | `log₂ r` | base | points | configurations | repeats |
|---:|---|---:|---|---:|---|---:|
| B0 | `K_0`, n = 13 | 11.0 | invariant subspace, dim 12 | 4,005 | fold / no fold × `mitm` (the `mitm-frobenius` arm was mis-named in this run and skipped) | 8 |
| B1 | `K_1`, n = 17 (`K_0` has composite `r`) | 16.0 | invariant subspace, dim 8 | 205 | fold / no fold × `mitm-frobenius` / `mitm` | 8 |
| B2 | `K_0`, n = 23 | 21.0 | invariant subspace, dim 11 | 2,025 | fold / no fold × `mitm-frobenius` / `mitm`; plus prefix subspace dim 11 (L4) | 8 |

All 8 × (2 + 4 + 4 + 3) = 104 runs recovered and verified the planted
logarithm; none exhausted its trial budget. The rho reference on each
instance is the counted matched walk (`A = 2` at n = 13, 17 where the
bench's rule picked the negation walk; `A = 2n = 46` at n = 23), 16 runs.

## L1. Signed Frobenius folding: landed, as predicted, engineering

Ratios are *unfolded control / folded baseline* at the same base and oracle
(`mitm-frobenius` at n = 17, 23; `mitm` at n = 13):

| rung | columns | relations needed (rows) | LA row ops | pair-table build adds | decomposition adds | trials per relation |
|---:|---:|---:|---:|---:|---:|---:|
| n = 13 | 12.9 | 6.1 | 5.7 | 1.0 (same `mitm` table) | 1.42 | 1.02 |
| n = 17 | 14.7 | 8.5 | 8.3 | 11.7 | 1.8 | 1.49 |
| n = 23 | 22.5 | 16.9 | 18.3 | 21.6 | 4.7 | 1.00 |

Readings against the pre-registered prediction (columns `÷ 2n`, rows
`÷ 2n`, LA `÷ (2n)²`, decomposition unchanged):

1. **The fold on columns is `n`, not `2n`, relative to the negation-folded
   control** (12.9 ≈ 13, 14.7 ≈ 17 less torsion merges, 22.5 ≈ 23): the
   control already folds `±P`, so the Frobenius adds the factor `n`.
   Against one column per signed point the fold is `2n` (`points_per_column`
   reads 25.8, 29.3, 45 against 26, 34, 46).
2. **Relations needed scale with columns but a little below them** (6.1, 8.5,
   16.9 against 12.9, 14.7, 22.5): the folded matrix needs more rows per
   column because orbit rows are denser and dependent rows occur; this is the
   realised relation saving, and it is what the ledger's `2n` should be read
   as at these sizes.
3. **LA is linear, not quadratic, in rows at this size** (ratio ≈ rows ratio,
   not its square): incremental Gauss on a matrix with ≤ 468 rows is far from
   the `(2n)²` regime; the quadratic saving is a statement about matrices
   that do not fit here.
4. **The pair-table build is the whole cost, and the Frobenius fold cuts it
   by `2n` exactly** (11.7 at `2n = 34` minus merges, 21.6 at `2n = 46`): the
   folded table stores one entry per orbit of pairs. At n = 23 the table is
   47,565 additions folded against 1,027,181 unfolded, and the decomposition
   phase is 424 additions. So on the `ic bench` ladder the fold is a
   **setup** saving first and a relation saving second.
5. **Decomposition is not unchanged**: the unfolded control needs
   1.4×–4.7× more additions in the decomposition phase because it needs
   `n`-fold more relations, each costing the same; per relation the cost is
   unchanged (trials per relation 1.0–1.5, within the binomial spread of
   27 relations).

Classification: **engineering** (`AGENTS.md` §3); the floor already credits
`A = 2n`. `B·D²/N` with `D` the online cost per target is 1.2·10³ at n = 23
on the folded baseline against 7.1·10³ unfolded: both far above 1, generic
territory, as expected.

## L1 twin: Frobenius-folded table versus negation-folded table on the same folded base

`mitm` (negation fold only) against `mitm-frobenius` at the same
`koblitz-orbit` base: setup adds × 11.7 (n = 17) and × 21.6 (n = 23), nothing
else moves (columns, rows, trials per relation within 4 %). This is the
clean measurement of what the *table* fold buys apart from the *column* fold:
`2n` on the table build and zero elsewhere.

## L4. Orbit-invariant base versus a prefix subspace of equal dimension (n = 23)

`binary-subspace:dimension=11` (prefix, not Frobenius-stable) against
`koblitz-orbit:divisor=1` (the invariant subspace), both with the
negation-folded `mitm` table because the Frobenius table refuses a
non-invariant base:

| metric | invariant subspace | prefix subspace | ratio |
|---|---:|---:|---:|
| signed points | 2,025 | 2,079 | 1.03 |
| columns | 45 | 1,040 | 23.1 |
| trials per relation | 2.61 | 4.28 | 1.64 |
| relations needed | 26.5 | 362 | 13.7 |
| LA row ops | 251 | 3,150 | 12.5 |
| table build adds | 47,565 | 1,081,000 | 22.7 |

The prefix base loses the whole fold (columns `|F|/2`), as predicted, **and**
its targets decompose 1.64× less often per trial at equal point count. The
second effect was not predicted; the likely cause is that the prefix
subspace's points are less uniform across the cofactor classes than the
invariant subspace's (the exact counting ceiling of `ic_boundary` would
decide this and is the next thing to read off the report). On the
challenge curve no invariant subspace exists, so the orbit-union base is
the only way to recover the first effect; the second effect is a property
of the base's point distribution and should be re-measured on
`compact-orbit-scan` before it is attributed to invariance.

## Reference rows

| rung | IC baseline `S` | rho `S` (matched) | IC / rho | rho plain `A = 1` |
|---:|---:|---:|---:|---:|
| n = 13 | 3.69 | 2.67 (`A = 2`) | 1.4 | 19.2 |
| n = 17 | 6.75 | 1.65 (`A = 2`) | 4.1 | 5.75 |
| n = 23 | 33.9 | 1.06 (`A = 46`) | 32 | 2.71 |

`S` here includes the measured-price GAE of the base build and LA (see the
methodology's calibration caveat), so the IC column is host-dependent in its
non-addition part; the ordering is not in doubt. The ladder loses to rho at
every rung and the gap widens with `n`, which is the known picture for the
`m = 3` pair-table shape at fixed subspace dimension; the point of the
baseline is the per-lead ratios above, not this column.

## Not measured in this run

- L2 (canonicalisation cost): the counted oracle reports `canonicalisations`
  (62 per 27 relations at n = 23, one per probe) but `ns_per_canon` is
  `null` in the calibration, so its share is uncharged and unknown. The
  first L2 step is to price it.
- L3, L5, L6, L7: no toggle exists yet (methodology §4).
- The orbit-union arm (`compact-orbit-scan`, the only arm that reaches
  n = 131) was not run: `ic bench --list` does not expose it on this build,
  and the compact-orbit pipeline's baseline is the landed `K = 600`
  configuration in `BOUNDARY_TARGETS.md`.
- The n = 23 floor-curve twin (L8): the CM-route construction
  (`e2_floor_curves_polclass.sage`) was still running when this file was
  written; the Vélu-built T23 family already exists in the parallel
  workspace and E1 measured the pipeline on it.
