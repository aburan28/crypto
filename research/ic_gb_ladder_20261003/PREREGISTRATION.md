# Does the refutation degree keep rising past the sharp first-fall bound? An external-engine `m = 3` ladder

Registered before any registered cell runs. §2 lists everything that ran before
registration and what was seen.

## 1. Why

- **The crux the program keeps returning to.** Every subexponential small-characteristic
  ECDLP claim rests on the first-fall-degree assumption, `D_reg = D_ff + o(1)`
  (Petit–Quisquater 2012; survey of the stalemate in
  [RESEARCH_ECC2K130_IC_LITERATURE.md §1–2](../notes/ecc2k130/RESEARCH_ECC2K130_IC_LITERATURE.md)).
  Kousidis–Wiemers (JMC 2019, Thm 3.2) prove `D_ff ≤ m² − m + 1` for the Weil descent of
  `S_{m+1}` over `F_{2^n}`, sharp at **7 for `m = 3`**. If the solving degree tracked that
  bound, the per-target solve would be polynomial in `n` and the decomposition route would
  be alive (survey §2: `c → 0`). If instead the solving degree grows with `ℓ = dim V`, the
  route is dead at every scale the growth persists. The repository's measurements all read
  growth, but at `ℓ ≤ 4` (`x4`) and `ℓ ≤ 5` (`rr`), where the gap over 7 is at most 2:
  [ic_rr_norm_ladder_20260930](../ic_rr_norm_ladder_20260930/RESULTS.md) read the direct
  `S₄` descent (`x4`, `3ℓ` unknowns, Boolean degree 6) at 6, 8, 9 for `ℓ = 2, 3, 4`
  and the Riemann–Roch norm form (`rr`, `3ℓ + 1` unknowns, degree 4) at 4, 5, 6, ≥ 8 for
  `ℓ = 2 … 5`.
- **Why it stopped there: the engine.** The in-tree instrument is a sparse Macaulay
  scan (`solving_profile_sparse`). At `ℓ = 5` (16 unknowns) it did not resolve degree 8
  in 6,000 CPU-s; at `ℓ = 6` it scanned to 6 only. The ladder ended two rungs short of
  anything that separates "`+o(1)`" from "grows".
- **What changes here.** The same systems, byte-identical from the same builders and
  seeds, are handed to an **external, independent Gröbner engine**: Singular 4.3.2's
  `slimgb` on the homogenisation, truncated at degree `D` (`degBound`). For a
  homogeneous ideal the truncated basis spans exactly the degree-`D` Macaulay row space
  (Lazard), so "`h^D` in the truncated basis" is exactly "`1` in the Boolean Macaulay matrix
  at degree `D`" — the same observable as the in-tree ladder, from different code. In the
  smoke runs (§2) it resolved the `ℓ = 5` `rr` cell in 67 s where the in-tree scan
  failed at 6,000 s, and it reads **below the system degree**, which the in-tree scan
  never tried (its loop starts at the system degree). That floor hid two things (§2):
  `x4` at `ℓ = 2` actually refutes at 5, and some `rr` draws are refuted at degree 1 by a
  **constant equation**, the trace condition `Tr(e₁) = 0` on a draw with `V ⊂ ker Tr`
  and `Tr(x_R) = 1`. Those read as 4 at the floor and were mistaken for degree-4
  refutations.
- **What is being asked.** Over `ℓ = 5, 6, 7, 8` in the regime `n ≥ 3ℓ`, does the
  refutation degree of `x4` keep rising by about one per rung (law `ℓ + 5`, excess over
  the sharp bound growing 3, 4, 5, …), or does it stop (an excess that stays put across
  two consecutive resolved rungs)? Same for `rr` (law `ℓ + 2`). A plateau at these sizes
  would be the first measured saturation signal in the program and would reopen the
  algebraic route's exponent audit; continued growth extends the measured boundary by two
  to three rungs and puts the first-fall assumption 4–5 degrees behind the data on the
  canonical object.
- **The prior is growth.** Every ladder so far has read about one degree per unit `ℓ`,
  and the literature's own refutation (Huang–Kosters–Yeo) argues the assumption fails.
  This round is cheap enough that the asymmetric payoff (a plateau would matter far more
  than another rung of growth) justifies it, and the independent engine audits the
  earlier readings for free.

## 2. What ran before registration (disclosed; not evidence)

Smoke runs on this container, 2026-10-03, with the `--dump-dir` export added to
`examples/rr_degree_ladder.rs` and the driver [refute.sing](refute.sing):

- **Calibration on `K₁/2¹⁷`, draw 0** (the registered ladder's own draws, same seed):
  `rr` refuted at 4, 5, 6 at `ℓ = 2, 3, 4` — identical to the in-tree readings; `x4` at
  8, 9 at `ℓ = 3, 4` — identical; `x4` at `ℓ = 2` reads **5** (in-tree: 6, its floor).
  `std` with `degBound` gave the same degrees, 10–30× slower than `slimgb`.
- **All four draws of `K₁/2¹⁷ ℓ = 2`, all arms** (the plumbing check `SMOKE=1`): `rr`
  4, triv, 4, 4; `x4` 5, 5, 5, 6; control pinned at 5 (three draws) and unresolved to 12
  (one). "triv" is draw 1, whose `rr` system contains the constant equation `1`.
- **Constant equations.** Of the `rr` systems dumped so far (`K₁/2¹⁷ ℓ = 2…6`,
  `K₀/2¹³ ℓ = 5`, `K₁/2¹⁹ ℓ = 5, 6`, `K₀/2¹⁹ ℓ = 6`), three contain a constant
  equation: `K₁/2¹⁷ ℓ = 2` draw 1, `K₁/2¹⁷ ℓ = 5` draw 1 and `K₁/2¹⁹ ℓ = 5` draw 0.
  The last two are exactly the "degree-4 refutations at `ℓ = 5`" that
  [ic_rr_norm_ladder_20260930/RESULTS.md](../ic_rr_norm_ladder_20260930/RESULTS.md)
  attributed to the Kosters–Yeo regime. That sentence is wrong and gets an erratum with
  this round's results; the cell medians it fed are unchanged (4 was never the median).
- **Reach.** `rr` at `K₁/2¹⁷ ℓ = 5` draw 0 (16 unknowns): refuted at **8** with
  `degBound = 9`, 67 s, under 1 GB. `x4` at the same draw (15 unknowns) with
  `degBound = 11` in one shot: **out of memory at 4.5 GB after 349 s** (no reading).
  `rr` at `K₁/2¹⁷ ℓ = 6` draw 0 (19 unknowns) with `degBound = 10` in one shot: **out of
  memory at 4.5 GB after 426 s**. Both one-shot runs compute every degree up to the bound,
  which the registered scan (one process per `D`, stopping at the first refutation) does
  not, so they overstate the memory the scan needs; memory, not time, is the limit, which
  is why §3 has a second phase.
- **Macaulay2 1.22** (`gb` with `DegreeLimit` on the same homogenisation) was tried on
  `K₁/2¹⁷ ℓ = 4` `x4` and killed after 600 s where `slimgb` took 8.6 s; it is not used.

Nothing above is cited as a result. Every registered cell is run again from scratch by
[run.sh](run.sh), including the calibration cells.

## 3. Object

- **Systems.** Exactly those of the earlier ladder, from the unchanged builders in
  `examples/rr_degree_ladder.rs` (`build_rr`, `build_direct_x_system` with `m = 3`,
  `random_control_system` of `rr`'s shape), same seed `20260930`, same draw sequence,
  four rootless draws per cell (`--unsat 4`); exported with `--dump-dir` as Singular
  ideals. A draw is measured on an arm only when that arm's root count is 0 (`x4`: field
  roots in `V³`; `rr`: ordered base-point triples summing to `R`).
- **Observable.** For each (cell, draw, arm), the least `D` at which the truncated
  homogeneous `slimgb` basis contains `h^D` (**refutation degree**), scanned upward one
  process per `D` from 3 (`rr`, `ctrl`) or 4 (`x4`) to a cap of 12 (`rr`, `ctrl`) or 14
  (`x4`). A satisfiable control that is never refuted is **resolved by pinning** at the
  first `D` whose basis has a leading term `x_i h^k` for every variable (the in-tree
  `resolves()` rule); it is reported as `Dp`. A system refuted at degree 1 contains a
  constant equation and is reported as **`triv`**.
- **Cells**, in the regime `n ≥ 3ℓ`, `n` prime, the smallest the tooling builds:

  | cell | `n ≥ 3ℓ` | `rr` unknowns | `x4` unknowns | role |
  |:--|:--|--:|--:|:--|
  | `K₁/2¹⁷ ℓ = 2, 3, 4` | 17 ≥ 12 | 7, 10, 13 | 6, 9, 12 | calibration against the in-tree ladder |
  | `K₁/2¹⁷ ℓ = 5` | 17 ≥ 15 | 16 | 15 | first rung past the in-tree reach |
  | `K₁/2¹⁹ ℓ = 6`, `K₀/2¹⁹ ℓ = 6` | 19 ≥ 18 | 19 | 18 | the discriminating rung for `x4` (law says 11) |
  | `K₁/2²³ ℓ = 7` | 23 ≥ 21 | 22 | 21 | reach |
  | `K₁/2²⁹ ℓ = 8`, `K₀/2³¹ ℓ = 8` | 29, 31 ≥ 24 | 25 | 24 | reach; `x4` expected to censor |

  `K₀` does not exist in the tooling at `n = 17, 29`; `K₁` not at 31.
- **Budget and censoring.** Each Singular process: 3,600 CPU-s (`ulimit -t`) and 4.5 GB
  of address space (`ulimit -v`); three lanes pinned to CPUs 1–3, cheap cells first,
  `rr` before `x4` before `ctrl` within a cell. A killed process censors the draw at
  `≥ D`; a process that finishes unrefuted and unpinned advances to `D + 1`. **Phase 2:**
  after the lanes finish, every draw censored by memory (Singular's "no more memory") is
  retried **once**, alone, with 12 GB, at the degree it died at and upward under the same
  CPU limit. Censoring is never negative evidence. No other retries.
- **Isolation.** Degrees are not timings; the lanes run without the benchmark lock.
  Singular's own `ms` and the wall per process are recorded and reported as advisory.

## 4. Analysis (`examples/gb_ladder_analyze.rs`, fixed now)

1. Per draw: the reading (`D`, `Dp`, `≥D`, `triv`), the summed wall, and the in-tree
   reading where the draw is shared.
2. Per cell and arm: the readings, and the **median** by the earlier ladder's rule
   (majority exact → median of the exact readings; otherwise `≥ min`), with `triv` draws
   **excluded** from medians, slopes and pairs (they are reported, and their `x4`
   readings are kept as an observation on the same draws).
3. Per arm and curve: the OLS slope of exact cell medians against `ℓ`.
4. Paired `rr − x4` on draws exact on both arms.
5. **Excess** of the `x4` median over the sharp first-fall bound 7, per rung.
6. Calibration: agreement with the in-tree reading on every shared draw where both are
   exact and the in-tree reading is above its floor.

## 5. Predictions (pass/fail)

- **P1, calibration.** On shared draws with both readings exact and the in-tree reading
  above its floor (`rr` > 4, `x4` > 6), the two engines agree on **every** draw. On
  `K₁/2¹⁷ ℓ = 2`, `x4` reads 5 on at least three of four draws.
- **P2, growth of `x4`.** The `x4` median is exact and ≥ 10 at `K₁/2¹⁷ ℓ = 5`, and, if
  exact at `ℓ = 6` on either `n = 19` curve, ≥ 11 there. The excess over 7 then reads
  ≥ 3 and ≥ 4.
- **P3, growth of `rr`.** `rr` medians read 8 at `ℓ = 5` and 8 or 9 at `ℓ = 6`
  (both `n = 19` curves); if `ℓ = 7` resolves, 9 or 10. No two consecutive resolved rungs
  on the same curve are equal.
- **P4, control.** No random control is refuted below the `rr` reading of its draw; the
  controls resolve by pinning or censor, as before.
- **P5, pairs.** `rr` is below `x4` on every draw exact on both, by 1 to 3.

## 6. Decision rule, fixed now

- **Growing**: P2 and P3 hold on the resolved rungs (at least `ℓ = 5` for `x4` and
  `ℓ = 5, 6` for `rr`). Then the solving degree is measured 3–4 degrees above the sharp
  first-fall bound on the canonical object, still rising, at `n = 17…19`; the stop
  decision's item 3 stays closed and the first-fall assumption is behind the data by
  that margin. No exponent is claimed in either direction.
- **Saturation signal**: for either arm on one curve, two consecutive resolved rungs
  with `ℓ ≥ 5` whose medians are equal or decreasing, with P4 holding. Then the
  algebraic route's exponent audit is reopened as a Coordinator decision, and the next
  round is the same ladder one rung further, on a second curve, before anything else.
- **Inconclusive**: fewer than two new rungs (`ℓ ≥ 5`) resolve for `x4` and fewer than
  three for `rr`. Then the reach is reported with its censoring and the next step is an
  engine question, not a mathematical one.
- A `triv` majority in any cell drops that cell from the fit and is reported.

## 7. Scope

Two Koblitz curves, `n` 17–31, `ℓ` 2–8, `m = 3`, four rootless draws per cell, one
external engine (Singular 4.3.2 `slimgb`, degree-truncated on the homogenisation). The
observable is the Macaulay refutation degree of one presentation of the ideal; it is a
stage diagnostic (AGENTS.md §2, §5): no end-to-end cost, no yield, nothing at `n ≈ 83` or
131, and nothing asymptotic. At `m = 3` no reading can tie rho even with a free oracle
(survey §2); the reading bears on `m ≥ 4` only through the assumption's fate.

## Amendment 1 (2026-10-04, written while the first `ℓ = 6` retry is still running): CPU budget of the retry phase

Phase 1 is complete: every `ℓ ≥ 6` run and every `ℓ = 5` `x4` run was censored by the 4.5 GB
lane limit. In phase 2 (12 GB) the four `ℓ = 5` `x4` draws resolved (10, 10, 10, 10), and
`K₀/2¹⁹ ℓ = 6` draw 0 `rr` finished degree 7 unrefuted in 535 s and has been running degree 8
for 52 minutes against the 3,600 CPU-s limit. The limit that binds the decisive rung is now
CPU time, which was set with the in-tree engine's costs in mind, not this one's. This
amendment is written before that process ends and before any `ℓ = 6` reading exists.

- **Change.** Phase-2 retries of the `rr` and `x4` arms at `ℓ ≤ 6` run under 14,400 CPU-s
  per process (four hours), still at 12 GB and still once per draw. A phase-2 process that
  was already killed by the 3,600 s limit on one of those arms is re-run once under the new
  limit (`PHASE=3` in `run.sh`); every other limit, cell, arm and rule is unchanged. The
  controls and the `ℓ ≥ 7` cells keep 3,600 s.
- **Why it cannot bias the reading.** A CPU budget decides only whether a process finishes;
  the refutation degree is a property of the system and the engine's truncation, so a larger
  budget can turn `≥ 8` into an exact reading, never move an exact reading.
- **Order.** Retries run `rr` at every `ℓ`, then `x4`, then the controls, smallest `ℓ`
  first (the first phase-2 launch was ordered by `ℓ` and had reached an `ℓ = 5` control,
  which was taking 17 minutes per degree; it was stopped during that control's degree-8 run
  with no line written, and re-ordered; that run is repeated from degree 8 when its turn comes).

## Amendment 2 (2026-10-04, after the first `ℓ = 6` retry): one memory death per cell and arm

`K₀/2¹⁹ ℓ = 6` draw 0 `rr` ran degree 8 for 55 minutes under Amendment 1 and died at the
12 GB limit. The remaining memory-killed draws of that cell and arm would repeat the same
death after about an hour each, adding at most one degree to a lower bound (`≥ 8` for
`≥ 7`) and no exact reading.

- **Change.** Within phase 3, once a retried draw of a (cell, arm) dies by memory at 12 GB,
  the remaining retries of that same (cell, arm) are not run; those draws keep their
  phase-1 bound. Every (cell, arm) still gets at least one retry at 12 GB, so each rung
  and arm is tested once at the larger limit. Nothing else changes.
- **Why it cannot bias.** A skipped retry leaves a censored reading censored. It can
  neither create nor move an exact reading, and the medians are formed by the registered
  rule from the bounds that exist.

## Amendment 3 (2026-10-04, after the second `ℓ = 6` retry): the controls' retries

`K₁/2¹⁹ ℓ = 6` draw 0 `rr` died at 12 GB at degree 8 after 44 minutes, as `K₀/2¹⁹`'s had;
the `ℓ = 6` rung is beyond this engine on this host for both curves. What remains of the
retry phase is dominated by the four `ℓ = 5` controls, which keep the 3,600 s limit and take
about 17 minutes per degree each, and which can only add bounds (a control is informative
only through P4, "not refuted below `rr`", and `rr` at `ℓ = 5` is 8 while every control
there already reads `≥ 8`).

- **Change.** Amendment 2's one-death rule is extended, for the control arm only, to a death
  of either kind (memory or CPU): after the first retried control draw of a cell dies, that
  cell's other control retries are skipped and keep their phase-1 bound. The `rr` and `x4`
  arms are unchanged (Amendment 2 already governs them).
- **Why it cannot bias.** As in Amendment 2: skipped retries stay censored; no exact reading
  is created or moved, and P4 is evaluated on the bounds that exist.

## Amendment 4 (2026-10-04, after the first `ℓ = 5` control retry): the controls' degree cap

`K₁/2¹⁷ ℓ = 5` control draw 0 scanned to the registered cap of 12 under 12 GB without
dying: it pins 15 of its 16 variables from degree 8 on and never the last one, at about 20
minutes per degree. Amendment 3 ends a cell's control retries only on a death, so the other
three `ℓ = 5` controls would each spend two hours producing the same `≥ 13`.

- **Change.** In the retry phase, a control is scanned only up to degree 7 (one process per
  retried draw). Degree 7 is the degree below `rr`'s reading in every cell where `rr`
  resolved (8 at `ℓ = 5`), so an unrefuted degree-7 run is exactly what P4 needs
  ("not refuted below `rr`"), and nothing higher is used by any prediction or rule.
- **Why it cannot bias.** The controls enter the verdict only through P4, which asks for the
  absence of a refutation below `rr`'s degree; a cap at 7 tests exactly that and leaves the
  reading censored above it.

## Amendment 5 (2026-10-04, at the end of the retry phase): no control retries at `ℓ ≥ 7`

Every `rr` and `x4` retry has run, every `ℓ ≤ 6` control has had its retry, and the phase
had reached the `ℓ = 7` control, which was 31 minutes into degree 6 under 12 GB. The `ℓ ≥ 7`
cells have no exact `rr` reading, so P4 ("not refuted below `rr`") has nothing to compare a
control against there, and no other rule reads a control.

- **Change.** The `ℓ ≥ 7` controls are not retried; they keep their phase-1 bound (`≥ 6`).
  The `ℓ = 7` control's degree-6 process was stopped without a line being written.
- **Why it cannot bias.** As in Amendments 2–4: a censored control stays censored, and no
  prediction or rule at `ℓ ≥ 7` uses one.
