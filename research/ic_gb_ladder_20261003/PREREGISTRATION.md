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
