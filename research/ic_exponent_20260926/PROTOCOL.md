# Ledger §20 protocol, v1: the Koblitz thread's exponent against batch rho

Declared 2026-09-26, before the pricer existed and before anything below
ran. The only computation made first is `predict.py`, which reads frozen
files that predate this protocol and is committed with it
(`prediction.json`).

## Question

At `k = 32` targets, with every phase priced on current `main`, how does
the Koblitz collection thread's `m = 3` method compare with batch rho at
the same `k`, and how does that ratio scale with `r`? Is there a size
where it crosses one?

§19 answered this at the thread's own sizes. At `r = 2^39` the closest
figure was `6.38×` batch rho, with every rho step and stored pair priced
(ledger §19.5). But it answered it on the commit those figures were frozen
at, and two things have changed since:

- **The prices.** `main` has made the Koblitz primitives faster since
  §19's commit (§19.8). The canonical step, the stored pair and the scan
  all cost less there, and none of them has been re-priced.
- **The descent.** It was priced per trial. The frozen ledger charges it
  as trials times the `1.19` units of one walked probe: 46 to 52 trials,
  or 54 to 62 units, for 32 targets at the headline. The holdout arm's
  total in `docs/ic/runs/koblitz-collection-aim-20260922.json` sits 7
  units below the main arm's, for 6 fewer trials. But an `m = 3` descent
  trial is a full scan of the base, stopping at the first witness in
  blocks of 1,024. So its cost per target is the probes it scans, and
  those do not amortise over the targets.

## Boundaries, stated before measuring

- **Floor, per target at `k = 32`:** Kuhn–Struik's batch share times
  the single-target generic floor, `L(32) · √(π/4n)` with
  `L(32) = Σ_{i<32} C(2i,i)/4^i / 32 = 0.19869`.
- **Reference:** batch rho at `k = 32` (§19, `rho_batch_with` on
  `SignedFrobeniusClasses`), run on **the same 32 targets** as the
  index-calculus run, in the same process.
  - Its counted group operations are priced at a canonical step: one
    batched addition plus the table-driven canonicalisation, the step
    §19.1 declared. The step is measured in the same process against the
    same unit.
  - Bailey et al.'s step is reported beside it as a model, because its
    step count is still unmeasured (§19.7).
- **Cold (secondary):**
  - The index-calculus run's shared phases plus one mean descent,
    derived from the `k = 32` run and labelled so.
  - Against single-target rho run on each of the first measurement set's
    32 targets separately.

## Sizes

Nine curves with a prime subgroup, spanning `2^18`–`2^47.2`.

| curve | `log₂ r` |
|:--|--:|
| `K_1/GF(2^19)` | 18.0 |
| `K_1/GF(2^23)` | 22.0 |
| `K_1/GF(2^45)` | 24.8 |
| `K_0/GF(2^37)` | 27.8 |
| `K_1/GF(2^43)` | 32.1 |
| `K_1/GF(2^47)` | 36.6 |
| `K_0/GF(2^41)` | 39.0 |
| `K_0/GF(2^53)` | 44.3 |
| `K_0/GF(2^61)` | 47.2 |

The `r` values come from `examples/koblitz_degree_census.rs`. The spacing
brackets both predictions' features: the law's crossing near `2^23` and
the model's minimum near `2^28`.

## The recipe, and how its free parameter is chosen

**Fixed at every size**, by the thread's own rules:

- `summands = 3`, `solver = pair_table`, pair-table tier `auto`;
- `collection_aim = true`;
- one planned unit, then units until the relations determine every column
  (`max_units = 100000`);
- 32 random known-answer targets;
- default linear algebra.

**Rules, not tuned per size:**

- **Window** `w = max(1, round(|F|/32))`, the rule of the thread's v5 and
  `n = 61` rungs.
- **Unit** `unit_trials = clamp(round(2r/(n|F|²)), 16, 65536)`, about
  sixteen units a run.
- **Descent cap** `max_trials = max(10^5, 64 · ⌈2.02 r/|F|²⌉)`.

**Swept on a seed set of its own, never on a measurement set:**

- The column count `c` (`|F| = 2n·c` points requested):
  - the grid is `prediction.json`'s `sweep_grid_columns`, the model's
    optimum times `2^{j/2}`, `j = −2 … 4`;
  - if the minimum falls on an edge, the grid extends one step at a time,
    at most three steps, and never below one column.
- The descent's summands, `2` or `3`. The table tier depends on it, so
  each `(c, m)` is a separate whole run.

Each size's recipe is the `(c, m)` with the least `S` per target on the
sweep set. It is fixed before any measurement set runs.

**Seed sets:**

- **Sweep `W`:** factor-base and workflow seed 101, targets
  `random_seed` 10100–10131.
- **Measurement `M1`–`M4`:** seeds 201–204, targets `100·seed` to
  `100·seed + 31`.
- **Rho:** batch rho on each set's own targets, with seed `0x200000 + seed`.
- **Cold rho:** on `M1`'s 32 targets, one walk each, seeds
  `0x210000 + i`.

## Accounting

**Unit.** One batched affine addition, `add_many` over 1,024 subgroup
points, on the same curve.

- It is measured in the same process immediately before and after each
  repetition.
- That repetition's conversion is the mean of the two measurements.

**Index calculus.** Every phase is timed on one thread, exclusively and
in the order the workflow runs them, and divided by the unit:

1. setup (the curve);
2. selection, including the projection the workflow charges to it;
3. pair-table build;
4. collection, planned and extension units alike;
5. relation verification;
6. linear algebra, every attempt, including its column-log check;
7. descent setup;
8. the descent of each target, including its recovery check;
9. the final `[d]G = Q` verification.

Target construction is excluded, because both sides receive `Q`.

Each phase's native counts are recorded beside its time:

- points, orbits and selection draws;
- stored pairs and tier;
- probes, summands scanned and relations;
- rejected and duplicate relations;
- solve attempts and the linear algebra's report;
- descent trials per target.

The conversion each count implies is reported, so that another host can
re-price the counts.

**Rho.** Counted group-addition equivalents (walk and setup) times the
canonical step's price in the unit, measured in the same process.

**Repetitions.**

- The index-calculus pipeline runs `R` times per set, each time
  rebuilding everything from nothing:
  - `R = 3`;
  - `R = 15` when one repetition takes under 50 ms.
- The figure is the median over repetitions.
- The first repetition is reported separately, since it carries the
  one-time costs.
- Counts must be identical across repetitions.

**`S` per target** is total units over `k·√r`. The ratio at a size is:

- the mean over `M1`–`M4` of the index-calculus `S`,
- over the mean of the priced rho `S`,
- with a 95% interval from the four per-set ratios (`t`, three degrees of
  freedom).

**Controls.**

1. **Replay against the workflow.** Every measurement run's counts must
   equal those of `ic workflow` on the same parameter file: per-unit
   relations, trials and summands scanned, stored pairs and tier,
   per-target descent trials, and recovered logarithms.
2. **The frozen headline, re-priced.**
   - `docs/ic/params/k0n41-least-on-u150.json` runs as it stands (seed 1,
     targets 100–131), with its baseline off.
   - Its counts are compared with the frozen run's
     (`koblitz-collection-aim-20260922.json`, arm
     `aimed_at_least_mentioned`): 197 relations, 3,450 trials, 883,200
     summands, 1,519,296 stored pairs and 52 descent trials. It is then
     priced like every other row.
   - This is the accounting row that links §20 to §19.5's `6.38×`.
3. **The thread's own recipes,** re-priced on `M1`–`M4` as diagnostic
   rows:
   - `n = 41`: `|F|` 15,300 requested, `w = 256`, `u = 150`, `m = 3`;
   - `n = 53`: `|F|` 15,000, `w = 477`, `u = 2048`, `m = 2`.

## Predictions (from `prediction.json`)

| curve | `log₂ r` | law | model | model's `(c, m)` | model's largest phase |
|:--|--:|--:|--:|:--|:--|
| `K_1/GF(2^19)` | 18.0 | 0.83 | 4.89 | 2, 2 | descent 53% |
| `K_1/GF(2^23)` | 22.0 | 1.20 | 3.42 | 5, 2 | selection 37% |
| `K_1/GF(2^45)` | 24.8 | 1.18 | 3.35 | 5, 2 | selection 39% |
| `K_0/GF(2^37)` | 27.8 | 1.84 | **2.85** | 10, 2 | descent 30% |
| `K_1/GF(2^43)` | 32.1 | 2.81 | 2.89 | 23, 2 | collection 34% |
| `K_1/GF(2^47)` | 36.6 | 4.53 | 3.56 | 55, 2 | collection 48% |
| `K_0/GF(2^41)` | 39.0 | 6.38 | 4.51 | 104, 2 | collection 56% |
| `K_0/GF(2^53)` | 44.3 | 10.3 | 6.76 | 263, 2 | collection 62% |
| `K_0/GF(2^61)` | 47.2 | 13.5 | 8.68 | 448, 2 | collection 64% |

- **The law:** `6.38 · (r/2^39)^{1/6} · (n/41)^{−1/2}`. It crosses one
  near `r = 2^23` at `n = 41`, and at `n = 19` it is already below one at
  `2^18`.
- **The model:** frozen constants plus a descent priced per probe. Its
  minimum is `2.85×` near `2^28`, with no crossing. Its local exponent
  over the four largest sizes, after multiplying the ratio by `√n` to
  remove rho's `n^{−1/2}`, is `0.141`, against the law's `1/6`.

## Targets

1. **Correct.** Every target of every run is recovered and verified,
   index calculus and rho alike. No relation fails verification, and
   counts are identical across repetitions.
2. **Controls hold.**
   - Control 1: every count identical to `ic workflow`'s on every
     measurement run.
   - Control 2: the frozen headline's counts reproduced on current `main`,
     or each difference explained.
3. **Exponent.**
   - Fit `ln(ratio · √n)` against `ln r` over the four largest sizes and
     give the slope `β` with its 95% interval.
   - The law's `1/6` is *consistent* if the interval contains it, and
     *falsified at these sizes* if it does not.
   - The model's `0.141` is read the same way.
4. **The small end.**
   - The law's prediction there (`0.83`, `1.20`, `1.18` at the three
     smallest sizes) is *falsified* if each measured ratio's interval
     lies above twice the law's value.
   - The model's minimum is *confirmed* if the least measured ratio falls
     at `K_0/GF(2^37)` or a neighbour (`K_1/GF(2^45)`, `K_1/GF(2^43)`),
     within a factor of two of `2.85`.
5. **Crossing.**
   - A size whose ratio interval lies wholly below one is a crossing of
     the matched reference at `k = 32`, with every phase priced.
   - It is not claimed until two fresh seed sets (`M5`, `M6`, seeds 205
     and 206) repeat it.

## Inadmissible

- Choosing `(c, m)` on a measurement set.
- Tuning the window, unit or cap rules per size.
- Changing `k`.
- Dropping a phase, a failed run or a repetition.
- Multi-threaded times.
- Pricing rho's step at the implemented single-inversion walk.
- Using the Kuhn–Struik formula where the batch measurement exists.
- Reusing a factor base, table or log database across runs or
  repetitions.
- Quoting the first repetition's time as the figure, or discarding it
  unreported.

## Stop and abandon

- **A verification failure:** stop that size, diagnose, and report it
  whatever the outcome. No figure is given for the size until it is
  resolved.
- **Control 1 fails:** the pricer is wrong. Fix it before any figure is
  quoted.
- **The repetition spread** (max over min of the index-calculus total) is
  above `1.25` at a size: rerun that size once with double the
  repetitions, and report both.
- **The table budget refuses a grid point:** record it and continue.

## Host and noise (AGENTS.md §10)

**Recorded with the results:** commit and tree state, `rustc --version`,
CPU model and flags (`popcnt`, `avx2`, `avx512*`, `pclmulqdq`), logical
cores, memory, OS, and the `uptime` load before and after each size.

**Controls:**

- Every timed process runs with `RAYON_NUM_THREADS=1` under
  `taskset -c 2`.
- Nothing else runs while it does.
- The per-set repetitions are the A/A spread, and the unit is interleaved
  with them.

**Hardware class.** One x86-64 cloud container. It says nothing about
Arm64, GPUs or other x86-64 hosts.

**Suite.** The AGENTS.md §8 frozen WDSat suite does not apply: no code
path of the method changes, and the round adds a measurement tool.

## Class

**Accounting.** No index-calculus algorithm changes. The thread is
re-priced on current `main`, its descent is priced at the probes it
scans, and the measurement is extended to new sizes.

## Amendment 1, before any declared run

Committed with the pricer, before any declared size or seed set ran. It
comes from two smoke tests of the pricer on curves that are not among
the declared sizes:

- `K_1/GF(2^17)`, seed 999, 8 columns, `m = 2`;
- `K_0/GF(2^39)`, seed 998, 24 columns, `m = 3`.

Both passed Control 1 against `ic workflow`, run on one thread and on
four. Four things came out of them.

1. **The base grows eight orbits at a time.**
   - `build_subgroup_orbit_factor_base_with_cost` adds eight
     representatives a round and stops at the first base with at least
     the requested points. So a request for `c` columns builds
     `8⌈c/8⌉` of them, and never fewer than eight. At `n = 17`, a request
     for 68 points (2 columns) built 272 (8).
   - The grid is therefore read in actual columns. Requests that build
     the same base are one grid point, run once.
   - The edge rule, at most three steps as declared:
     - a minimum on the lowest grid point extends the grid down by one
       batch (eight columns) at a time, never below eight;
     - a minimum on the highest extends it up by `×√2`, rounded up to a
       whole batch.
   - At the three smallest sizes the model's optima (2 to 5 columns) are
     below what the thread's recipe can build. That is the recipe's
     property, and it stays in the measurement.
2. **The workflow's constructions are priced, and now on their own
   clocks.**
   - The workflow builds the base's projected orbit map in the general,
     big-integer arithmetic several times a run: for the selection's
     column count, for each collector and coverage it creates (two of
     each when extension units run), for the log solver, and for the
     descent solver.
   - At `K_0/GF(2^39)` with 1,872 points these constructions came to
     674K of 970K units, 69%. The work the model prices (selection,
     build, collection, linear algebra, descent) came to 296K.
   - The pricer now keeps those constructions on clocks of their own:
     `select_projection`, `collect_setup`, `logs_setup` and
     `descent_setup`, beside `setup`. The clocks still sum to the whole
     run, and every phase stays in `S` as declared.
   - The write-up reports three groups at every size, as read-outs of
     the same total: the work, the constructions, and the verification.
     The work alone is a stage diagnostic and is labelled as one.
3. **`setup` is the curve's construction** (`experiment::curve`: the
   group order, its factorisation and a generator). The rho side
   receives this already built. It stays in `S` as declared, and it is
   reported as its own phase so that its weight is visible.
4. **One smoke test's spread exceeded the threshold.** Over 15
   repetitions at `n = 17` the spread was 1.29, above 1.25. The declared
   rule applies (rerun with double the repetitions), unchanged.

No target, size, seed set, rule or prediction changes.
