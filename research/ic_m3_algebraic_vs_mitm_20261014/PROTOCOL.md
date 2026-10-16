# Protocol: algebraic and SAT decomposition for m ≥ 3 against the meet-in-the-middle floor

Registered 2026-10-14. This file is committed before any claim run of the
grid below. The pilot (§8) came first, on a different seed and with fewer
targets. Its records are kept and labelled as pilot data, and they enter no
fit.

## 1. Question

Binary Koblitz curves `K_a / F_{2^n}`, `n` prime, with the standard subspace
factor base `F = {P : x(P) ∈ V}`, where `V = span(1, z, …, z^{d−1})` and both
signs are included.

The repository's pair-table root-finders are exact. They also sit on a
generic meet-in-the-middle floor:
- `m = 3` costs `|F|` table probes per target;
- `m = 4` costs `|F|²` probes per target;
- the `|F|²` table is built once.

The question is whether any algebraic or SAT method that solves the
summation-polynomial system on `V^m` (Weil descent to `F₂`, then Gröbner or
SAT) has a cost per **relation** that grows more slowly in `n` than the
floor's on the same `(n, d)` cells. This applies to `m = 3` and `m = 4`.

**Prediction (declared before data).** No arm grows more slowly than the
floor. Prior evidence for this:
- SAT behaves as an enumerator, at about 1 conflict per candidate triple
  (RESEARCH_SAT_SEMAEV).
- The F4 deciding degree is about `ℓ/2` (RESEARCH_PDP_SAT_SMT_NAGAO_SURVEY).
- In `research/pdp_root_finder_20261012/`, the prime-field F4 slopes are 3.2
  and 6.6 bits per unit of `t`, against the root-finders' 1 to 2.

We expect F4 on the S₃ chain to be the cheapest algebraic arm, with a slope
above the floor's.

## 2. Cells (frozen)

**Curve rule.** For each `n`, use the Koblitz curve with the smaller cofactor
(`koblitz_instance_best`). This gives `a = 0` at n = 13 and 31, and
`a = 1` at n = 11, 17, 19, 23 and 29.

**Dimension rule.** `d = ⌈n/m⌉`. This is the classical `d ≈ n/m` regime, in
which a target has `Θ(1)` expected decompositions. The pilot measured
`p = P(decomposable)` between 0.2 and 0.5.

| m | cells `a:n:d` |
|---|---|
| 3 | `0:13:5`, `1:17:6`, `1:19:7`, `1:23:8`, `1:29:10`, `0:31:11` |
| 4 | `1:11:3`, `0:13:4`, `1:17:5`, `1:19:5`, `1:23:6`, `1:29:8`, `0:31:8` |

For `m = 4` the chain system has `4d + 2n` Boolean unknowns. The library's
Boolean polynomial type is capped at 64 variables (`koblitz_groebner::MAX_VARS`),
so n = 23, 29 and 31 are expected to be **unsupported** for the F4, F6-IC and
SAT arms. That is an implementation limit, recorded as such; it is not a
censored cost. `n = 11` is added at `m = 4` so that four cells remain fittable.

## 3. Targets (frozen)

- **Pool.** Per cell, seeded uniform scalars `k ∈ [1, r)` give targets
  `R = [k]G`. Seed `20261014`, mixed with `(m, n, d)`. The pool is 2,048
  targets at `m = 3` and 1,024 at `m = 4`.
- **Floor run.** The floor decides every pool target and counts every
  decomposition (§4).
- **Arm targets.** In pool order, take the first 6 decomposable and the first
  6 non-decomposable targets. They are run interleaved (yes, no, yes, …).
  Every algebraic arm sees the same 12 targets.

## 4. Arms (frozen)

Instrument: `examples/ic_algebraic_floor/` (native Rust, one binary).

| arm | system and solver | unknowns |
|---|---|---|
| `floor-mitm` | `PairSumTable` sorted scan, stopping at the first witness. `m = 3`: `R − P_k` for each `k`, one probe each; `m = 4`: `R − (P_k + P_l)`. Exact and complete, so it is also the ground truth. | — |
| `f4-chain` | `groebner_decompose`, `SolverEngine::InheritedF4 { max_degree: 3 }` (the repository's default oracle), S₃ chain | `md + (m−2)n` |
| `f6-ic` | `groebner_decompose_f6_ic` (inherited F4 plus factor-base support pruning and group closure, #1333), same system | same |
| `sat` | `sat_decompose_with`: in-house CDCL with native XOR, chain system, `restrict_to_factor_base`, every non-lifting model blocked, so a final UNSAT refutes | same |
| `s4-hybrid-f10`, `s4-hybrid-f12` (`m = 3`) | unsplit `S₄(x₁,x₂,x₃,x_R)`, Weil descent (`pdp_split`), F4 over `F₂[v]/(v²−v)`. Guesses `g = 3d − K` summand unknowns and keeps `K = 10` or `12` free; branches run in order, stopping at the first verified decomposition | `3d` |
| `s22-hybrid-g8` (`m = 4`) | `S₃(x₁,x₂,u)` and `S₄(u,x₃,x₄,x_R)` (`pdp_split` plan `[2,2]`), `g = 8` | `4d + n` |

Every arm answers the same decision: does `R` decompose, and as which
summands. Each stops at its first decomposition. A witness is re-added in the
group before it counts.

## 5. Measurement, estimator, statistic and decision rule (frozen)

**Unit.** Wall-clock nanoseconds on one core. Every algebraic solve runs in
its own child process with `RAYON_NUM_THREADS=1`. The floor runs in the
orchestrator.

**Secondary counted units, per arm** (exponents only, never compared across
arms):
- floor: group additions;
- F4 and F6-IC: 64-bit word XORs;
- hybrid: F4 word XORs;
- SAT: conflicts.

**Isolation.** L0: a shared cloud container with 4 vCPU. The `m = 3` and
`m = 4` grids run concurrently. Wall time is descriptive at this level, and
slopes are read on a log scale spanning many octaves.

**Budgets.**
- Per target: 120 s wall, and an address-space cap of 4 GiB.
- A solve that hits a cap, crashes, or ends without a verdict (budget spent,
  or unresolved branches) is **censored**. Its cost is recorded as at least
  `max(wall, 120 s)`.
- **Early stop.** If the first 4 targets of an (arm, cell) are all censored,
  the cell stops. That arm's larger cells are then not run, and are recorded
  as skipped.

**Cost per relation of an arm on a cell.** One relation is one decomposed
target. The cost is

`C = (p·c_yes + (1 − p)·c_no) / p`

where:
- `p` is the decomposable fraction of the cell's pool, from the exact floor;
- `c_yes` and `c_no` are the arm's mean wall over its decomposable and
  non-decomposable targets.

A cell is a **measurement** only when all 12 targets finished. A cell with
any censored target, or one stopped early, is a **lower bound**.

**Floor cost per relation.**
- Primary (amortized): `mean decision wall / p + table build wall / |F|cols`.
  This spreads the `|F|²` table over the `|F|cols` relations a solve needs,
  where `|F|cols` is points up to sign, i.e. matrix columns.
- Secondary (online): the same without the build term.

**Slope.** OLS of `log₂ C` on `n` across cells, per arm.

**Intervals.** 2,000 bootstrap replicates, seed `20261014`.
- Each cell resamples its pool, which changes `p` and the floor.
- Each arm resamples within its two strata.
- The floor slope is refitted on the same cells in the same replicate.
- Δ = arm slope − floor slope, read as a 95% percentile interval.

**Decision rule, per arm and per `m`.**

1. With at least 3 measured cells:
   - Δ interval entirely below 0: **grows more slowly than the floor**. This
     is a positive result and triggers §7.
   - Δ interval entirely above 0: **grows faster than the floor**.
   - Otherwise: **indistinguishable on this range**.
2. With fewer than 3 measured cells, lower-bound cells are substituted (the
   lower-bound fit). This is valid only when every lower-bound cell lies
   above the mean `n` of the fitted cells, which makes the fitted slope a
   lower bound on the true slope.
   - Δ interval of that fit entirely above 0: **grows faster than the floor
     (lower-bound fit)**.
   - Otherwise: **undetermined**.
3. With fewer than 3 usable cells: **undetermined**.

**Exactness gate.**
- An (arm, cell) whose verdict disagrees with the floor on any target, or
  that returns a witness that does not add up, is **inexact**. It is excluded
  from fits and reported.
- The floor itself is checked against exhaustive enumeration
  (`pdp_split::exhaustive_decompositions`) on the 12 arm targets, wherever
  that enumeration has at most 5·10⁶ multisets.

**Reported per cell.**
- `C` with its interval, and `C ÷ floor`.
- Status, and censoring counts.
- The cheapest arm.

## 6. What does not count

- A cell that is a lower bound is never presented as a cost.
- A per-attempt cost is never presented as a per-relation cost.
- The floor's online cost is never presented as its full cost.
- No end-to-end, `vs_rho` or exponent claim follows from this stage
  diagnostic (AGENTS.md §8), except through §7.

## 7. Conditional follow-up

If any arm's verdict is **grows more slowly than the floor**, run an ecbench
complete-DLP session on the largest such cell's curve:
- that arm as the decomposition oracle;
- against `rho.signed_frobenius` and `rho.signed_frobenius_strong`;
- under `ecbench-measure`.

Otherwise the follow-up is not run, and this file says so.

## 8. Pilot (before freezing; not claim data)

**Setup.** Seed 1, 256-target pool, 2 decomposable and 2 non-decomposable
targets, budget 60 s, multithreaded solver. The `m = 4` pilot was stopped after
SAT at n = 17 (2 of 4 targets, both censored); its split-hybrid rows at
n = 13 had already censored on all targets.

**Cells and arms.**
- `m = 3`, n = 13 and 17: `f4-chain`, `f6-ic`, `sat`, and `s4-hybrid` with
  `g` = 0, 4 and 8.
- `m = 3`, n = 19: `s4-hybrid` with `g` = 9, 11 and 13.
- `m = 4`, n = 13 and 17: `f4-chain`, `f6-ic`, `sat`, and `s22-hybrid` with
  `g` = 4 and 8.

**Used to fix:**
- **d rule:** `p` at `⌈n/3⌉` was 0.31–0.48.
- **Hybrid guess rule:** the best `g` was 8 at n = 17 (10 free) and 9 at
  n = 19 (12 free), so the frozen arms keep 10 and 12 free.
- **Budget:** F4 and F6-IC decided in 0.004–0.3 s; SAT censored from n = 17
  (`m = 3`) and n = 13 (`m = 4`); the split hybrid censored at n = 13
  (`m = 4`).
- **64-variable limit at `m = 4`.**

Records: `pilot/`.
