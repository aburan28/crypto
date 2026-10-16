# Nagao's Riemann–Roch system against the symmetrised Semaev system on `E / F_{p^n}`, odd characteristic

Registered before any registered cell runs. §6 lists everything that ran before
registration and what it showed.

## 1. Question

Every Nagao-vs-Semaev comparison in this repository is binary at prime extension degree
(`research/nagao_relations/`, the RR solver panel, `ic_rr_norm_ladder_20260930`,
`ic_gb_ladder_20261003`, `ic_dense_ladder_20261004`). The PDP survey
(`research/notes/index-calculus/RESEARCH_PDP_SAT_SMT_NAGAO_SURVEY.md` §5.2) says the gain
Nagao offers "needs a subfield". The survey's open list names "composite-degree or
odd-characteristic Nagao vs Semaev under one engine" as having no data. This study
supplies that data.

**Setting.** `E : y² = x³ + ax + b` over `F_q`, `q = p^n`, `p = 65521`, with `a, b`
uniform in `F_q` and not both in `F_p`. Factor base `F = {P : x(P) ∈ F_p}` (Gaudry,
Diem). The task: decide whether a target `R` is `Σ_{i≤m} P_i` with `P_i ∈ F`, and list
every such decomposition.

**Question.** Under one Gröbner engine and one zero-dimensional solver, which system
costs less per target, and by how much, at each `(n, m)`?

- `sem`: Semaev's `S_{m+1}` in the elementary symmetric functions `e_1..e_m ∈ F_p`,
  specialised at `x_R` and Weil-descended. `m` unknowns, `n` equations, total degree
  `2^{m−1}`.
- `rr`: Nagao's form. `f ∈ L((m+1)O)` monic, `f(−R) = 0`, `m − 1` coefficients
  Weil-restricted. `n(m−1)` unknowns, `m(n−1)` quadratic equations requiring the norm
  cofactor `g = N_f/(x − x_R)` to lie in `F_p[x]`.

Both systems have the same Bézout number when `m = n`: `2^{m(m−1)}`.

## 2. Cells, targets and seeds

| cell `(n, m)` | `sem` (vars, eqs, degree) | `rr` (vars, eqs, degree) | regime |
|:--|:--|:--|:--|
| (2, 2) | 2, 2, 2 | 2, 2, 2 | square |
| (3, 2) | 2, 3, 2 | 3, 4, 2 | over-determined (`m = n − 1`) |
| (3, 3) | 3, 3, 4 | 6, 6, 2 | square (Gaudry/Diem `m = n`) |
| (4, 3) | 3, 4, 4 | 8, 9, 2 | over-determined (Joux–Vitse `m = n − 1`) |
| (4, 4) | 4, 4, 8 | 12, 12, 2 | square |
| (5, 4) | 4, 5, 8 | 15, 16, 2 | over-determined |

- **Seed** 20261016 for every cell. The curve and modulus come from the seed and the
  cell. Targets come from their own stream, so the arm list never moves them.
- **Targets.** For `m ≤ 3`: `T = 8` planted (`R = Σ P_i` with `P_i` drawn uniformly
  from `F`) and `T = 8` uniform (`R` a uniform point). For `m = 4`: `T = 4` of each,
  because one `sem` solve at (4,4) takes about six minutes and every `rr` solve at
  `m = 4` ran out of memory in pre-flight. Both kinds redraw while `x_R ∈ F_p` or
  `y_R = 0`.
- **Why `T = 8` is enough.** The cost metric is a deterministic count of `F_p`
  multiplications. In pre-flight it varied by at most a few hundred multiplications
  between targets of the same cell (§6). The question is a ratio of one to two orders
  of magnitude, not a small effect.

## 3. Instrument (identical for both arms)

- **Harness:** `examples/nagao_semaev_fpn.rs`, at this file's commit.
- **Solver:** `crypto_lib::cryptanalysis::jv_cover::solve_zero_dim` with
  `F4Options::new(Grevlex, 60).with_budget(B)`, built fresh for each solve so the
  deadline is per solve. Inside: `f4_fp::f4`, then the multiplication matrices on the
  standard monomials, a Krylov minimal polynomial, its roots in `F_p` and the
  eigenvectors. The solver refuses an ideal with more than 4096 standard monomials.
- **Budget:** `B = 1800 s` per solve (F4 deadline), with a wall cap of `B + 600 s` per
  process.
- **Process and memory.** One harness process per (cell, kind, target, arm)
  (`--only`), under `prlimit --as=6 GB` on a 7.7 GB, 2-core host, run one at a time by
  [`run.sh`](run.sh). A process that dies leaves a `failure` row with its exit code and
  reason (`oom`, `wall_timeout`, `crash`). Such a row is an unknown cost, never a zero.
- **One verifier** turns every solution into `g = ∏(x − x_i) ∈ F_p[x]` and classifies
  it as:
  - no split over `F_p`;
  - twist (some `h(x_i)` is a non-square in `F_q`);
  - no sign pattern;
  - a relation (`Σ ±P_i = R`, checked by curve arithmetic).
- **Not timed:** the per-curve `S_{m+1}` symmetrisation is setup, reported as its own
  phase (`sem_build_ms`). Per-target system construction is reported as `build_ms`.

## 4. Metrics

1. **Primary:** `f4_muls`, the `F_p` multiplications spent in F4 row reduction (the
   engine's own counter).
2. **Secondary, same unit:** `lin_muls`, the multiplications of the post-F4
   linear-algebra stage. Their sum is the per-target solve cost.
3. **Structural:**
   - `solving_degree`;
   - `max_cols`, the widest Macaulay/F4 matrix;
   - `delta`, the degree of the ideal (number of standard monomials);
   - algebraic solutions, broken down by verifier outcome.
4. **Wall time:** `f4_ms` and `lin_ms`. Recorded, not isolated (shared two-core host),
   and not used for any claim.

## 5. Decision rule and predictions

- **Correctness gate.** On every target where both arms complete, the verified relation
  sets must be equal. Every planted decomposition must be found by every completing arm.
  A violation invalidates the cell; it is never averaged away.
- **Cost ratio.** Per cell, `ρ = median(f4_muls + lin_muls)_rr / median(…)_sem`, with
  the min–max range. A cell where either arm times out on any target reports that arm's
  cost as `> B`. Cost then appears as a lower bound or "unknown", never as a number.
- **Prediction (from the pre-flight in §6).**
  - `ρ ≈ 1` at (2,2).
  - `ρ > 1` at every cell with `n ≥ 3`: the symmetrised Semaev system is cheaper under
    this engine. Pre-flight read about 2× at (3,2), about 2.9× at (3,3) on F4 plus
    extraction (16× on F4 alone) and about 255× at (4,3).
  - At `m = 4`, `rr` fails by memory on every target while `sem` completes.
- **Falsifiers.**
  - `ρ ≤ 1` at any cell with `n ≥ 3` where both arms complete all targets.
  - Any `m = 4` target on which `rr` completes within the cap while `sem` fails.
- **Not established by this study, whatever the outcome:**
  - anything about an F5 or msolve engine;
  - asymptotics in `n`;
  - anything about prime-field curves or about end-to-end cost against rho.

  This is a stage diagnostic of the decomposition step (AGENTS.md "Research work
  belongs in pull requests"). End-to-end cost and speedup stay unset.

## 6. What ran before registration (not data)

- **Correctness gate.** Output in [`checks/`](checks/).
  - `validate_p37_m2_m3.txt` (`--validate --p 37 --targets 6`, cells (2,2), (3,2),
    (3,3), (4,3), both arms). On every target, both arms' verified relation sets equal
    brute force over all `m`-multisets of the factor base, and every planted
    decomposition is found. Passed.
  - `validate_p23_n5m4_sem.txt` (`--validate --p 23 --targets 6 --cells 5,4
    --arms sem`). The same check for the `m = 4` Semaev path. Passed.
  - **Not obtained.** The `m = 4` `rr` arm cannot be checked by solving on this host,
    because it runs out of memory (below). A first `--validate --p 23` at (4,4) also
    returned `undetermined` for both arms within its 600 s budget.
  - `rr_plant_check_p65521.txt` (`--rr-plant-check`) is the `rr` construction's check
    for every cell, `m = 4` included. On 8 planted targets per cell it interpolates the
    true `f ∈ L((m+1)O)` through the planted points and checks two things: its Weil
    coordinates satisfy every `rr` equation; and moving one coordinate by 1 breaks
    some equation (negative control). Passed 48/48.
  - Every validation run also checks the symmetrisation identity
    `T(e(x), X) = S_{m+1}(x, X)` at four random points of `F_q` per cell.
- **Pre-flight timing.** Seed 1, `p = 65521`, two planted targets per cell, before this
  file existed:

  | cell | arm | vars | eqs | solving degree | max cols | `f4_muls` | `lin_muls` | `delta` | outcome |
  |:--|:--|--:|--:|--:|--:|--:|--:|--:|:--|
  | (2,2) | sem | 2 | 2 | 3 | 10 | 129 | 1 066–1 978 | 4 | planted found |
  | (2,2) | rr | 2 | 2 | 3 | 10 | 84 | 1 059–1 324 | 4 | planted found |
  | (3,2) | sem | 2 | 3 | 3 | 8 | 155 | 48 | 1 | planted found |
  | (3,2) | rr | 3 | 4 | 3 | 15 | 339 | 52 | 1 | planted found |
  | (3,3) | sem | 3 | 3 | 10 | 108 | 62 728 | 520 316 | 64 | planted found |
  | (3,3) | rr | 6 | 6 | 7 | 390 | 1 028 000 | 653 453 | 64 | planted found |
  | (4,3) | sem | 3 | 4 | 3 | 81 | 127 719 | 52 | 1 | planted found |
  | (4,3) | rr | 8 | 9 | 3 | 1 514 | 32 586 718 | 72 | 1 | planted found |
  | (4,4) | sem | 4 | 4 | 29 | 5 507 | 5.88·10⁹ | 1.07·10¹¹ | 4096 | planted found (F4 2.6 s, extraction 354 s) |
  | (4,4) | rr | 12 | 12 | — | — | — | — | — | killed at 4.2 GB resident after 3 min (out of memory; no row) |
  | (5,4) | sem | 4 | 5 | 3 | 3 679 | 1.11·10¹⁰ | 56 | 1 | planted found (7 s, 0.32 GB) |
  | (5,4) | rr | 15 | 16 | — | — | — | — | — | killed (exit 137) after 5 min, 5.5 GB resident |

  The `rr` arm at `m = 4` first overflowed the main thread's 8 MB stack inside F4. The
  harness now runs on a thread with a 1 GB stack, and the out-of-memory kills above
  were measured after that fix.

  The first pre-flight used `f4_fp::solve` and was discarded. That solver finds roots
  by trial over all of `F_p`, and when a grevlex basis has no univariate element it
  enumerates every value of the first variable and reruns F4 each time. It also built
  its deadline once per process, so every target after the first had no budget. Both
  defects are instrument defects, not properties of either system. The registered
  instrument is `solve_zero_dim` with a per-solve deadline.

## 7. Amendment 1 (2026-10-16): external engines, an F5-family solver included

Registered before any external-engine data. The F4 run of §1–§6 was in progress and
nearly complete when this amendment was written. Its finished cells were seen; its
outcome is not changed by this amendment. The owner asked for an F5 run. The F4
result is specific to one engine (`f4_fp`, F4 without signatures), and §5 says so.
This amendment asks whether the ranking survives a stronger engine.

- **Systems.** The registered systems themselves. The harness at this amendment's
  commit, with `--dump-dir … --dump-only` and the same seed, cells and per-kind target
  counts as §2, writes every `sem` and `rr` system. The target streams are those of
  the registered run. `systems.sha256` hashes every file.
- **Engines.** One process per system, one at a time, single-threaded, under the §3
  caps: `prlimit --as=6 GB`, wall `3600 s`.
  - **Singular 4.4.1 `sba(i)`**, the signature-based Gröbner basis algorithm (F5
    family), under `dp` (grevlex). Defaults: `sbaOrder` and the rewrite criterion
    are Singular's own.
  - **msolve 0.11.0** (commit `a18d0de`, built from source against the Sage
    environment's FLINT/GMP/MPFR), `-g 2 -t 1`, the reduced grevlex basis. msolve is
    F4, not F5. It is included as the reference state-of-the-art grevlex engine, so
    that `sba`'s reading is not compared with a weak F4 alone.
- **Scope.** Gröbner basis only; no solution extraction. Extraction is the
  `lin_muls` stage of §4, which depends on `delta` and not on the arm. Singular's
  `vdim` must equal the §3 solver's `delta` wherever both exist (a consistency check).
- **Metrics.**
  - Wall seconds and peak RSS (`/usr/bin/time`).
  - Each engine's own time: msolve's "overall time"; for `sba`, Singular's `rtimer`
    around the `sba` call, in ms.
  - msolve's zero reductions, rows reduced and widest matrix.
  - Exit code: 124 = wall timeout; out of memory under the cap is an allocation
    failure.

  Wall time is recorded on the shared, non-isolated two-core host. It is used only as
  the within-engine ratio `rr/sem` on the same system pair, run back to back, never
  as an absolute speed. Times below 0.05 s are reported as "below resolution" and
  carry no ratio.
- **Decision rule.** Per engine and cell: `ρ_engine = median(engine_time)_rr /
  median(engine_time)_sem` with min–max ranges. A failed arm is "> cap" (memory) or
  "> 3600 s".
- **Predictions.**
  - Under `sba` the `rr` disadvantage shrinks relative to `f4_fp`, because F5 removes
    the reductions to zero that a square, many-variable quadratic system produces.
  - `ρ_sba > 1` remains at (4,3).
  - `rr` at `m = 4` stays out of reach for at least one engine.
- **Falsifiers.**
  - `ρ_engine < 1` (rr cheaper) at any cell with `n ≥ 3` where both arms complete
    every target. This would reverse the §5 ranking for that engine.
  - `rr` completing at `m = 4` within the caps under either engine. This would turn
    the F4 run's out-of-memory rows from "beyond reach" into "beyond `f4_fp`'s
    reach".
