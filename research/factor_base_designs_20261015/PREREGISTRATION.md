# Factor-base designs: preregistration (frozen before any data)

**Frozen:** 2026-10-15, before any run in this directory. Later changes go in
an appended, dated "Amendments" section. They never edit the text above it.
**Study:** `research/factor_base_designs_20261015/`. The design register and
the Lean proof map are in [`README.md`](README.md).
**Kind of work:** solver-axis stage diagnostics (AGENTS.md §2). Nothing here
prices a pipeline or reports `S`. No run here can show a speedup.

## 1. The quantity every experiment measures

Every arm is scored by the **slope of the solving degree** against the
factor-base entropy per summand. Write `b = log₂|V|` per summand. That is
`ℓ` for a binary subspace and `t` for a depth-`t` tower or a `t`-bit exponent
interval. The score is `σ = ΔD / Δb`, fitted by least squares over the
completed cells of one (arm, `m`) series.

Why the slope and not `D` at a point: under the repository's yield
conservation (`AR:KN-FIND-007`), yield per trial depends only on `|V|` and
`m`. One more bit per summand buys `m` bits of yield. At deployed sizes,
one more unit of `D` costs several bits of Gröbner work (README,
"Break-even", with integer instances checked in Lean, `BreakEven.lean`).

**Decision thresholds (fixed now):**

| reading | meaning |
|:--|:--|
| `σ ≤ 0.25` on at least 3 completed cells | **candidate**: degree grows slowly enough for an interior optimum. Replicate on a second prime before any further claim |
| `0.25 < σ < 0.5` | **weak**: reported. Follow up only if it beats the incumbent by ≥ 1 in `D` at equal `b` |
| `σ ≥ 0.5` | **no middle ground at this arity**. Kept and reported |
| fewer than 3 completed cells | **inconclusive** |

The **incumbent** for prime fields is the Kummer tower (PKM), measured with
`examples/pkm_tower_pilot.rs` (`research/pkm_tower_pilot_20260924/`).

## 2. Instrument

The instrument is the repository's degree-bounded F4 (`src/cryptanalysis/f4_fp.rs`)
in grevlex, with `solving_degree_max` as `D`, exactly as the PKM pilot defines
it. Every system is checked to vanish at its planted solution before it is
measured. A failing check aborts the run and is never recorded as a degree.
Limits: 300 s per system, `ulimit -v 13000000`. A system that hits a limit is
recorded with its reason and counts as a missing cell.

## 3. Experiments in this PR (prime field, `p₁ = 786433 = 3·2^18 + 1`)

`p₁` is the pilot's prime, so arm A0 must reproduce committed pilot rows,
which serves as the instrument check. Seeds are derived from `0x46424453`
("FBDS") combined with (arm, `m`, `t`, target index). Each cell has 2 planted
targets: the target is the sum of `m` random base points, so a solution
exists. `m ∈ {2, 3}`. `t` rises from 2 until a limit is hit.

| arm | design (README #) | `V` per summand | membership presentation |
|:--|:--|:--|:--|
| **A0** Kummer tower (incumbent) | — | `g·μ_{2^t}` | `y_{j+1} = y_j²`, `y_{t−1}² = 1` (`t` unknowns) |
| **A1z** exponent interval, ζ | #1 | the same set `g·μ_{2^t}` | `x = g·∏_j(1 + (ζ^{2^j} − 1)b_j)`, `b_j² = b_j`, partial products (`2t − 1` unknowns) |
| **A1g** exponent interval, generic | #1 | `{h·g^k : 0 ≤ k < 2^t}`, `ord(g)` not a power of 2 | as A1z |
| **A4** parallel tower | #4 | `{u + v}`, `u ∈ h·μ_{2^a}`, `v ∈ h'·μ_{2^b}`, `a = ⌈t/2⌉`, `b = ⌊t/2⌋` | two Kummer towers and `x = u + v` (`t + 1` unknowns) |
| **A11** norm chain | #11 | `g·μ_{2^t}` | Riemann–Roch function through `−R`, `G(gX) ∣ X^{2^t} − 1` as `t` squaring steps mod `G` |
| **N** null | control | as A0 | as A0, with `S₃`'s support given random coefficients and a planted zero |

A1z isolates the **presentation** from the **set**: it has the same `V` as
A0. A1g tests whether dropping the smoothness requirement costs extra
degree. A4's `|V|` is recorded exactly (`u + v` can collide). A cell with
`|V| < 2^t` reports its true `b`.

**Predictions recorded now:**
- **A1z**: `D` within ±1 of A0 at equal `b`, and `σ` no smaller. The extra
  variables make a win unlikely.
- **A1g**: equal to A1z within ±1. If so, the smoothness restriction on
  towers can be removed at no extra degree.
- **A4**: if `D` tracks elimination depth, A4 at `t` matches A0 at `⌈t/2⌉`
  plus at most 1. If `D` tracks the Bezout number, A4 matches A0 at equal `t`.
- **A11**: no prediction. This is the arm with the most room to differ.

## 4. Designed, not run in this PR

Each one gets its own amendment before running. The arms are in the README
register.

- E2 Lucas-ladder torus interval; E3 elliptic interval on an auxiliary curve;
  E5 radix-3 auxiliary towers; E6 torus tower at analogue primes
  `p ≡ −1 (mod 2^k)` (the P-521 shape); E7 fan-in roots; E8 lopsided depths;
  E9 weighted grevlex on towers; E10 split chain + towers, with treewidth
  recorded.
- E12 / E13: binary norm chain and algebraic large prime, through
  `research/nagao_relations/` (its open item 1).
- E14 / E15: learned frames and family-slope selection, in the
  ml-cryptanalysis CSD pipeline with a held-out field family.
- E16: arity portfolio (engineering).
- **E-K1**: Gauss-period window in `s = x + 1/x` on Koblitz analogues
  `n ∈ {11, 23}` (`2n + 1` prime, type-II normal basis), against the
  `x`-window and the `u`-frame bases of `koblitz_symmetrised.rs`.
- **E-K4** (MinRank modelling) and **E-K6** (shared-trunk towers): reopened,
  because no Lean-checkable statement disposes of them.

## 5. Stop rules

- If A0 fails to reproduce the pilot's committed `D` on the same cells, the
  study stops (instrument fault).
- An arm reaching **candidate** on `p₁` is replicated on `p₂ = 2013265921`
  before it is described as anything more than a single-prime reading.

## Amendments

### Amendment 1, 2026-10-15: implementation details, before the first data row

Registered by the experiment lane (`examples/fb_designs_slope.rs`) before any
row exists in `runs/`. Each item settles a detail §2–§3 leave open. None
changes an arm, a threshold or a prediction.

1. **Instrument version.** `f4_fp` has applied the Gebauer–Möller pair
   criteria by default since commit `3d4634596` (2026-10-11), after the pilot.
   This study runs the current default. Reason: before any data, the 16 planted
   systems of the pilot's Kummer cells (`m = 2`, `t ≤ 6`; `m = 3`, `t ≤ 3`)
   were rebuilt with the pilot's seed rule. Their `D`, `x_R`, curve, tower,
   unknown and equation counts, `inconsistent` and `planted_ok` equal every
   committed pilot row. Only `degree_reached` (the certification tail, not `D`)
   is lower. With `F4_FP_GM=0` it matches too. `degree_reached` is recorded,
   but it is not comparable with the pilot's.
2. **Instrument check (§5).** "The same cells" means A0 rebuilt with the
   pilot's seed rule and draw order (`repro`), compared row by row. This
   study's own seeds give other towers and curves, so the A0 series is also
   compared with the pilot cell by cell at matching (`m`, `t`). That
   comparison is reported but does not stop the study.
3. **Seeds.** The cell seed is `0x46424453 ⊕ arm≪40 ⊕ m≪32 ⊕ t≪16`, with arm
   index A0 = 0, N = 1, A1z = 2, A1g = 3, A4 = 4 and A11 = 5. It draws the set
   and the target curve, which both targets of a cell share. Target `k` is
   planted from `cell seed ⊕ (k+1)≪8`. N's random coefficients come from the
   same stream.
4. **Completed cell and fit.** A cell is completed when both of its planted
   systems finish within the budget and the planted solution satisfies F4's
   basis. The cell's `D` is the mean of the two. `σ` is the least-squares slope
   of the cell `D` on `b`. Any verdict other than "inconclusive" needs at
   least 3 completed cells.
5. **Limits and the stop.** `t` rises from 2 to at most 18, because
   `2^18 ∥ p₁ − 1`. Each system runs in its own child process under
   `ulimit -v 13000000`, with F4's 300 s budget. The child is killed at 420 s
   if F4 has not returned to its deadline check. A series stops after the
   first cell in which any system hits a limit. That cell is recorded with its
   reason and counts as missing.
6. **N at `m = 3`.** The pilot's null is defined only at `m = 2`. Here every
   `S₃` of the chain gets random coefficients on its own support, shifted to
   vanish at the planted solution.
7. **A1.** `x = π_{t−1}(1 + β_{t−1}b_{t−1})` is substituted into `S₃`; it is
   not an unknown. That gives `2t − 1` unknowns per summand. For A1g, `h` and
   `g` are drawn uniformly. `g` is redrawn when its order is a power of 2 or
   when the `2^t` elements repeat. `ord(g)` is recorded.
8. **A4.** `x` is an unknown, with the linear equation `x − h·y_0 − h'·z_0`,
   giving `t + 1` unknowns as in §3.
9. **A11.** The function passes through `−R`:
   `f = y − P(X)` with `P(X) = y(−R) + Σ_{i≥1} c_i (X^i − x_R^i)`. At `m = 2`
   this is the line `y = λ(X − x_R) + y(−R)` with `λ = c_1`. `G` comes from
   exact division of `N(X) = X³ + AX + B − P(X)²` by `X − x_R`, and the zero
   remainder is checked symbolically. The coefficient equations are
   `lead(G)·g^{m−i}·s_i − G_i = 0`. At `m = 2`, `lead = 1`, so they read
   `g·s_1 = e_1` and `g²·s_0 = e_0`. **`m = 3` is run:** `deg P = 2`, `G` is
   a cubic with leading coefficient `−c_2²`, and each `q_j` is linear. That
   gives `2 + 3 + 5t` unknowns and `3 + 5t + 3` equations. `m = 2` has
   `3 + 3t` unknowns and `3t + 4` equations. `b = t`. The presentation is
   checked in two ways. Every element of `V` is `g` times a `2^t`-th root of
   unity. On 128 random `m`-subsets, half of them with one element outside
   `V`, `X^{2^t} ≡ 1 mod ∏(X − v/g)` holds exactly when every element is in
   `V`.

### Amendment 2, 2026-10-15: A11 follow-up (F1–F4), before any of its rows

Registered after the p₁ run (`reports/run-20261015T1657Z`, review in
`reports/review-20261015T2030Z`) and before any row of F1–F4 exists. It adds
experiments. It changes no p₁ row, verdict, threshold or prediction. The
instrument, the per-system limits (300 s F4 budget, 420 s hard kill,
`ulimit -v 13000000`), the stop rule (Amendment 1 item 5), the completed-cell
rule and the σ fit are those of Amendment 1. The exceptions are stated below.

**Why.** On p₁, A11 at `m = 3` was the one arm below A0 on both D (4 vs 6) and
matrix width (≈ 4× narrower). That rests on 3 completed cells, one prime and
no null of its own shape (review, point 4). §5 requires replication. The
review lists F1–F4 as owed.

**Seeds.** As in Amendment 1 item 3, with a prime index:
`cell seed = 0x46424453 ⊕ prime≪48 ⊕ arm≪40 ⊕ m≪32 ⊕ t≪16`, where p₁ = 0
(so every p₁ seed is unchanged) and p₂ = 1. New arm indices: N11 = 6.

- **F1, replication on p₂ = 2013265921 = 15·2^27 + 1.** Arms A0 and A11, at
  `m = 2` and `m = 3`, with `t` rising from 2 (at most 27) until the stop
  rule. Before any p₂ row, check that `f4_fp` and the field arithmetic are
  exact at p₂ (`p₂ < 2^31`): rebuild one A0 and one A11 planted system at
  `t = 2`, check the planted solution, and record the check. If the
  instrument cannot run at p₂, F1 stops and the reason is recorded.
- **F2, A11-shape null (N11), p₁, `m = 2` and `m = 3`.** N11 is A11 with the
  curve removed, exactly as N is A0 with the curve removed. The equations
  that define `s_i` from the function's coefficients (Amendment 1 item 9,
  `lead(G)·g^{m−i}·s_i − G_i = 0`) are replaced by polynomials on the same
  monomial support, with independent uniform random coefficients, shifted
  by a constant so they vanish at the planted solution. The chain equations
  and `r_t = 1` are kept unchanged. The planted solution satisfies every
  equation, and this is checked. Unknown and equation counts equal A11's.
- **F3, the `m = 3`, `t = 5` cell, p₁.** Arms A0 and A11, 2 planted targets
  each, with F4 budget 3600 s and hard kill 3900 s (the only budget change).
  Memory limit unchanged. A system that still hits a limit is recorded as such.
- **F4, `m = 4`, p₁.** Arms A0 and A11, with `t` rising from 2 until the
  stop rule.
  - A0 at `m = 4` uses the `S₃` chain with two intermediate unknowns, as the
    pilot does at `m = 3` with one.
  - A11 at `m = 4` uses `f = Q(X)·y − P(X)` with `Q = X + q₀` and
    `deg P = 2`, passing through `−R` (3 free coefficients, `= m − 1`).
    `N(X) = Q²(X³ + AX + B) − P²` has degree 5. Exact division by `X − x_R`
    gives a quartic `G` with leading coefficient `1 − p₂²`, where `p₂` is the
    leading coefficient of `P`. Rows where that coefficient vanishes are
    recorded and redrawn. The chain uses `r_j` of degree 3 and `q_j` of
    degree 2.
  - Unknown and equation counts are recorded. The excess must be 1
    (`norm_chain_excess`), and this is checked.

**Decision rules, fixed now.**
- *A11 replicates* if, on p₂ at `m = 3`, A11's cell D ≤ A0's cell D − 1 and
  A11's `max_cols` < A0's at every shared completed `t`, with at least 3
  shared completed cells. *Fails to replicate* if any shared completed cell
  has A11's D ≥ A0's. Otherwise: *inconclusive*.
- *Shape, not structure* if N11's D equals A11's in every completed cell of
  F2. *Structure beyond shape* if N11's D exceeds A11's by ≥ 1 in at least 2
  completed cells. Otherwise: *mixed*, reported cell by cell.
- *F3*: reported as recorded (D, width, time). It adds a fourth A11 and A0
  cell to the `m = 3` p₁ σ fits only if both of the cell's systems complete.
- *F4*: σ verdicts as in §1. The threshold was derived at `m = 4`, so here it
  is used as designed.
- **Secondary, descriptive (not a verdict):** the per-bit growth of the
  matrix width, `Δ log₂ max_cols / Δb`, by least squares over completed
  cells, per (arm, `m`, prime). It is recorded because D alone misread A1 on p₁.

**Predictions.**
- F1: A11 replicates.
- F2: shape, not structure. The p₁ null N equalled A0, so the same is
  expected here.
- F3: both systems of A11 complete within 3600 s. A0's may not.
- F4: A0's σ ≥ 0.5. No prediction for A11.

### Amendment 3, 2026-10-15: implementation details of F1–F4, before any of their rows

Registered by the experiment lane (`examples/fb_designs_slope.rs`) before any
row of F1–F4 exists. It changes no arm, threshold, decision rule or
prediction. Items 1 and 2 settle what Amendment 2 leaves open; item 3
corrects a statement in it that cannot hold.

1. **N11's random coefficients** come from the target stream after the
   target is planted, as N's do (Amendment 1 item 3). The `m` replaced
   equations are the first `m` of the A11 system, in the same order.
2. **A cell's `max_cols`** in the F1 rule and the secondary slope is the mean
   over the cell's two systems, as a cell's D is (Amendment 1 item 4).
   Amendment 2 does not say how a cell's width is formed.
3. **The leading coefficient of `G` at `m = 4`.** Amendment 2 says `G` has
   leading coefficient `1 − p₂²` (`p₂` the leading coefficient of `P`). With
   `Q = X + q₀` monic and `deg P = 2`,
   `N(X) = Q²(X³ + AX + B) − P² = X⁵ + (2q₀ − p₂²)X⁴ + …`: `Q²X³` reaches
   `X⁵` and `P²` only `X⁴`. So `N` is monic and `G = N/(X − x_R)` is a monic
   quartic. The harness checks symbolically that `lead(G)` is the constant 1
   in every built system. Consequences: the coefficient equations read
   `g^{4−i}·s_i − G_i = 0`, and the "leading coefficient vanishes" redraw
   cannot trigger. The redraw that can trigger is the other degeneracy: when
   no `f` with `Q` monic passes through the four summands (the 4 × 4 linear
   system for `(q₀, P)` is singular), the target is replanted from the same
   stream, as A11 replants on a vanishing lead at `m = 3`. The free
   coefficients are ordered `q₀, c₁, c₂`, with
   `P(X) = Q(x_R)·y(−R) + c₁(X − x_R) + c₂(X² − x_R²)`. The excess
   (`norm_chain_excess`) is asserted to be 1 at every `m`.

### Amendment 4, 2026-10-16: the `m = 4` series to three cells; distinct targets

Registered by the experiment lane before any G1 row exists. It finalises the
draft committed as `AMENDMENT-4-DRAFT.md` (`1b44b48cc`). The draft's three
open questions were settled with the study's owner, as items 1–3 below. It
adds one experiment (G1) and one construction rule (G2). It changes no
committed row, verdict, threshold or decision rule.

**Why.** F4 left both `m = 4` series inconclusive with 2 completed cells each.
Every `t = 4` system stopped at the 300 s budget, A0 at partial D = 7 and
A11 at partial D = 6. §1's threshold was derived at `m = 4`. F3 showed that a
3600 s budget completes cells 300 s cannot, and that 13 GB is the limit that
binds next (A0, `m = 3`, `t = 5`).

**G1. `m = 4` from `t = 4`, p₁, arms A0 and A11.**
- F4 budget 3600 s, hard kill 3900 s (as F3).
- **Memory limit `ulimit -v 15000000`** (item 1; raised from 13 GB). That is
  the physical memory of the host that runs G1 (15 GB, no swap, no cgroup
  limit). A system that outgrows it is killed, by the allocator's abort or
  by the kernel's OOM killer, and recorded as aborted with its reason. A
  larger limit needs a larger host and a new amendment.
- Seeds as in Amendment 2. With G2 the `t = 4` systems are the F4 systems
  that timed out, unless G2 redraws one; the rows record `target_redraws`.
- `t` rises from 4 until the stop rule, per arm. Files are
  `runs/G1-{arm}-m4.jsonl`.
- An arm's `m = 4` series is its F4 rows at `t = 2, 3` plus its G1 rows. A
  cell enters the σ fit only if both of its systems complete. F4's `t = 4`
  timeout rows stay committed, and G1 supersedes them in the fit.

**G2. Distinct targets (construction; harness flag `--distinct-targets`).**
- *Rule:* a planted target whose **system** (every equation, term for term)
  equals an earlier target's system in the same cell is redrawn from its own
  stream (item 2). Targets are built in index order, so the rule gives the
  same result in every process.
- *Why the system and not `x_R`:* N's and N11's systems take random
  coefficients from the target stream, and A11's depends on `y(−R)`. So equal
  `x_R` does not mean equal systems for every arm, but equal systems is what
  makes a cell hold one system twice.
- *Scope:* without a collision every draw is the first one. Committed rows
  are therefore reproduced unchanged: checked before registration with
  `regress` (flag off: A0 `m = 3` `t = 3`, A11 `m = 3` `t = 3`, N `m = 2`
  `t = 4`), and with the flag on for A11 `m = 3` `t = 3` and A0 `m = 3`
  `t = 4` (target 1, equal to the committed `x_R`, abscissae and D). The
  flag changes the collided cells, for example A0 `m = 2` `t = 2` target 1:
  one redraw, a new `x_R`.
- G2 applies to G1 and later runs. Committed collided cells are not re-run.

**More targets per cell: not added.** All 97 completed cells so far have
equal D across their two systems, 76 of them with distinct `x_R`. 96 have
equal `max_cols`; the exception is A1g at `m = 3`, `t = 3`.

**Decision rules.**
- *G1:* §1's σ verdict per arm at `m = 4`, over the arm's completed `m = 4`
  cells. With fewer than 3 the verdict stays inconclusive, its reason is
  recorded, and no further budget or memory change is made without a new
  amendment.
- *G1, secondary (descriptive):* Δlog₂ `max_cols`/Δb per arm, and A11's cell
  `max_cols` over A0's at each shared completed `t`.

**Predictions.**
- G1: A0's σ ≥ 0.5 at `m = 4`. This is F4's registered prediction, still
  undecided.
- No prediction for A11's σ or width (item 3). The draft's width prediction
  rested on widths of systems the budget cut off, so it is dropped.
- A0's `t = 4` systems may hit the memory limit (F3).
