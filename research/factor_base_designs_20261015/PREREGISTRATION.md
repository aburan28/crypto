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
