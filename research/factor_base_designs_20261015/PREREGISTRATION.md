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

(none yet)
