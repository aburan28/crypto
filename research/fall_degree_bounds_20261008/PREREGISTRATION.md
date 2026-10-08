# Preregistration: proving or refuting the fall-degree bounds for `m = 2` index calculus

**Written:** 2026-10-08, while the E1, E1b and E2 runs were in flight and
before any of their output was read. The bounds being tested are ledgered
in `research/notes/index-calculus/RESEARCH_FALL_DEGREE_BOUNDS.md`.

**What was seen beforehand, disclosed:**

- The `lfd.sage` cells for `n = 8 … 20`, `TABLE.md`.
- An exploratory fit, `fit.sage`, of the constant term of the
  Kosters–Yeo trace relation on 189 draws at `n = 6 … 10`. Its seeds are
  `31000·n + s`, disjoint from every seed below.
- Three `n = 12` syzygy draws (seeds `60000–60002`) that motivated E1b.
  E1b uses seeds `5000·n + s`, so `n = 12` reuses `60000–60007`; those
  three draws are excluded from scoring.

## System

- Curve `y² + xy = x³ + b` over `F_{2^n}`, minimal-weight modulus.
- Target `x_3`: uniform random for `rand` cells. For `sat` cells it is a
  root of `S_3(v_1, v_2, T)` with `v_1, v_2 ∈ V`.
- Factor base `V ⊂ F_{2^n}` of dimension `n' = ⌊n/2⌋`. It is either
  random (`random`) or `F_{2^{n/2}}` (`subfield`).
- The system `F`: the `n` Weil-descended equations of `S_3(X_1, X_2, x_3)`
  with `X_i ∈ V`, in `N = 2n'` Boolean unknowns.
- `d_F`: the Huang–Kosters–Yeo last fall degree, computed by the mutant
  closure in the Boolean ring (definition and equivalence argument in the
  ledger, §0).
- `d_ff`: the operational first fall degree.

## E1: the trace relation, exactly

**Hypothesis H1 (literature-derived, sharpened).** For every draw with
`x_3 ≠ 0`, the polynomial

```
L = Σ_j Tr(v_j)·(y_{1j} + y_{2j}) + Tr(b / x_3²)
```

lies in `span_{F_2}(F)`, the degree-2 Macaulay space with no multipliers.
Here `v_j` is the basis of `V`, `y_{ij}` are the bits of `X_i`, and `Tr`
is the absolute trace.

Consequences, which become theorems once H1 is proved (§P1):

- **C1.** If `Tr|_V ≢ 0`, then `d_ff = 2`.
- **C2.** If `Tr|_V ≡ 0` and `Tr(b/x_3²) = 1`, then `1 ∈ span(F)`. The
  target is then not decomposable over `V`, and `d_F ≤ 2`.
  - The subfield base at `n = 2n'` always has `Tr|_V ≡ 0`.
  - So C2 predicts that about half of all subfield targets are refuted
    at degree 2 by one trace bit.
- **C3.** If `Tr|_V ≡ 0` and `Tr(b/x_3²) = 0`, then `L = 0` and the trace
  gives no fall. Any `d_ff = 2` there has another cause.

**Run.** `sage exp.sage trace n n/2 8 {random,subfield} rand 6` for
`n ∈ {8, 10, 12, 14, 16}`, and the random base at odd `n ∈ {9, 11, 13, 15}`.

**Falsified if** any of these holds:

- some draw has `L ∉ span(F)`;
- some draw with `Tr|_V ≢ 0` has `d_ff ≠ 2`;
- some subfield draw with `Tr(b/x_3²) = 1` is not refuted at degree 2.

**Registered side prediction.** Among draws with `Tr|_V ≢ 0`, the
degree-≤1 part of `span(F)` has dimension exactly 1 in every draw at
`n ≥ 12`.

## E1b: the "bounded-defect" syzygy is the trace relation

`RESEARCH_FFD_PROOF_COMPLEXITY.md` §3.4 found, in EXP-J, exactly one
excess degree-3 syzygy `ℓ·(Σ_{i∈S} f_i) ≡ 0`, with `ℓ` symmetric. It
proposed this as the antecedent of its conditional lower bound.

**Hypothesis H1b.** That syzygy is the trivial Boolean identity
`(L + 1)·L = 0` applied to the trace relation. Concretely:

- the space of degree-3 linear syzygies `Σ_i ℓ_i f_i = 0` has dimension
  exactly 1;
- its multiplier is `ℓ = L + 1`, and `Σ_{i∈S} f_i = L`.

**Run.** `sage exp.sage syz n n/2 8 {random,subfield} rand 0` at
`n ∈ {12, 14, 16, 18, 20}`, and the random base at odd `n ∈ {13 … 19}`.

**Falsified if** a random-base draw with `Tr|_V ≢ 0` and `n ≥ 12` has
`syz_dim ≠ 1`, or its multiplier is not `L + 1`.

**Consequence if it holds.** The bounded-defect lemma has a proof of
existence: Kosters–Yeo Cor. 4.11 plus the Boolean identity. Its
uniqueness half reduces to C1's dimension statement.

## E2: is the low `d_F` structural or generic?

**Controls.** Random Boolean quadratic systems with the same shape: `n`
equations in `n` unknowns. They come in two kinds:

- `plain`;
- `planted_linear`: the last equation is the sum of the others plus a
  random affine-linear form. This mimics the trace fall.

Each kind is run unsatisfiable-or-random (`rand`) and with a planted root
(`sat`). The run is `sage exp.sage control n n/2 4 kind {rand,sat} 6` at
`n ∈ {8, …, 18}`.

**H2a.** `d_F(Semaev) ≤ d_F(planted_linear control)` holds at every `n`,
comparing cell medians with the same satisfiability.

**H2b.** The gap `d_F(control) − d_F(Semaev)` is nondecreasing in `n`.

- The literature anchors this. The Kousidis–Wiemers F4 top degree, which
  is at least `d_F`, is 5 at `n = 48`. The semi-regular `D_reg` there is
  8.
- So a generic control must eventually sit above Semaev.

**Falsified if** Semaev's median `d_F` exceeds the control's at any `n`
(H2a fails), or the gap shrinks between consecutive `n` (H2b fails).
Equality at every `n ≤ 18` would not falsify H2a. It would mean the
structure becomes visible only beyond toy size, and E3 would decide.

## E3: scaling of `d_F` for `m = 2` (needs E5's engine)

Cells: `n ∈ {22, 24, 26, 28, 32}`, random base, 8 draws, rand and sat.

- **P3a.** `d_F ≤ 5` for every `n ≤ 32`. This follows from Kousidis–
  Wiemers' F4 data together with `d_F ≤ sd` (Caminata–Gorla). A violation
  would mean a definitional or engine bug, not a discovery.
- **P3b, the scientific question.** Between `n = 18` and `n = 32`,
  `d_F` steps from 4 to 5 at most once.
- Two models compete:
  - `LOG`: `d_F = a + b·log₂ n`, fitted on `n ≤ 22`;
  - `LIN`: `d_F = a + n/c`, fitted on the same points.
- The model with the smaller absolute error on `n ∈ {24 … 32}` is
  preferred. A tie is reported as a tie.

**Kill condition for "d_F bounded".** `d_F` reaches 5 by `n = 32` *and*
the step widths `n(3→4)`, `n(4→5)` grow less than geometrically.

## E4: `m = 3` (chained `S_3` with one auxiliary variable; needs E5)

The cells are the `dreg_ladder` cells `(n, ℓ) ∈ {(5,2), (7,2), (5,3), (7,3)}`,
on the same draws as their committed runs.

- **P4a (THEOREM check).** `d_F ≤ D*` on every draw, where `D*` is the
  ladder's plain refutation degree. A violation is a bug.
- **P4b.** `d_F < D*` on at least half the draws. The mutant closure is
  then strictly cheaper than the ladder's plain tower.
- **P4c.** The Kousidis–Wiemers bound `d_ff ≤ m² − m + 1 = 7` holds for
  the direct `S_4` descent at `n ≤ 15`.

## E5: native engine (prerequisite for E3, E4)

Port the closure `V_{F,c}` to Rust on the existing M4RI path. The port
must reproduce the per-degree dimension of `V_{F,c}` on every Sage cell
in `runs/` and `exp/`, with exact equality. Any disagreement blocks E3
and E4.

## Proof obligations

| id | statement | route | status at registration |
|---|---|---|---|
| P1 | H1 as a theorem for `a_1 = 1, a_3 = 0` curves | Kosters–Yeo Cor. 4.11 with `a_2` chosen so that `x_3` lies on `E_{a_2}`. Their constant is `Tr(x_3 + a_2)`, and `Tr(x_3 + a_2) = Tr(b/x_3²)` exactly when `(x_3, ·) ∈ E_{a_2}(F)`. `S_3` does not depend on `a_2`. | proof sketch, ledger §5 |
| P2 | C1–C3 | linear algebra from P1 | follows from P1 |
| P3 | `dim(span(F) ∩ R_{≤1}) = 1` when `Tr|_V ≢ 0`, `n' ≥ 6` | genericity of the remaining `n−1` quadratic parts over a random `V` | open, E1 and E1b evidence |
| P4 | `d_F ≥ 3` for random `V` with `Tr|_V ≢ 0`, `n ≥ 10` | show that `V_{F,2}` misses some degree-2 element of the ideal | open |
| P5 | `d_F = ω(1)`, or a uniform bound | PC-degree lower bound on the trace-quotiented system (Alekhnovich–Razborov immunity), or an explicit low-degree certificate | open, decided empirically by E3 at toy size only |

## Scoring

`score.py` (to be written after the runs, from committed JSONL only)
reports each hypothesis as **held**, **falsified**, or **not testable**,
quoting the cell that decides it. Results go in `RESULTS.md`. This file
is not edited after commit except by dated amendments.
