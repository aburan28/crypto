# Is the `m = 4` tree-split set cut out at low degree? A cut-out degree screen

Registered before any registered cell runs. §5 lists everything that was run before
registration and what was seen.

## 1. Why

- **The route.** The top pick of the 2026-09-30 deep-research shortlist, and the one
  reopening condition of [DECISION-20260929-stop.md](../ic_candidate_tournament_20260915/campaign_20260916/DECISION-20260929-stop.md)
  that no one has tested: "a decomposition that is not a table and not the chained `S₃`
  system".
- **The tree split.** Write a 4-decomposition `R = P₁ + P₂ + P₃ + P₄` as
  `(P₁ + P₂) + (P₃ + P₄) = R`. With `e = x(P₁ + P₂)` and `f = x(P₃ + P₄)`, both `e` and
  `f` lie in `T2 = {x(P + Q) : P, Q ∈ F}`, and `S₃(e, f, x_R) = 0`.
- **Why it might not be a table.** In characteristic 2,
  `S₃(e, f, x_R) = (E₂ + E₁·x_R)² + x_R·E₂ + b` with `(E₁, E₂) = (e + f, e·f)`. Squaring is
  GF(2)-linear, so this is **affine-linear over GF(2)** in the bits of `(E₁, E₂)`: for each
  target it defines an affine subspace `A_R` of dimension `n` in `GF(2)^{2n}`. A
  decomposition exists exactly when `Z ∩ A_R ≠ ∅`, where `Z = Sym²(T2) = {(e + f, e·f)}`.
  All the difficulty sits in `Z`, and `Z` does not depend on the target.
- **What would make it a lever.** Suppose `Z` (or `T2`) is cut out by equations whose
  degree stays flat as `ℓ` grows. Then `Z ∩ A_R` is a low-degree system in `n` unknowns
  for every target. Its solving degree would not climb with `ℓ` the way the chained
  system's did (`ĉ ≈ 1` at `m = 4`,
  [head-engine audit](../ic_m4_head_engine_20260929/RESULTS.md)).
- **What would close it.** If `T2` and `Z` need as high a degree as a random set of their
  size, then the tree split has no algebraic shortcut. It is then a table method, and the
  survey's §3.5 bound (`r^{j/(2j−1)}`) applies.
- **The prior is closure.** `T2` is the image of `2ℓ` bits under a map of algebraic degree
  well above 1 in `n` bits, and such images usually look random at low degree.

## 2. Objects and controls

Curve `K_a : y² + xy = x³ + a·x² + 1` over `GF(2^n)`, with the field polynomial chosen by
`find_irreducible(n)`.

- **Factor base.** `V` is a random `ℓ`-dimensional GF(2)-subspace of `GF(2^n)`, seeded.
  `F` is the set of points with `x ∈ V \ {0}`.
- **`T2`** is `{x(P + Q) : P, Q ∈ F, P + Q ≠ O}`, a subset of `GF(2)^n`.
- **`Z`** is `{(e + f, e·f) : e, f ∈ T2}`, unordered with `e = f` allowed, a subset of
  `GF(2)^{2n}`. The map is injective on unordered pairs, so `|Z| = |T2|(|T2| + 1)/2`.

Controls have the same size and the same ambient dimension:

- **`R`** is a uniformly random subset.
- **`ZR`** (for `Z` only) is `Sym²` of a uniformly random set of `|T2|` field elements. It
  separates what `Sym²` does on its own from what the curve adds.

## 3. Metric

For a set `S ⊂ GF(2)^N` and a degree `D`, form the evaluation matrix of all multilinear
monomials of degree `≤ D` (`M(D)` columns) at the points of `S`. Its GF(2) rank is `h_S(D)`.

- **`deficit(D)`** is `min(|S|, M(D)) − h_S(D)`. For a random set it is 0; a positive
  value means `S` satisfies extra polynomials of degree `≤ D`.
- **Sampled cut-out degree `D_s`** is the smallest `D ≤ 6` at which none of `K = 2000`
  points, drawn uniformly from `GF(2)^N \ S`, satisfies every degree-`≤ D` polynomial
  vanishing on `S`. A point satisfies them all exactly when its monomial vector lies in
  the row space.
- **What `D_s` resolves.** It is a screen, not an exact cut-out degree. With 2000 samples it
  detects a degree-`D` closure larger than `S` by at least 0.15% of the space with
  probability 95%, and can miss a smaller excess.
- **For a random set**, `D_s` is the counting threshold: the smallest `D` with `M(D)`
  comfortably above `|S|`. A structured set can be cut out below it.
- **Censoring.** A degree is skipped when `M(D) > 90,000` columns. An arm that stops before
  reaching `D_s` is **censored** (`D_s > D`), and a cell killed by its CPU or memory limit
  (`run.sh`) is censored too. **An `ℓ` with any censored arm is excluded from the
  verdict** and listed. Censoring is never evidence either way.

## 4. Cells and decision rule

[run.sh](run.sh) runs every cell once; [analyze.py](analyze.py) is the readout.

| object | `n` | `a` | seeds | `ℓ` | arms |
|:--|:--|:--|:--|:--|:--|
| `T2` (primary) | 17, 19 | 0, 1 | 20260930, 20261001 | 4, 5, 6, 7, 8 | curve, `R` |
| `Z` (secondary) | 13, 15 | 0, 1 | 20260930 | 4, 5 | curve, `ZR`, `R` |

That gives eight `T2` series and four `Z` series. A **series** is one `(object, n, a, seed)`.
For each series and each `ℓ`, the **gap** is `D_s(curve) − min over the controls of D_s`.

**`T2` verdict:**

- **ALIVE** if the gap at the largest uncensored `ℓ` is `≤ −2` in **every** series, **and**
  the median over series of the OLS slope of `D_s(curve)` on `ℓ` is `≤ 0.15`.
- **CLOSED** if the gap is `≥ −1` at every uncensored `ℓ` in every series.
- **INCONCLUSIVE** otherwise.

**`Z` verdict:** the same rule, without the slope clause. Two values of `ℓ` do not make a
slope, and a `Z` result on its own neither opens nor closes the route.

Reported, but no decision rests on them: the full deficit and false-zero profile of every
arm, and the curve arm's false-zero fraction at `D_s(R)`, the "closure excess".

### What each verdict would mean

- **ALIVE** means the tree-split set has a low-degree description at `n ≤ 19`, `ℓ ≤ 8`.
  - It motivates a pre-registered solving-degree audit of `Z ∩ A_R`, against the
    `c* = 0.25` gate of the `m = 4` audits.
  - It is **not** an exponent claim.
- **CLOSED** means that at these sizes neither `T2` nor `Z` is cut out below the degree a
  random set of its size needs.
  - The tree split then has no low-degree shortcut, and the reopening condition stays
    unmet for this decomposition.
  - It says nothing about other models of `Z ∩ A_R`: linearisation with auxiliary
    variables, or a structured intersection algorithm. It also says nothing about larger
    `n` or `ℓ`.

## 5. What ran before registration

- **Correctness checks**, run on non-registered data:
  - the field inverse;
  - commutativity and associativity of the group law on 30 random points at each of
    `(a, n) = (0, 13), (1, 17), (0, 7)`;
  - `#K_0(GF(2^7)) = 116`, which matches the Frobenius trace recurrence, with `116·P = O`
    on five points;
  - `S₃(x₁, x₂, x(P + Q)) = 0` in the affine `(E₁, E₂)` form on up to 200 pairs per curve,
    with no failures.
- **Linear-algebra checks:**
  - a dimension-8 affine subspace of `GF(2)^20` gives `D_s = 1`;
  - the full quadric `x₀x₁ + x₂ = 0` in `GF(2)^14` gives `D_s = 2` with deficit 3 at
    `D = 2`, which is the expected three multilinear multiples;
  - 1000 random points in `GF(2)^20` give `D_s = 3`, the counting threshold.

  These were run first with a numpy elimination and then with `examples/gf2_span.rs`. Both
  gave the same results.
- **Smoke, `n = 11` (not a registered size), seed 20260930, `a = 0`. Seen in full:**
  - `T2` at `ℓ = 2, 3, 4`: `|T2|` = 1, 4, 36, and `D_s` equal to `R`'s at every `ℓ`.
  - `Z` at `ℓ = 2, 3, 4`: `|Z|` = 1, 10, 666, and `D_s` = 1, 1, 3 against `ZR` 1, 2, 3
    and `R` 1, 1, 3.
  - At these sizes the sets are too small to separate from random, which is why the
    registered `ℓ` starts at 4 for `T2` and the registered `n` sit above 11.
- **Timing runs (not registered cells):**
  - `T2` at `n = 21`, `ℓ = 7`, seed 1: the output was discarded unseen. It took 6 s.
  - `Z` at `n = 11`, `ℓ = 5`, seed 1, curve arm only. `|F| = 28`, `|T2| = 178`,
    `|Z| = 15,931`.
    - At `D = 4`: rank 9,109, full column rank.
    - At `D = 5`: rank 15,931, full row rank, with **10 of 2000 samples in the span**.
    - That is, this curve set was **not** cut out at `D = 5`, where a random set of its
      size is expected to be.
- **Written after that timing run.** The CLOSED clause counts a curve arm that sits
  **above** its controls (gap `> 1`) as closed, not inconclusive. The lever needs a
  **lower** degree, so a set that is harder to cut out than random cannot be it. The
  shortlist's draft rule said "within ±1 of the controls". This clause is disclosed as an
  extension of that draft.

## 6. Implementation and scope

- **Code.** [cutout.py](cutout.py) builds the sets and the evaluation matrices.
  `examples/gf2_span.rs` does the GF(2) elimination:
  `cargo build --release --example gf2_span`.
- **Environment.** Python 3.11 with numpy 2.4.6.
- **Not timed.** The metric is a degree, so the AGENTS.md §10 isolation rule for timed runs
  does not apply. The per-degree `secs` fields are logged and are not measurements.
- **Scope.** Two Koblitz curves, `n` 13–19, `ℓ ≤ 8`, one random subspace per seed, and a
  sampled cut-out degree capped at 6 and at 90,000 monomials. Anything beyond that is
  untested.
