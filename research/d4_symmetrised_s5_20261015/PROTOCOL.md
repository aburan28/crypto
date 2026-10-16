# Protocol: the `D₄`-symmetrised `S₅` (lead 5b), as four weighted unknowns

Registered 2026-10-15, before any measured run of this study.  Lane `d4-s5`
of the five-lane follow-up to `docs/ic/factor-base-encyclopedia/README.md`
§8 item 5.

## 0. What already exists (found before building)

The encyclopedia lists "the `D₄`-symmetrised `S₅` (Bézout count 4096 → 512,
never run)".  That entry is stale.  It was copied from
`research/pdp_root_finder_requirements_20261012` (candidate C8, written
2026-10-12).  C8 was built and run the next day in PR #1785 (commit
`0846e6a12`, an ancestor of the encyclopedia's `a80b15bbc`):
`research/gaudry_quartic_d4_20261013/`, `src/cryptanalysis/gaudry_quartic_d4.rs`.
Its result: `C₄^{D₄} = 2.535–2.552·10¹⁰` `F_p` multiplications at `p = 269`,
one curve, against the `S₄` form's `1.213–1.227·10¹²`.  That is **`47.5–48.4×`,
not the `≈ 630×` the GLV note §9.4 projected**.

C8 wrote the system as **five unweighted unknowns** `(σ₁, σ₂, σ₃, σ₄, π)`
with `π² = σ₄` adjoined.  That gives degrees `4, 4, 4, 4, 2`, Macaulay degree
14 and a `13,412 × 11,628` matrix.  The `630×` projection is
`8^{3.1}`: it assumes the solve costs what a generic system of quotient
dimension 512 costs in this solver (`C ≈ D^{3.1}`, §11.17).  The GLV note
§9.4 named the build that would test that projection, "a weighted Macaulay
solve", and nobody built it.  This study builds it.  It also adds what C8
left open: several field sizes, exhaustive ground truth on tiny fields, the
charged lift-back, and the handover at the cofactor E18 measured.

## 1. The construction (derived; checked in code where marked)

- **Curve and base.**  C8's, unchanged.  The curve is over `F_{p⁴}`
  (`p ≡ 1 mod 4`).  It has `T = (x₀, 0)`, `x₀ ∈ F_{p⁴} \ F_p`, `c ∈ F_p^*`,
  `a = c² − 3x₀²`, `b = −x₀³ − a x₀`, so `f'(x₀) = c²`.  The coordinate is
  `Y = (x − x₀ − c)/(x − x₀ + c)`, with `Y(P + T) = −Y(P)` and
  `Y(−P) = Y(P)`.  The base is `{P : Y(P) ∈ F_p \ {±1}}`, which is E18's
  `τ_T` line `x ∈ x₀ + F_p` and E13/E17's FGHR line at `k = 4`.
- **The group.**  `D₄ = W(D₄) = (Z/2)³ ⋊ S₄`, of order 192.  It acts on the
  four summand coordinates `(Y₁, …, Y₄)` by permutations and by sign changes
  of an even number of them.  The even sign changes are the translations of
  an even number of summands by `T`, which leave the sum unchanged because
  `2T = O`.
- **Invariance.**  `H = S₅(x(Y₁), …, x(Y₄), x(Y_R)) · Π(1 − Y)⁸` has degree
  `≤ 8` in each `Y_i`.  It vanishes on a `D₄`-stable hypersurface, so
  `H∘g = ±H`.  The interpolation onto invariant monomials is the check that
  the sign is `+`: it succeeds only if `H` lies in the invariant span
  (checked on fresh points in every run).
- **Invariants.**  `F[Y]^{D₄} = F[σ₁, σ₂, σ₃, π]` with `σ_k = e_k(Y₁², …, Y₄²)`
  and `π = Y₁Y₂Y₃Y₄` (degrees 2, 4, 6, 4; `192 = 2·4·6·4`).  A monomial
  `σ₁^a σ₂^b σ₃^c π^d` has degree `2(a+b+c) + d` in each `Y_i`, so `H` lies
  on the 105 monomials with `2(a+b+c) + d ≤ 8` (unit-tested).
- **The rewritten system.**  Weil restriction of `H(·; Y_R)` gives four
  equations over `F_p` in `(σ₁, σ₂, σ₃, π)`, each of **weighted degree 8**
  under the weights `w = (2, 2, 2, 1)`.
- **Bézout, derived two ways.**
  1. Weighted Bézout: `Π dᵢ / Π wⱼ = 8⁴ / (2·2·2·1) = 4096 / 8 = 512`.
  2. Orbits: the `S₄` form has `8⁴ = 4096` solutions (§11.16, measured
     4,096 in §11.17), and each is an `S₄`-orbit of ordered tuples.  A
     generic `D₄`-orbit holds `192/24 = 8` `S₄`-orbits, so `4096/8 = 512`.
  C8's five-unknown form gives the same number: `4⁴·2 = 512`.
- **Weighted Macaulay bound.**  The homogenised quotient has Hilbert series
  `(1 − t⁸)⁴ / ((1 − t²)³ (1 − t)²)`.  Its numerator over `(1 − t)` is a
  polynomial of degree `32 − 7 = 25`, so the bound is **weighted degree
  26**.  There are 4,200 columns of weighted degree `≤ 26`, against C8's
  11,628 and the `S₄` form's 40,920.  `π` has weight 1, so `M_π` closes at
  26.  A degenerate eigenspace of `M_π` needs a weight-2 unknown, whose
  multiplication needs degree 27.  The solve then reruns at 27, and both
  runs are charged.
- **Lift-back (charged).**  For each rational `(σ₁, σ₂, σ₃, π)`:
  1. split `T⁴ − σ₁T³ + σ₂T² − σ₃T + π²` for `Y_i²`;
  2. take four square roots in `F_p`;
  3. take the 8 sign vectors with product `π`, which are the 8 even
     `T`-translates;
  4. compute `x(Y)` and lift the points in `E(F_{p⁴})`;
  5. check `±P₁ ± P₂ ± P₃ ± P₄ = R` by group arithmetic for every translate.

  Every `F_p` multiplication is charged, including the `F_{p⁴}` arithmetic
  at `gaudry_quartic`'s rates.

## 2. Arms (same instance, same targets, same counted engine)

| arm | what | status |
|:--|:--|:--|
| **w** (new) | four weighted unknowns, weighted degree 26 (27 on a degenerate eigenspace), `src/cryptanalysis/d4_weighted_s5.rs` | measured |
| **c8** | C8's five-unknown form, degree 14, `gaudry_quartic_d4::solve_d4` | measured, re-run on the same targets |
| **plain** | Gaudry's `S₄` form, `D = 4,096`, Macaulay degree 29 | **out of budget**: its `41,780 × 40,920` `u16` matrix is `3.42 GB`, above this lane's `2.5 GB` cap at every `p` (the matrix does not depend on `p`). Quoted, not re-measured: `C₄ = 1.213·10¹²` (§11.17, `experiments/24_gaudry_quartic_c4.json`), `1.213–1.227·10¹²` (`research/gaudry_quartic_replay_20261013`). §11.17 showed the cost is what a generic solve of four octics costs, `≈ 4 %` from the `D^{3.1}` fit. |

Every run uses one thread (`RAYON_NUM_THREADS=1`), `nice -n 10` and
`prlimit --as=2684354560`.

## 3. Cells

| cell | `p` | instance seed | targets | ground truth |
|:--|:--|:--|:--|:--|
| tiny | 13, 17, 29, 37 | 1 | 6 each (3 planted, 3 uniform), arms w and c8 | brute force over every 4-subset of the base with every sign, and meet in the middle |
| measured | 269, 1033, 4093 | 1 | 8 each (4 planted, 4 uniform), arms w and c8 | meet in the middle |
| w-only | 269, 1033, 4093 | 2 | 12 each, arm w only | meet in the middle |

Planted targets are a signed sum of four distinct base points.  Uniform
targets are random points of `E(F_{p⁴})`.

**Controls.**
- *Positive, solver:* planted generic systems, four equations with random
  coefficients on all 105 monomials of weighted degree `≤ 8` and a common
  root.  Three seeds at each of `p = 269, 1033, 4093`.  Expected: quotient
  512, root found.
- *Positive, pipeline:* the planted targets.  Their planted quadruple must be
  among arm w's decompositions.
- *Negative, invariance detector:* the same curve with a wrong `c`, one whose
  `c²` is not `f'(x₀)`.  `Y` is then not negated by `τ_T`, and the invariant
  interpolation must fail.  Run at `p = 13` and `269`.
- *Ground truth:* brute force (tiny cells) and meet in the middle (all cells).

## 4. Smoke runs used to size the cells (disclosed; not data)

- Generic planted system, `p = 269`, seed 1: `3,979 × 4,200` matrix,
  quotient 512, root found, `1.45·10⁹` multiplications, 1.0 s.
- `p = 13`, seed 1, two targets, both arms:
  - planted target: arm w went to degree 27 on a degenerate eigenspace,
    `2.82·10⁹`;
  - uniform target: `1.17·10⁹` at degree 26;
  - c8: `2.18·10¹⁰` and `2.17·10¹⁰`;
  - all three decomposition sets equal brute force and meet in the middle.
  - Wall time: about 75 s a target for c8 single-threaded, 3–8 s for w.

## 5. Statistics

- `C_w`: arm w's `F_p` multiplications per solve, all phases, reruns
  included.  This is the same unit and the same phase set as C8's `C₄` and
  §11.17's: Weil restriction, triangularisation, elimination, normal forms,
  characteristic polynomial, roots, eigenvectors, verification, split.
- `L_w`: lift-back multiplications per target.
- `R_plain = 1.213·10¹² / C_w` per target, against the quoted plain `C₄`.
  Reported as the median and min–max per size.
- `R_c8 = C_c8 / C_w`, per target and paired.
- Per target, the counts:
  - rational invariant solutions;
  - solutions whose `Y_i²` split;
  - solutions whose `Y_i` are rooted in `F_p`;
  - candidate representatives and verified decompositions;
  - degree-27 reruns.

## 6. Predictions and decision rule (fixed now)

- **W1 (structure).**  Quotient dimension 512 on every target whose matrix
  closes, with 4,200 columns at degree 26.  A target with a different
  dimension is reported as degenerate, as C8's target 8 (dimension 517) was.
- **W2 (exactness).**  Arm w's decompositions equal brute force on every
  tiny target and meet in the middle on every target.  Every planted
  quadruple is found.  Any disagreement voids arm w's cost numbers.
- **W3 (flat in `p`).**  The median `C_w` at the three measured sizes agrees
  within `±25 %`, since only root finding depends on `p`.
- **W4 (the projection).**  The generic `D^{3.1}` scaling projects
  `C ≈ 1.9·10⁹`, `R_plain ≈ 630`.
- **Decision on lead 5b** (median `R_plain` at every measured size; W2
  must hold):

  | outcome | condition |
  |:--|:--|
  | **projection met** | `315 ≤ R_plain ≤ 1260` (within `2×` of 630) |
  | **short of projection** | `R_plain < 315` |
  | **above projection** | `R_plain > 1260` |
  | **accounting** | the median `L_w` exceeds `10 %` of the median `C_w`: the lift-back eats the gain, and the net ratio `1.213·10¹² / (C_w + L_w)` is the headline |

- **Against C8:**

  | outcome | condition |
  |:--|:--|
  | **advances on C8** | paired `R_c8 ≥ 2` on every paired target |
  | **no gain over C8** | median `R_c8 < 1.2` |
  | **mixed** | otherwise |

- **Handover (derived, extrapolation; not a measurement).**  E18's `τ_T` row
  (`RESEARCH_GLV_INVARIANT_FACTOR_BASES.md` §9.3) is
  `p* = y·(cols/p)·C·√h / (S_rho·(1 − r))`.  It uses the measured
  `y = 106` residuals per relation, `r = 0.341` against the negation walk,
  `h = 4`, and `S_rho = 90.8`.  At `C₄ = 1.213·10¹²` it gives
  `n* = 2^{157.9}`.  Since `n* ∝ C⁴`, the solve moves it to
  `157.9 − 4 log₂(1.213·10¹² / C)` bits.  It is reported for `C_w`,
  `C_w + L_w` and C8's `C`, beside:
  - the plain prime-order figures `2^{151.1}` (§11.19, unfolded walk) and
    `2^{156.1}` (negation walk);
  - E18's `ψ`-folded `2^{146.0}`.

  `r∞` does not depend on the solve, so it is unchanged.

**Class if any of this holds: engineering** (a constant of one phase).  No
end-to-end `S` is measured, and none will be claimed.
