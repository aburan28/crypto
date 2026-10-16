# Protocol: is the class-level `full_success` deficit a granularity artefact?

Frozen before any study data. Follows up the "also open" item of
[`problems/smith-kernel-isogeny-class-walk.md`](../../problems/smith-kernel-isogeny-class-walk.md):
why do 560 of 1,855 (`p = 7`) and 1,210 of 7,431 (`p = 11`) tractable
`L`-classes hold no `full_success` model, so that the class-level share `C`
lies below the independence null of
[`research/smith_prop24_native_20261014`](../smith_prop24_native_20261014/REPORT.pdf)?

## Hypothesis under test

The 2026-10-14 null (`smith_prop24::summarize_census`) scatters successes
uniformly over *models* (monic octics with zero `x^7` coefficient, and their
twists by the least non-residue) within each `#S(H)` stratum. But
`full_success` is a property of the curve up to `F_p`-isomorphism, and one
isomorphism class of curve accounts for many census models. A small
`L`-class may hold only one or two isomorphism classes, so it receives one or
two independent coins, not dozens. **H_gran:** with the isomorphism class as
the unit, the observed `C` is what independent coins give.

## Definitions

- *Model*: `(f, t)` with `f` a squarefree monic octic over `F_p` with zero
  `x^7` coefficient and `t ∈ {1, c}`, `c` the least non-residue; the curve is
  `y^2 = t f(x)`. These are exactly the 2,823,576 / (p = 11) census curves of
  the 2026-10-14 study.
- *Isomorphism class*: two models are equivalent when
  `t' f'(x) = μ^2 (γx+δ)^8 · t f((αx+β)/(γx+δ))` for some
  `[[α,β],[γ,δ]] ∈ GL_2(F_p)` and `μ ∈ F_p^*`. Orbits are enumerated
  exactly by applying all of `PGL_2(F_p)` to the binary octic form,
  normalising the leading coefficient to `1` or `c` by a square, and keeping
  images with zero `x^7` coefficient.
- Per isomorphism class: `L`-key `(a1, a2, a3)`, stratum `#S(H)`,
  `full_success`, and orbit size (number of census models).

## Measurements

1. Analyse every census model with `smith_prop24::analyze_curve` (as the
   census did) and group by isomorphism class.
2. **Controls (must all pass, else the run is invalid):**
   - C1: per `L`-class, the model counts, tractable counts and
     `full_success` counts equal `census_p{7,11}_classes.csv` exactly.
   - C2: `full_success`, stratum and `L`-key are constant on every
     isomorphism class (every model checked).
   - C3: orbits partition the census models (sizes sum to the census total).
3. **Primary null (iso-level):** within each stratum, permute the
   `full_success` labels uniformly among isomorphism classes (count of
   successful classes fixed); recompute `C` (share of models in `L`-classes
   with a `full_success` member). 999 permutations, seed `20261015`.
4. **Secondary null:** the same, stratified by (stratum, orbit size), since
   orbit size reflects the automorphism group.
5. Descriptive: per tractable `L`-class, the number `k` of tractable
   isomorphism classes; the share of `L`-classes with no `full_success`
   member as a function of `k`, against `(1 − q_s)^k` with `q_s` the
   per-stratum iso-level success rate.

## Decision rule (per prime, primary null)

- Observed `C` within the central 95% of the null → **"granularity explains
  the deficit"**.
- Observed `C` below the 2.5% quantile → **"residual deficit beyond
  granularity"** (a structural class effect remains to be named).
- Above the 97.5% quantile → **"excess over granularity"**.

Overall: the common verdict if both primes agree, else "prime-dependent".

## Scope and limits

Exact enumeration at `p = 7, 11` only, degree-8 models only (curves with a
rational non-Weierstrass point at infinity; curves whose eight `P^1(F_p)`
points are all Weierstrass are absent from the census and from this study).
It does not price or build any isogeny walk, and says nothing at
cryptographic sizes or in characteristic 2. A pass of H_gran says the deficit
needs no structural explanation at these primes; it does not prove that none
exists.
