# Protocol: which invariant finer than the `L`-polynomial carries Smith's `full_success`?

Frozen before any invariant was computed. Follows up
[`research/smith_iso_granularity_20261015`](../smith_iso_granularity_20261015/PROTOCOL.md)
(PR #1928), whose verdict at `p = 7` and `p = 11` was *residual deficit
beyond granularity*: with the `F_p`-isomorphism class of curve as the unit,
the observed class share `C` (0.4304 / 0.4750) lies below all 999
iso-level permutations (null 95% bands [0.4545, 0.4611] / [0.4858, 0.4871]),
so `full_success` clusters inside `L`-classes beyond what independent coins
per isomorphism class give. A post-hoc probe found `v2(#J(F_p))` (a function
of `L`) does not explain it.

## Question

Is there a Galois-module invariant of `J = Jac(H)`, finer than the
`L`-polynomial, such that, once the iso-level permutation null is
conditioned on it, the observed `C` is what independent coins give?

## Frozen inputs

- `research/smith_iso_granularity_20261015/data/iso_p{7,11}.csv.gz`
  (rows `rep_code,rep_twist,a1,a2,a3,stratum,full,size`, one per
  isomorphism class; 34,227 / 324,562 rows), read-only. The curve of a row
  is `y^2 = t f(x)`, `f` monic of degree 8 with zero `x^7` coefficient whose
  coefficients `f_0..f_6` are the base-`p` digits of `rep_code`, and `t = 1`
  or the least non-residue `c` (`rep_twist = 1`).
- `#J(F_p) = L(1) = 1 + a1 + a2 + a3 + p a2 + p^2 a1 + p^3`, the convention
  of the frozen `summarize` code.
- Seeds: permutations `20261015` (as in PR #1928); relation sampling for the
  class group `SplitMix64(rep_code * 2 + rep_twist)`; placebo buckets
  `SplitMix64(20261016)`.

Invariants are computed for every *tractable* isomorphism class
(`stratum != 0`, 16,823 / 160,010 classes). Non-tractable classes have
`full = 0` and never move in any null.

## Candidates

Let `σ` be the Frobenius permutation of the eight Weierstrass points (the
roots of `f`; all eight are finite in these models).

**X1. `J[2]` (cycle type of `σ`).** `J[2]` is the even subsets of the eight
roots modulo the full set, and Frobenius acts through `σ`, so the
`F_2[Frob]`-module `J[2]` with its Weil pairing is a function of the
`S_8`-conjugacy class of `σ`, i.e. its cycle type. The symplectic class of
Frobenius in `Sp_6(F_2)` is therefore *at most* as fine as the cycle type
(proposed in the task as possibly finer; it cannot be, for a hyperelliptic
curve). X1 is the finest `J[2]`-level invariant: the cycle type. The stratum
`#S(H)` is a function of it (Smith Lemma 1).

**X2. 2-Sylow of `J(F_p)`.** The isomorphism type of `J(F_p)[2^∞]`, the list
of exponents `v_1 >= v_2 >= ...` of `⊕ Z/2^{v_i}`. Not a function of `L`
(only `Σ v_i = v2(#J(F_p))` is).

**X3. Frobenius on `J[4]` (fixed-point profile).** With `m` the order of
`σ`, Frobenius `F` acts on `J[4] ≅ (Z/4)^6` with order dividing `2m` (`F^m`
is trivial on `J[2]`, so `F^m = 1 + 2N` and `F^{2m} = 1` on `J[4]`). X3 is
the cycle type of `σ` together with `(#J(F_{p^k})[4])_{k = 1..2m}`. These
counts are `#ker(F^k - 1 | J[4])`, so X3 determines the orbit structure of
Frobenius acting on the 4,096 points of `J[4]`. The characteristic
polynomial of Frobenius mod 4 is a function of `L`; X3 is not.

**X3+. Frobenius on `J[4]` (halvability pattern).** The cycle type of `σ`
with, for every `k = 1..2m`, the *set* `H_k ⊂ J[2]` of 2-torsion points
rational over `F_{p^k}` that are halvable over `F_{p^k}`, put in canonical
form under relabelings of the roots that commute with `σ` (lexicographic
minimum over the centralizer of `σ` in `S_8`). X3+ refines X3. It is an
operational invariant of the Galois module `J[4]`; whether it is a
complete invariant of the `GSp_6(Z/4)`-conjugacy class of Frobenius is not
established here and is not claimed.

## Methods (native Rust, no Python)

- **Cycle type** from the factor degrees of `f` (`smith_prop24::factor_squarefree`).
- **X2 by the class group of the function field.** Generators: all places
  of degree `<= 3` (over `x = ∞`: `∞±` when `t` is a square, else one place
  of degree 2; over monic irreducible `a` of degree `<= 3`: split, ramified,
  or inert places, inert only for `deg a = 1`). Relations: `div(a)` for every
  such `a`, and `div(y - v)` for polynomials `v` (`deg v <= 3`, and
  `v = ±√t x^4 + w`, `deg w <= 3`, when `t` is a square) whose norm
  `v^2 - t f` is squarefree with all factors of degree `<= 3`. Places of
  degree `<= 3` generate `Pic(H)` (for `E` of degree 0 and an effective
  rational `D` of degree 3, `E + D` is linearly equivalent to an effective
  divisor of degree 3 by Riemann–Roch). With `e = v2(#J(F_p))` and
  `E = e + 1`, compute the Smith form of the relation lattice over
  `Z/2^E`. Because `Pic = Z ⊕ J(F_p)`, the quotient surjects onto
  `Z/2^E ⊕ J(F_p)[2^∞]`; the result is **certified** only when its order
  is exactly `2^{E+e}` (then the surjection is an isomorphism). Relations
  are added until certified (cap: 40 rounds; uncertified classes are
  recorded and count against control K3).
- **X3 / X3+ by 2-descent** over `K = F_{p^k}`. For a `K`-rational
  `T_S ∈ J[2]` (`S` an even set of roots with `σ^k S ∈ {S, S^c}`) the
  descent map to `L_K^* / K^* L_K^{*2}` (`L_K = K[θ]/(f)`) has
  component at root `e_j`:
  `δ_j(T_S) = ∏_{i ∈ S, i ≠ j} (e_i - e_j) · (f'(e_j) if j ∈ S)`
  (evaluation against fixed uniformizers: `y` at Weierstrass points, `1/x`
  at infinity; the ambiguity from generalized Weil reciprocity is a single
  constant in `K^*`). `δ(T_S)` is trivial iff all quadratic characters
  `χ_{K(e_j)}(δ_j)` are `+1`, or all equal `(-1)^{[K(e_j):K]}`. We count
  `T` with `δ(T)` trivial as halvable and set
  `#J(K)[4] = #J(K)[2] · #{T ∈ J(K)[2] halvable}`. The kernel of `δ` on
  `J(K)/2J(K)` may be nontrivial for even-degree models; control K4 tests
  this at `k = 1` against X2 on every class. Roots are taken in the
  splitting field `F_{p^m}` (`m <= 15`).

## Controls (all must pass, else the affected candidate is invalid)

- **K1** (reproduction): the unconditioned iso-level null (strata only,
  seed 20261015) reproduces `C`, `lo`, `median`, `hi` of
  `summary_p{7,11}.json` exactly.
- **K2** (class invariance): for every tractable isomorphism class, X1,
  X3 and X3+ recomputed on a second orbit member (the last element of the
  sorted orbit) equal the representative's; for X2, on a fixed sample of
  1,000 classes per prime (seed 20261015).
- **K3** (X2 certification): every tractable class is certified;
  `Σ v_i = v2(#J(F_p))`; the number of `v_i` equals `dim J(F_p)[2]`
  computed from `σ` (count of `T` with `σT ∈ {T, T^c}`).
- **K4** (descent vs class group): at `k = 1`, the descent count of
  halvable `T` equals `2^{#{i : v_i >= 2}}` from X2 on every class.
- **K5** (profile sanity): for each `k`, `#J(F_{p^k})[4]` divides
  `2^{v2(#J(F_{p^k}))}` (with `#J(F_{p^k})` from `L`), and is
  `#J(F_{p^k})[2]` when `v2(#J(F_{p^k})) = log2 #J(F_{p^k})[2]`.

If K4 fails, X3/X3+ are reported as invalid (the descent kernel is then
not trivial) and the failure pattern is recorded; X1 and X2 stand on their
own controls.

## Nulls and statistics (per prime, per candidate X)

- **N_X:** permute `full` labels among iso classes within cells
  `(stratum, X)` (non-tractable classes: one cell); 999 permutations,
  RNG `SplitMix64(20261015)` freshly seeded per candidate; recompute `C`
  (share of models in `L`-classes with a `full_success` member).
- **Placebo P_X (L-class proxy of equal cardinality):** each `L`-class gets
  a pseudo-random bucket `b ∈ {1..n_s}` per stratum `s`, where `n_s` is the
  number of distinct X values among stratum-`s` tractable classes; permute
  within `(stratum, b)`. A candidate that "explains" `C` only because it
  partitions iso classes along `L`-class lines also passes its placebo.
- **Gap share:** `g_X = (med N_0 - med N_X) / (med N_0 - C)` with `N_0` the
  unconditioned null of K1.
- Descriptive, not part of the verdict: number of cells; share of tractable
  iso classes in single-class cells; per-cell success rates.

## Decision rule (per candidate, per prime)

- `lo(N_X) <= C <= hi(N_X)` and `C < lo(P_X)` → **"explains"**.
- `lo(N_X) <= C <= hi(N_X)` and `C >= lo(P_X)` → **"confounded with L-class"**
  (the conditioning is too close to the `L`-class to discriminate).
- `C < lo(N_X)` and `g_X >= 0.5` → **"partial"**.
- `C < lo(N_X)` and `g_X < 0.5` → **"does not explain"**.
- `C > hi(N_X)` → **"over-explains"**.

Overall per candidate: the common verdict if both primes agree, else
"prime-dependent". An invalid control makes the candidate's verdict
"invalid" at that prime.

## Stop conditions and cost

Exact enumeration at `p = 7, 11` only; 2 threads. If X2 certification has
not finished for `p = 11` after 6 hours wall clock, record the uncertified
remainder, report X2 at `p = 11` as "incomplete", and do not substitute a
smaller sample silently.

## Scope and limits

Degree-8 models only (as in PR #1928); no isogeny is built; nothing is
claimed at cryptographic sizes or in characteristic 2. "Explains" means
"the conditioned null admits the observed `C`", not that the candidate
causes `full_success`; a negative bounds only these candidates at these
primes.
