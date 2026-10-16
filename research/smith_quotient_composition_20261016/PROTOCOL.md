# Protocol: is Smith's `full_success` decided by the make-up of the target class?

Frozen before any study data. Stage 1 of the staged follow-up to
[`research/smith_class_invariant_20261015`](../smith_class_invariant_20261015/REPORT.pdf)
(PR #1955), which found that no `J[2]`, 2-Sylow or `J[4]` invariant of
`J_H` explains why `full_success` clusters inside `L`-classes at
`p = 7, 11`.

## Background (from Smith 2009, §4–§6, arXiv:0806.2995)

For a tractable subgroup `S` with a rational trigonal map, the trigonal
construction gives a genus-3 curve `X` whose `F_p`-isomorphism class depends
only on `(H, S)`. Proposition 4 decides whether the isogeny `J_H → J_X` is
rational; when it is not, `J_X` is the quadratic twist of `A = J_H/S`
(Remark 4). So the Proposition 4 bit of `(H, S)` says whether the
principally polarized threefold `A = J_H/S`, which lies in the `L`-class of
`H`, is a Jacobian over `F_p` (bit 1) or the quadratic twist of a Jacobian
(bit 0). Remark 4 also gives `s' = c^3 s` for the twist `H'` of `H`, so the
bits of `(H, S)` and `(H', S')` are complementary.

## Hypothesis under test (H_comp)

The bit is a coin whose bias is a property of the `L`-class: a random target
`A` in class `L` that is geometrically a plane-quartic Jacobian is an actual
Jacobian with probability

`q(L) = m(L) / (m(L) + m(L^tw))`,

where `m(L)` is the mass (`Σ 1/#Aut`) of smooth plane quartics over `F_p`
with `L`-polynomial `L`, and `L^tw` is the quadratic-twist class
(`(a1, a2, a3) → (−a1, a2, −a3)`). Classes in which plane quartics are
scarce relative to the twisted class then hold few successes, which would
produce the observed clustering. `q(L)` is estimated independently of every
hyperelliptic datum, from a Monte Carlo census of plane quartics.

## Frozen inputs

- `research/smith_iso_granularity_20261015/data/iso_p{7,11}.csv.gz`
  (isomorphism-class rows; `full`, `size`, `L`-key, stratum), read-only.
- Convention for `L`-keys as in that study: power sums `P_k = p^k + 1 −
  #C(F_{p^k})`, `e1 = P_1`, `e2 = (e1 P_1 − P_2)/2`,
  `e3 = (e2 P_1 − e1 P_2 + P_3)/3`, key `(−e1, e2, −e3)`.
- Seeds: quartic sampling `SplitMix64(20261016 · 1000 + p)`; simulations
  `SplitMix64(20261016)`.
- Quartic samples: 4,000,000 at `p = 7`, 2,000,000 at `p = 11`.

## Measurements

1. **Subgroup census (every tractable isomorphism class).** On the
   representative `y^2 = t f(x)`, run `smith_prop24::analyze_curve` and, for
   each tractable subgroup, record whether it is *constructible* (Prop. 2
   holds and a valid trigonal map was built) and its Prop. 4 bit on this
   model. Subgroups where Prop. 2 holds but no map is built carry no bit in
   either twist (they cannot give `full_success`) and are counted
   separately. Per isomorphism class: `n_c` constructible subgroups, `n_1`
   with bit 1. The twin of a class is the class of `y^2 = c t f(x)`
   (looked up through its orbit).
   `ρ(L) = Σ n_1 / Σ n_c` over the isomorphism classes of `L` (each class
   counted once).
2. **Quartic census.** Uniform random ternary quartic forms over `F_p`
   (15 coefficients, not all zero). Keep the nonsingular ones (test below),
   compute `#C(F_{p^k})` for `k = 1, 2, 3` and the `L`-key. `n_Q(L)` =
   number of smooth samples with key `L`; this estimates `m(L)` up to a
   common factor because a uniform form lands in the `PGL_3`-orbit of `C`
   with probability proportional to `1/#Aut(C)`.
   Nonsingularity: after a random linear change of coordinates, the
   resultants `Res_y(Q_x, Q_y)` and `Res_y(Q_x, Q_z)` (chart `z = 1`) have
   no common root and `Q_x, Q_y, Q_z` have no common zero on `z = 0`;
   three independent coordinate changes are tried before a form is declared
   singular (a singular point survives every projection).
3. **Nulls (999 simulations each, per prime).** Keep every isomorphism
   class, its size, `L`-key and `n_c`. For each twin pair `{H, H'}`
   (`H ≠ H'`) draw `n_c` bits `b_i ~ Bernoulli(q)`; `H` is a success iff
   some `b_i = 1`, `H'` iff some `b_i = 0`. A self-twin class gets
   independent fair bits. Recompute `C` (share of models in `L`-classes with
   a success).
   - **N_q:** `q = q̂(L_H)` drawn per simulation and class pair from
     `Beta(1 + n_Q(L_H), 1 + n_Q(L_H^tw))` (sampling error of the census
     included).
   - **N_0:** `q = ρ̄`, the pooled bit rate over all constructible subgroups
     (no class effect).

## Controls (all must pass, else the stage is invalid)

- **Q1:** `[n_1 > 0]` equals the frozen `full` of every isomorphism class,
  and the observed `C` equals the frozen `C` (0.4304 / 0.4750 to all
  printed digits).
- **Q2 (Remark 4, exact):** for every twin pair, `n_c(H) = n_c(H')` and
  `n_1(H) + n_1(H') = n_c(H)`; violations are counted and must be 0.
- **Q3 (census sanity):** for 500 smooth samples per prime, the `L`-key
  predicts `#C(F_{p^4})` and a direct count agrees; every sampled key
  satisfies the Weil bound `|P_1| <= 6√p`.
- **Q4 (census symmetry):** recorded only — `n_Q` summed over all classes
  equals the smooth-sample count.

## Decision rule (per prime)

`C` = observed class share; gap share `g = (med N_0 − med N_q)/(med N_0 − C)`.

- `C` inside the central 95% of `N_q` → **"composition explains"**.
- `C` below the 2.5% quantile and `g >= 0.5` → **"partial"**.
- `C` below and `g < 0.5` → **"does not explain"**.
- `C` above the 97.5% quantile → **"over-explains"**.

Overall: the common verdict if both primes agree, else "prime-dependent".

Secondary (association, not a verdict): weighted Pearson correlation of
`ρ(L)` with `q̂(L) = (1 + n_Q(L)) / (2 + n_Q(L) + n_Q(L^tw))` over
`L`-classes with `Σ n_c >= 1`, weights `Σ n_c`, with a 999-shuffle
permutation p-value (shuffling `q̂` across classes).

## Stage 2 trigger

If the verdict is not "composition explains" at both primes, stage 2
proceeds: explicit isomorphism classes of `X` per `(H, S)` (Smith §6 model,
canonical plane quartic, `PGL_3` canonical form) to measure target sharing,
under its own frozen protocol.

## Stop conditions and cost

2 threads. If the quartic census at a prime exceeds 6 hours wall clock, the
samples completed so far are used and the shortfall recorded. Exact
enumeration of hyperelliptic classes as in PR #1928; degree-8 models only.

## Scope and limits

`q(L)` assumes targets are spread over the class's plane-quartic Jacobians
and their twists in proportion to mass; hyperelliptic targets (both twists
are Jacobians) are ignored. "Composition explains" means the class-level coin
predicted by an independent quartic census reproduces `C`; it does not
identify individual targets (stage 2) and says nothing at cryptographic
sizes.
