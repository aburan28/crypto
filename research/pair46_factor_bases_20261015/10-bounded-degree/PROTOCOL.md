# 10-bounded-degree — preregistered protocol (frozen before any data)

Item 10 of the pair-46 factor-base campaign: bounded-degree algebraic factor
bases over prime fields. Source proposal: `IDEA-20261009-a7733f`
(crypto-autoresearcher `ledger/proposals/IDEA-20261009-a7733f.yaml`), design
audit `coordination/goals/GOAL-ECDLP-001/proposals/ALGEBRAIC-FACTORBASE-DIMENSION-20261010/design-audit.md`.

## Step 0: repair of the design-audit errors (before any data)

1. **Field table (audit "Immediate protocol error").** The proposal's
   `p in {2^16, 2^20, 2^24}` are powers of two, not primes. Replaced by:
   `p_k` = the smallest prime `p > 2^k` with `p = 3 (mod 4)`, for
   `k in {16, 20, 24}`. Primality is certified by deterministic Miller-Rabin
   with the 12 bases 2..37 (deterministic below 3.3e24) and recorded with its
   exact bit length. `p = 3 mod 4` is chosen so that square roots are one
   exponentiation `z^((p+1)/4)`.
2. **Extension control (audit "positive comparison").** `q_k` = smallest prime
   `q > 2^k`, `q = 3 mod 4`, `k in {8, 10, 12}`, so `Q = q^2` spans the same
   sizes 2^16..2^24. `F_{q^2} = F_q[i]/(i^2+1)`; `i^2+1` is irreducible because
   `(-1)^((q-1)/2) = -1` for `q = 3 mod 4` (checked and recorded). Factor base
   `F_ctl = {(x,y) in E(F_{q^2}) affine : x in F_q}` (x has zero `i`-component).
3. **Factor base defined as a subset of E(F_p) (audit "scope").** For every
   family below, the factor base is a SET of affine points of `E(F_p)`
   (counted once each: the union, never per component or with multiplicity);
   the point at infinity is excluded.
4. **Hasse-Weil hypotheses / p^0.4 claim (audit).** The proposal's "genus
   about p^0.4" exponent is not used. The exact bound tested here (claim class
   X2) is stated with its hypotheses: `C_f` is the curve `y^2 = h_f(t)`
   below, geometrically irreducible (h_f not a constant times a square) with
   squarefree `h_f` (smooth affine model), genus `g_f` computed exactly from
   the odd-multiplicity part of `h_f`.
5. **Representation size model (audit).** Dense coefficient encoding:
   `size_bits = (#coefficients) * bitlen(p)` (`d+1` for a polynomial,
   `2d+1` for a rational map with monic denominator). Membership cost is
   measured on the 53-bit pair-46 crater curve as counted F_p multiplications
   per membership test (gcd(t^p - t, f(t) - x) for polynomials).
6. **Relation and linear-algebra cost (audit).** Not measured: see "Cut
   scope". The exact size lower bound alone already exceeds `sqrt(p)`; that
   implication (X3) is stated as exact, not measured.

## Claim classes

- **X1 (exact).** For any map `f: P^1 -> P^1` of degree `d` over F_p,
  `|f(P^1(F_p))| >= (p+1)/d` (fibers have <= d rational points).
- **X2 (exact, with hypotheses of Step 0.4).** Let `F_f = {P in E(F_p) affine :
  x(P) in f(F_p)}`. For polynomial `f` (resp. rational `f = P/Q`, deg P = d,
  Q monic of degree d-1), `C_f: y^2 = h_f(t)` with `h_f = g(f(t))` (resp.
  `h_f = Q (P^3 + a P Q^2 + b Q^3)`), `g = x^3+ax+b`. `(t,y) -> (f(t), y/Q^2)`
  maps affine points of C_f with Q(t) != 0 onto F_f with fibers <= d, so
  `|F_f| >= (p + 1 - 2 g_f sqrt(p) - 2 - 2 deg Q - 1) / d =: LB_f`. Every row
  is checked against `LB_f` (an implementation check, not evidence).
- **X3 (exact, derived).** With `g_f <= 2d - 1`, `|F_f| <= sqrt(p)` forces
  `d >= sqrt(p)/5 (1 - o(1))`, i.e. dense size `>= sqrt(p) log p / 5` bits;
  `|F_f| < p^0.9` forces `d >~ p^0.1 / 5`. More generally for any cover
  `C -> E` of degree d and genus g: `|F| >= (p+1-2g sqrt p)/d`, so a sub-p^0.9
  factor base needs `d >~ p^0.1` or `g >~ sqrt(p)/2` (corrects the proposal's p^0.4).
- **X4 (exact).** A zero-dimensional condition `G(x,y) = A(x) + y B(x) = 0`,
  `deg A <= d`, `deg B <= d - 2`, has pole order `<= 2d` at O, so `|F_G| <= 2d`.
- **E1 (empirical, toy).** For random polynomials / rational maps of degree
  `d`, `rho = |F_f|/p` is close to the Chebotarev-generic value
  `c_d = sum_{k=1..d} (-1)^(k+1)/k!` (image density for generic monodromy S_d)
  and never approaches the extension-field control exponent 1/2.
- **E2 (empirical, pair 46).** On the 53-bit crater curve, Monte Carlo
  estimate of `rho` for random polynomials d = 2..8 and the counted membership
  cost.

## Metric (same instrument on every arm)

Exact enumeration: for each x in F_p the image indicator `I_f(x)` is built by
evaluating f on every t in F_p (rational maps: batch inversion; poles go to
infinity, excluded); `|F| = sum_x I_f(x) (1 + chi(g(x)))` with
`chi` = Legendre symbol, computed once per curve. Reported: `|image|`, `|F|`,
`rho = |F|/p`, exponent `e = log|F| / log(field size)`, `g_f`, `LB_f`,
`size_bits`. The extension control uses the same formula with `chi_{q^2}(z) =
Legendre_q(Norm(z))` and field size `Q = q^2`.

Effect size: for each d and random family, mean of rho over 3 primes x 3 seeds
(n = 9) with a 95% t-interval, and `Delta = rho - c_d` with interval;
headline exponent gap `min_e(algebraic rows) - max_e(control rows)`.

## Arms

For each prime `p_k` one curve `y^2 = x^3 + a x + b`, `a, b` from splitmix64
seeded `1825564914 + 1000 k`, rejecting `4a^3 + 27b^2 = 0`. `d in 2..8`.

- **RP** random monic polynomial of degree d, seeds 1825564914, ...915, ...916.
- **RR** random rational map `P/Q`, P degree d (random leading coeff != 0),
  Q monic degree d-1, rejected unless gcd(P,Q) = 1; same 3 seeds.
- **PW** power map `t^d` (deterministic; extremal small-image family).
- **DK** Dickson `D_d(t,1)` (deterministic).
- **ZD** zero-dimensional `A + yB` (X4), same 3 seeds.
- **NULL / negative controls.** (i) `f(t) = t` (d = 1): must give
  `|F| = #E(F_p) - 1` exactly. (ii) RF: a random function F_p -> F_p
  (splitmix64 hash; non-algebraic, "degree" undefined), 3 seeds: rho ~ 1 - 1/e.
- **POSITIVE controls.** (i) Extension control `x in F_q` over `F_{q^2}` for
  the three `q_k`: valid iff `0.40 <= log_Q|F_ctl| <= 0.60` for every q.
  (ii) Enumerator check: PW image must equal `(p-1)/gcd(d,p-1) + 1` exactly;
  DK must be a permutation (image = p) when `gcd(d, p^2-1) = 1`.
  (iii) Point counts: `#E` by the Legendre sum, checked by `[#E] R = O` for
  4 random points R and Hasse `|#E - p - 1| <= 2 sqrt p`.
  Any failed control marks the affected rows `invalid` (never evidence).
- **P46** pair-46 crater `p = 5947698030651343`,
  `a = 5757116383569871, b = 3560173328417442`: RP with seed 1825564914,
  d = 2..8, 20000 uniform x per d (seed 1825564915 + d); estimate
  `rho_hat = mean(I_f(x)(1 + chi(g(x))))` with 95% normal interval; count
  F_p multiplications per membership test.

## Stop rule

Run every row above once (no reruns, no added arms). Wall-clock budget 35
minutes on 1 core. If exceeded, stop and report completed rows; cut order:
RR at p_24 first, then P46 samples down to 5000 per d. No row is dropped for
its value.

## Cut scope (preregistered)

Relation generation (Semaev summation) and linear algebra are not run: by X2
every bounded-degree family has `|F| >= p/d (1 - o(1)) >> sqrt(p)`, so any
complete index calculus over these bases has linear algebra in `|F|`
unknowns, exceeding rho's `sqrt(p)` before relations are generated. This
closes only the tested families (x-images of rational maps, zero-dim
conditions); it is not a universal theorem about all subsets of E(F_p).
