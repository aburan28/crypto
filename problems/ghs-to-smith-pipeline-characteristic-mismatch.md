# GHS (characteristic 2) cannot feed Smith's genus-3 trigonal construction

- **Status:** `open`. The characteristic-2 pipeline stays blocked. The
  odd-characteristic alternative exists; Smith's Props. 2 and 4 on it are
  unchecked.
- **Kind:**
  - proved obstruction for the pipeline as stated (characteristic 2),
    unchanged;
  - for the odd-characteristic alternative, no characteristic obstruction
    remains. It is an implementation gap: Smith's Props. 2 and 4 and the
    trigonal construction are not implemented.
- **Affected work:**
  - [weak isogeny classes study, 2026-10-13](../research/weak_isogeny_classes_20261013/REPORT.md),
    section "Where the leads meet";
  - [odd-characteristic descent study, 2026-10-14](../research/odd_char_descent_genus3_20261014/REPORT.md).
- **Next action owner:** unassigned
- **Last updated:** 2026-10-14

## Problem

The proposed end-to-end pipeline does not compose. It reads: walk
`E/F_{2^N}` to a curve whose GHS descent is a genus-3 hyperelliptic curve,
apply Smith's isogeny to a plane quartic, run Diem's attack. The GHS curve
`C` lives over `F_{2^l}`, and Smith's construction is stated and proved only
away from characteristics 2 and 3.

Since 2026-10-14 the odd-characteristic replacement for the first step is
known to exist (see the evidence below). What remains open is whether
Smith's isogeny is `F_q`-rational, and non-degenerate, on its output.

## Evidence and uncertainty

- **Observed (source):** Smith, *Isogenies and the discrete logarithm problem
  in Jacobians of genus 3 hyperelliptic curves*, J. Cryptology 22 (2009),
  arXiv:0806.2995. The introduction says "the trigonal construction (and hence
  our formulae) does not apply in characteristics 2 and 3". All its
  algorithms assume characteristic at least 5.
- **Derived (labeled argument, not machine-checked):** for an ordinary
  genus-3 Jacobian `J` over a field of characteristic 2, `J[2](k̄)` has order
  `2^3`, the 2-rank. It is not `2^6`, so Smith's `(Z/2)^3` kernels inside a
  6-dimensional `J[2]` of Weierstrass-point differences do not exist. The
  étale part of `J[2]` is the kernel of Verschiebung, and its quotient is a
  Frobenius twist of `J`, not a plane-quartic Jacobian. GHS curves are
  hyperelliptic of the form `y^2 + h(x) y = f(x)`, so their Weierstrass
  structure differs from that of `y^2 = f(x)`.
- **Observed (sources, 2026-10-14):** an odd-characteristic descent onto
  genus-3 hyperelliptic curves exists for `n = 3`:
  - Momose–Chao (ECC 2007, slide 24): `E/F_{q^3}: y^2 = e g(x)(x−α)(x−α^q)`,
    `α ∈ F_{q^3} \ F_q`, `g ∈ F_q[x]` of degree 1 or 2. The covering has
    degree 8 and there are `Θ(q^2)` classes.
  - Joux–Vitse (Eurocrypt 2012, §4.1 eq. (3), App. A.1) give the explicit
    `H` over `F_q` and call these "the only" such curves (citing Diem 2003,
    Momose–Chao 2005, Thériault 2003).
  - This holds for any odd `q`, so in particular for `p >= 5`.
  - Diem 2003 itself was **not read**: no copy is reachable from the
    session host. Its content is known here only through the papers that
    cite it.
- **Observed (repository, before this update, not previously linked):**
  - [`research/ghs_cubic_cover_census_20261012`](../research/ghs_cubic_cover_census_20261012/README.md)
    (PR #1759) and
    [`research/ghs_cubic_hyperelliptic_genus3_20261013`](../research/ghs_cubic_hyperelliptic_genus3_20261013/README.md)
    count and build this `m = 3, r = 5` type over `F_q`, `q ∈ {5,7,11,13}`.
  - They check `L(T) = 1 − tT^3 + q^3T^6`.
  - They execute the transfer `E(F_{q^3}) → Jac(C)(F_q)` with Cantor
    arithmetic.
- **Exact computation (2026-10-14, [study](../research/odd_char_descent_genus3_20261014/REPORT.md)):**
  - Joux–Vitse's model passes the L-polynomial check on 81 of 81 instances
    (`p ≤ 31`).
  - The Frobenius cycle type on the eight Weierstrass points of `H` is
    `(3,3,1,1)` for split or linear `h`, giving `#S(H) = 3`
    `F_q`-rational tractable subgroups (Smith, Lemma 1). It is `(6,2)` for
    irreducible `h`, giving `#S(H) = 1`.
  - This is the "third arm" screen. It used the study's own
    distinct-degree factorisation, equivalent to
    `smith_kernel_coverage::cycle_type`, and was checked against Smith's
    Lemma 1 table.
- **Derived (with proof in the study, not machine-checked):**
  - The GHS genus is `1 + 2^{m−2}(r−4)` for every degree-2 `x`,
    independent of `q`.
  - Genus 3 needs `(m, r) ∈ {(2,6), (3,5)}` and a branch-set period of 2
    or 3. So for `F_q(j(E)) = F_{q^n}` with `n >= 4` no odd-characteristic
    GHS cover has genus 3.
- **Missing capability:** a native test of Smith's Prop. 2 (rational
  trigonal map) and Prop. 4 (leading coefficient of `s(t)` a square) on this
  family. Also an executed trigonal construction showing that `J_H/S` is a
  smooth plane-quartic Jacobian. The family has three bielliptic
  involutions, so Smith's random-curve hypotheses (1/2 each) do not
  obviously apply. See also
  [`smith-kernel-isogeny-class-walk.md`](smith-kernel-isogeny-class-walk.md),
  which needs the same Prop. 2/4 tests.
- **What this does not establish:**
  - It shows nothing about characteristic 2 or 3.
  - It gives no cost, and it does not show that the pipeline beats rho.
  - It does not show that a plane quartic results.
  - The non-hyperelliptic genus-3 covers (Momose–Chao Types 1 and 2; the
    repository's `m = 2` build) already give plane quartics without Smith.
    Smith matters only for curves whose sole genus-3 cover is hyperelliptic
    (20 at `q = 5` and 42 at `q = 7` in the census).

## What would unblock it

- [ ] A cited construction of a rational isogeny from a genus-3 hyperelliptic
  Jacobian in characteristic 2 to a plane-quartic Jacobian, with a native
  implementation checked on toy curves (point counts equal on both sides).
- [x] Or an odd-characteristic descent producing genus-3 hyperelliptic
  curves (sources above; built in the repository; re-checked 81/81),
  - [x] screened by Frobenius cycle type (`#S(H) = 3` or `1`),
  - [ ] and by Smith's Props. 2 and 4, with the trigonal construction
    executed on toy instances and the quartic's `L`-polynomial equal to
    `L_H`.

Either completed route gives an executable pipeline. Neither shows that the
pipeline beats rho on the source curve.

## Next discriminating step

- **Action:** implement Smith's Props. 2 and 4 natively (shared with
  `smith-kernel-isogeny-class-walk.md`). Apply them to the three tractable
  subgroups of each Joux–Vitse `H` in
  `research/odd_char_descent_genus3_20261014/results/run-20261014T1800Z/jv_family.csv`.
- **Controls:** random genus-3 hyperelliptic curves of the same cycle type
  `(3,3,1,1)`, to see whether the family's success rate differs from
  Smith's 1/4 per subgroup.
- **Outcomes:**
  - If some subgroup passes and the quotient is a smooth quartic with a
    matching `L`, add the family as the third arm of the coverage study.
  - If every subgroup fails on every instance, record that as a bounded
    negative for this family. Types 1 and 2 remain the direct route to
    plane quartics.
- **Recheck trigger:** a cited char-2 or char-3 trigonal analogue; access to
  Diem 2003.

## Resolution and history

- 2026-10-13: opened from the weak isogeny classes study. Unresolved.
- 2026-10-14: proposal 5 run ([study](../research/odd_char_descent_genus3_20261014/REPORT.md)).
  An odd-characteristic genus-3 hyperelliptic family exists (`n = 3`), so
  the record is **not** closed as a proved obstruction. The cycle-type
  screen is done; Smith's Props. 2 and 4 remain. The characteristic-2
  statement is unchanged.
