# Protocol: the hyperelliptic genus-3 GHS cover (`m = 3, r = 5`) over `F_q`, its Frobenius, and the transfer `E0(F_{q³}) → Jac(C)(F_q)`

Frozen 2026-10-13, before `hyper_genus3.rs` was run. Conductor T-2463.
**Structural stage diagnostic on toy fields.** It builds a curve, checks its
Frobenius and executes the transfer as a group map with Cantor arithmetic. It
runs no index calculus, collects no relations, solves no discrete logarithm and
claims no cost or speedup. Follows `research/ghs_cubic_cover_census_20261012/`
(PR #1759), which counted this type without building it, and the plane-quartic
build `research/ghs_cubic_genus3_build_20261012/` (PR #1761), whose `m = 2`
curves have no class arithmetic in the repository
(`problems/plane-quartic-jacobian-class-arithmetic-missing.md`).

## Construction under test (derived here; verified by the run)

`q` an odd prime `> 3`, `K = F_{q³}`, `L = F_{q⁶}` (all arithmetic is done in
`L`; `F_q, F_{q²}, K` are its subfields), `σ` the `q`-power map. Choose

- `p ∈ K \ F_q`, `r0 = p`, `r1 = p^q`, `r2 = p^{q²}` (distinct);
- a monic squarefree `Q = x² + c1 x + c0 ∈ F_q[x]`: **split arm** `Q = (x−a)(x−b)`,
  `a < b` in `F_q`; **conjugate arm** `Q` irreducible over `F_q`.

`F0 = Q·(x−r0)(x−r1)`, `F1 = F0^σ = Q·(x−r1)(x−r2)`, `F2 = Q·(x−r2)(x−r0)`,
`N = (x−r0)(x−r1)(x−r2) ∈ F_q[x]`. `E0 : y0² = F0(x)` over `K`. In the census
language this is `τ = σ`, `B = {a, b, r0, r1}`, `τB = {a, b, r1, r2}`,
`τ²B = {a, b, r2, r0}`: `r = |B ∪ τB ∪ τ²B| = 5`, `m = 3`, genus
`2^{m−2}(r−4)+1 = 3`, and the rank-2 subspace `{B∆τB, B∆τ²B, τB∆τ²B}` has
union `{r0, r1, r2}` of size 3 (hyperelliptic test of `cover_census.rs`).

- **`C`**: function field `K(x, y0, y1, y2)`, `y_i² = F_i(x)`, a `(Z/2)³` cover of
  `P¹`. **Descent datum** `σ~ : y0 ↦ y1 ↦ y2 ↦ y0`, `x ↦ x`, `σ` on `K`;
  `σ~³ = 1`. Frobenius of the `F_q`-model on points:
  `Fr(x, y0, y1, y2) = (x^q, y2^q, y0^q, y1^q)`.
- **Genus-0 quotient `D`** (fixed field of `ι : y_i ↦ −y_i`): `u01 = y0y1/(Q(x−r1))`,
  `u12 = y1y2/(Q(x−r2))`, `u20 = y2y0/(Q(x−r0))` with `u01² = (x−r0)(x−r2)`,
  `u12² = (x−r1)(x−r0)`, `u20² = (x−r2)(x−r1)`, `u01u12u20 = N`. The
  anticanonical coordinates `X0 = u12u20, X1 = u20u01, X2 = u01u12` satisfy
  `X0² = N(x−r1)`, `X1² = N(x−r2)`, `X2² = N(x−r0)`, hence the conic
  `(r2−r0)X0² + (r0−r1)X1² + (r1−r2)X2² = 0`, and
  `x = (r0X0² − r1X2²)/(X0² − X2²)`. `σ~` cycles `X0 → X1 → X2`, so in
  `Z_j = Σ_i θ_j^{q^i} X_i` (`θ_j = γ^j`, `γ` a fixed generator of `K*`) the conic
  has `F_q` coefficients. A conic over a finite field has a rational point.
- **Hyperelliptic model.** Parametrize the `F_q`-conic by lines through an
  `F_q`-point `P*`, `z ∈ P¹(F_q)`; write `x(z) = xn(z)/xd(z)` reduced, `xd` monic.
  `w := y0y1y2·xd(z)/(Q(x)N(x))` is `σ~`-invariant, `ι`-anti-invariant, and
  `w² = h(z) := xn² + c1·xn·xd + c0·xd²` (because `F0F1F2 = Q³N²`). So
  `C : w² = h(z)` over `F_q`, with the explicit map `(x, y0, y1, y2) ↦ (z, w)`.
  Branch points of `C → D` are the 8 points of `D` over `x ∈ {a, b}`.
  - Split arm: exactly one `F_q`-point of `D` lies over `x = a` (derived from
    the norm condition `N_{K/F_q}(u01) = N(a)`); take it as `P*` with the
    tangent direction at `z = ∞`, so `h` has degree 7 (imaginary model, the
    `F_q`-rational Weierstrass point `P_a` at infinity).
  - Conjugate arm: no `F_q`-point of `D` lies over `a, b`; `P*` is the first
    `F_q`-point of the conic in enumeration order and `h` has degree 8.
- **Transfer** (split arm only; Cantor arithmetic is implemented for the
  imaginary model only). Base point of `E0` is `O = (r0, 0)` (group law through
  the Weierstrass model obtained by sending `r0` to infinity). Define
  `φ(P) = Σ_{j=0}^{2} Fr^j(π0^*[P] − 4∞)` in `Jac(C)(F_q)`, where `π0^*[P]` is the
  four points `(x_P, y0_P, ±√F1(x_P), ±√F2(x_P))` over `L`. This equals
  `(1 + Fr + Fr²)∘π0^*([P] − [O])` because `π0^*[O]` is the four points over
  `r0`, i.e. two `ι`-pairs, linearly equivalent to `4∞`.
- **Derived kernel bound.** `π0_*π0^* = [4]` (`deg π0 = 4`) and
  `π0_*∘Fr^j∘π0^* = π0_*π_j^*∘Frob^j`, which for `j = 1, 2` sends `[P] − [O]` to
  twice the difference of two `x`-fibres of `E0`, i.e. to `0`. Hence
  `π0_*∘φ = [4]` on `E0(K)` (the `m = 2` build had `[2]` because there
  `deg π0 = 2`) and `ker φ ⊆ E0(K)[4]`. Further, `K(E0)(u12)` is an unramified
  double cover of `E0` inside `C` (`u12² = (x−r1)(x−r0)` has even valuation at
  every point of `E0`), so `π0^*` kills the 2-torsion class
  `T1 = [(r1,0)] − [O]`: `T1 ∈ ker φ`, and `φ` is not surjective.

## Instances (frozen inputs)

Field `L = F_{q⁶}` from the first primitive polynomial found by the tool's
deterministic search; `K` is printed in the census basis
(`cover_census.rs`'s `F_{q³}`, transported by an explicit isomorphism found by
root search, first root in log order). Order: `p` over `K \ F_q` in increasing
census index `c0 + c1·q + c2·q²`, keeping only the minimum-index element of each
`σ`-orbit; for each `p`, `Q` in lexicographic order (`(a, b)`, `a < b`, for the
split arm; `(c1, c0)` for irreducible `Q`). Select the first cases with
`#E0(K) = h·ℓ`, `h ∈ {4, 8}`, `ℓ ≥ 5` prime:

| q | split arm | conjugate arm |
| --- | --- | --- |
| 5 | first 3 | first 2 |
| 7 | first 4 (the first is the explicit example) | first 2 |
| 11 | first 3 | first 2 |
| 13 | first 2 | first 2 |

Every candidate `(p, Q)` met before the last selected one is recorded with its
`#E0(K)`. Seed for random pairs and points: `20261013`, fixed.

## Pre-registered predictions and falsifiers

- **P1 (descent).** The conic matrix in `Z`-coordinates, `x(z)` (after
  reduction, `xd` monic) and `h(z)` have all coefficients in `F_q`; the map
  `(x, y0, y1, y2) ↦ (z, w)` commutes with Frobenius on random points of `C`
  over `L` (`map(Fr P) = (z^q, w^q)`); `x(z(P)) = x` and `w² = h(z)` there.
  The split arm has exactly one `F_q`-point of `D` over `a`.
  *Falsifier:* any coefficient outside `F_q`, any Frobenius or round-trip mismatch.
- **P2 (genus 3, hyperelliptic, smooth).** `deg x(z) = 4`; `h` squarefree
  (`gcd(h, h') = 1`) of degree 7 (split) or 8 (conjugate).
  *Falsifier:* other degree or a repeated root.
- **P3 (Frobenius).** Point counts on `w² = h(z)` give `N1 = q+1`, `N2 = q²+1`,
  `N3 = q³+1−3t_E`, `t_E = q³+1−#E0(K)`, i.e. `L(T) = 1 − t_E T³ + q³T⁶` and
  `#Jac(C)(F_q) = L(1) = #E0(K)`. Cross-checks: (a) the same `N_k` counted on
  the normalized `K`-model fibre by fibre (formulas in the README, special
  fibres over `a, b, r_i, ∞` derived by hand); (b) split arm: exhaustive
  enumeration of reduced Mumford pairs `(u, v)` over `F_q` gives
  `#Jac(C)(F_q) = #E0(K)`.
  *Falsifier:* any `N_k`, `L(1)` or enumeration mismatch.
- **P4 (transfer, split arm).**
  (a) every `φ(P)` is a reduced Mumford pair with `F_q` coefficients;
  (b) `φ(P ⊕ Q) = φ(P) + φ(Q)` on 400 random pairs per instance (seeded) with
  `P, Q, P ⊕ Q` eligible, plus 100 doublings;
  (c) for a generator `G` of `E0(K)[ℓ]`: `φ(G) ≠ 0`, `ℓ·φ(G) = 0`;
  (d) `π0_*(φ(P)) = 4P` for every eligible `P` whose reduced divisor's points
  avoid the special fibres (counted);
  (e) `ker φ ⊆ E0(K)[4]` and `T1 ∈ ker φ`, measured over all of `E0(K)`.
  *Eligible set:* `P` with `x_P ∉ {a, b, r0, r1, r2, ∞}`. Exceptional points are
  counted exactly; `O ↦ 0`; `(a,0) ↦ [u_a, 0]` and `(b,0) ↦ [u_b, 0]`
  (`u_a, u_b` the monic `z`-polynomials of the Weierstrass points over `a, b`)
  and `(r1, 0) ↦ 0` are derived directly; every exceptional `P` (these and the
  `x = r2, ∞` points) is also evaluated as `φ(P ⊕ R) − φ(R)` for three eligible
  `R` (a linear equivalence on `E0`, not an assumption about `φ`), which must
  agree with each other and with the direct values. Every point of `E0(K)` is
  thus accounted for.
  *Falsifier:* any failure of (a)–(e) or disagreement among routes.
- **P5 (census cross-check).** (a) For `q = 5, 7` (exhaustive census) and
  `q = 11` (sample, where present), the census record of the `K`-isomorphism
  class of the short Weierstrass model of `E0` (and of its two `σ`-conjugates)
  has `min_genus = 3`, includes the type `m = 3, r = 5, hyperelliptic`, and its
  `order` equals `#E0(K)`. (b) Census corollary: genus 3 needs `4 | #E(K)`, so
  every enumerated candidate in **both** arms has `4 | #E0(K)` (exhaustively
  over all `(p, Q)` at every `q` run). *Falsifier:* a record without the type,
  an order mismatch, or a candidate with `#E0(K) ≡ 2 mod 4` (a conjugate-arm
  `#E0(K)` odd is impossible: `Q | F0` gives a `K`-rational 2-torsion point).

## Not measured / out of scope

Cantor arithmetic for the degree-8 (conjugate-arm) model, so no transfer test
there. No index calculus, relation collection, DLP solver or cost. No claim
about deployed curves.

## Environment

`rustc --edition 2021 -O hyper_genus3.rs`, dependency-free, ≤ 6 threads (the
tool is single-threaded), one compile at a time on the shared host. Output
JSONL + `env.txt` (command, base commit, source and protocol sha256,
toolchain, host, uptime). Python only plots recorded values.

## Amendment 2026-10-13 (after development runs, before the recorded run)

Written after unrecorded development runs of the tool at q = 5, 7, 11, 13
(scratch output, not kept; the recorded run re-runs everything from scratch
with the same deterministic code plus these two changes). Nothing above is
changed.

1. **Conjugate arm at q = 7 has no admissible instance.** All 2352
   conjugate-arm candidates at q = 7 have `#E0(K) ∈ {320, 336, 352, 368}`, all
   divisible by 16, so the rule `h ∈ {4, 8}` selects nothing. Rule added: an
   arm with no admissible candidate takes its first 2 candidates in the frozen
   order, with `ℓ` the largest prime factor and `h = #E0(K)/ℓ`, flagged
   `selected_by_amendment`; P1–P3 and P5 are tested on them exactly as above
   (the cofactor plays no role there; no transfer is run on this arm anyway).
2. **P4(d) extension.** The development run showed that many reduced
   divisors `φ(P)` contain points over `x = ∞` (and some over `a, b, r_i`), so
   the preregistered "avoids the special fibres" set is about half of the
   eligible points at q = 5. Added, reported separately and also required by
   the P4 flag: `π0_*(φ(P)) = 4P` for **every** `P ∈ E0(K)`, using derived
   images of special-fibre points: over `x = ∞`, `v0 = w/(xn(z)·ρ1ρ2)` with
   `(ρ1, ρ2) = (X0/X1, X2/X1)` at `z` and image `(S, Y) = (0, e1·v0)`; over
   `a, b, r1`: `(x, 0)`; over `r0`: `O`; over `r2`:
   `u12 = (r2−r1)·X2/X0`, `y0 = w·(r2−r0)(r2−r1)/(u12·xd(z))`. The
   preregistered metric is still reported unchanged.
