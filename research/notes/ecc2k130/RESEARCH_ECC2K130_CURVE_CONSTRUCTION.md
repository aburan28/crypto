# Producing a curve that carries the ECC2K-130 subgroup

**Code:** `src/cryptanalysis/curve_construction/` (native Rust), binary
`src/bin/ecc2k130_curve_construction.rs`
**Protocol, frozen before the run:** `research/notes/ecc2k130/curve_construction_20261006/PROTOCOL.md`
**Frozen report:** `research/notes/ecc2k130/curve_construction_20261006/results/curve_construction.json`
**Related:** `research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md` (the question
this note answers, and the claim it corrects), `research/notes/cm-isogeny/RESEARCH_MESTRE_HOWE.md`
(the genus-2 CM construction), `research/notes/ecc2k130/RESEARCH_KOBLITZ_INDEX_CALCULUS.md`
(route 4).

The question.  The ECC2K-130 subgroup `⟨G⟩` is `A(F_2)` for the simple
130-dimensional trace-zero part `A` of `Res_{F_2^131/F_2}(E_0)`, where
`E_a : y² + xy = x³ + a·x² + 1` (`a = 0` is the challenge family, trace −1; `a = 1`
has trace +1).  A curve `C/F_2` of genus between 130 and about 300 with
`P_A | P_C` (Frobenius characteristic polynomials) and a correspondence anyone
can evaluate would put the discrete logarithm below Pollard rho: index
calculus on a genus-130 curve with `Jac ~ A` is modelled at `2^37.17` against
rho's `2^60.81`.  **Can such a curve be produced?**  Four routes were tried.

**Bottom line.  No construction reaches the window, but the picture is not the
one the hyperelliptic note drew, and that note's strongest claim was wrong.**

- **The correction first.**  The earlier note said the genus any construction
  reaches over `F_2` is 1, `2^129` or `2^130`.  That holds for **GHS/Hess
  descent** and for nothing else.  The Klein quartic is a classical non-GHS
  counterexample at `n = 3`: GHS assigns it magic number 1, yet its Jacobian is
  `Res_{F_8/F_2}(E_1)` (point counts over `F_{2^k}`, `k ≤ 8`, all match).
- **Route 1, lifting to characteristic 0, produces an explicit curve.**  The
  reduction mod 2 of the modular curve `X_H(3²·7²·263²)` carries `A`
  (Eichler–Shimura gives exactly `A`'s Weil numbers).  Its genus is
  **508,799,809 ≈ 2^28.92**, a hundred bits of genus below GHS's `2^129` but still
  `2^20.7` above the window, and index calculus on it is astronomically above
  rho.  At `n = 3` the same construction *is* the Klein quartic, `X_H(49)`.
- **Route 2, cyclic covers, is closed below genus 1300.**  Every geometrically
  cyclic degree-131 cover of a curve of genus ≤ 1 over `F_2` whose Jacobian
  carries `A` has genus **≥ 1300**: class field theory forces at least `130/d`
  branch points, and `μ_d`-stability of the Prym's Weil numbers at least `2d`
  (`d` the order of Frobenius on the covering group).  The Kummer toys show the
  `μ_d`-stability directly.  Fermat quotients of exponent `7·263` carry no factor
  of level 263 or 1841 of the required CM type.
- **Route 3, direct search, succeeds at `n = 3` and fails at `n = 5`.**  `A_n` is
  never itself a Jacobian at `n = 3, 5`; the genus-2 and genus-4 searches are
  complete.  At `n = 3`, genus-3 curves carrying `A` exist for both signs,
  e.g. `y² + (x³ + x + 1)·y = x⁷ + x⁶ + 1` with `Jac ~ A_2(E_0) × E_1`.  Route 1
  only reaches genus 121 there.  At `n = 5` nothing exists at
  genus 4, and nothing hyperelliptic at genus 5.  At `n = 131` the search space
  is `2^395` hyperelliptic models.
- **Route 4, no curve, is the repository's existing Koblitz index-calculus
  line**: bounded above rho or closed on every route tried.

What stays open is narrower than before but not closed: a curve of genus
130…~300 carrying `A` that is neither a GHS descent nor a cyclic cover of a
curve of genus ≤ 1.  At `n = 3` exactly such curves exist.  At `n = 5` none
was found where the search is complete.

## 0. Boundaries, stated before measuring

From the protocol, and replayed natively (agreement with the legacy Python-era
file: 20 values, max |diff| `0.0048`, every exact field equal):

| quantity | value |
|---|---|
| rho reference on `⟨G⟩` (`⟨−1⟩ × ⟨π⟩`, 262 automorphisms) | `2^60.8090`, `S = 0.0774` |
| genus-130 curve with `Jac ~ A`, exact zeta-function model | `2^37.17` (`b = 17`, `|FB| = 2^14.01`) |
| window where index calculus over `F_2` beats rho | genus 130 … between 290 and 300 |
| previous best construction carrying `⟨G⟩` | GHS: `2^129` |

**Falsification target** (protocol): a curve `C/F_2` of genus ≤ 300 with
`P_A | P_C` and a polynomial-time correspondence — **not met**.  A construction
below `2^129` counts as a boundary refinement, not an attack — **met** by route 1.

## 1. The table

One unit: log₂ operations, the modelled cost of index calculus on the curve,
against the rho reference.  Costs are model figures — the exact zeta-function
model where the genus allows it, `L_{2^g}(1/2, √2)` beyond that (an
extrapolation; on the calibrated range it *overestimates* the exact model by
4.1 bits at genus 130 and 10.9 at genus 1300, see §8).

| route | construction | genus | carries `A` | log₂ cost | vs rho |
|:--:|---|---:|---|---:|---:|
| — | **target**: genus-130 curve with `Jac ~ A` (hypothetical) | 130 | by definition; existence unknown | 37.17 | **2^−23.64** |
| 0 | GHS/Hess, ECC2K-130 itself (magic 1) | 0 or 1 | **no** — the transfer is the zero map on `⟨G⟩` | — | — |
| 0 | GHS/Hess, any curve over `F_2^131` (magic 130/131) | `2^129`, `2^130` | not established | ≥ 129 (floor) | 2^+68.19 |
| 1 | **modular curve `X_H(30,503,529)` mod 2** | **508,799,809 = 2^28.92** | **yes** (Eichler–Shimura) | ≈ 169,981 (extrap.) | 2^+169,920 |
| 2 | cyclic degree-131 cover of any curve of genus ≤ 1 | ≥ 1300 | necessary bound | ≥ 148.92 | 2^+88.11 |
| 2 | — with `F_2`-rational group (Klein mechanism, `d = 1`) | ≥ 8451 | necessary bound | ≈ 459.9 (extrap.) | 2^+399.1 |
| 2 | — Kummer `y^131 = f` (`d = 130`) | ≥ 16,901 | necessary bound | ≈ 675.9 (extrap.) | 2^+615.1 |
| 2 | — with two branch points | 131 | **excluded**, every twist | — | — |
| 2 | Fermat quotient `y^1841 = x^a(1−x)^b` | 920 | **excluded** (Stickelberger) | 121.08 | 2^+60.27 |
| 3 | direct search at genus 130 | 130 | — | `2^395` hyperelliptic models | — |
| 4 | no curve: index calculus on `A` | — | — | bounded above rho or closed | — |

None of these is an advance in the sense of AGENTS.md §3: no row moves a ratio
to rho below 1.  The modular row moves the *constructible genus* by 100 bits;
the cyclic rows bound a family.  The correction in §7 is **accounting**.

## 2. The toy ladder — route 3, run to completion where it can be

Every smooth curve over `F_2` of genus 2, 3 and 4 was enumerated, with no
isomorphism reduction, so each curve appears many times:

- hyperelliptic models `y² + h y = f`, smooth on both charts by an exact `gcd`
  criterion;
- all `2^15 − 1` plane quartics;
- all canonical genus-4 curves `Q ∩ K` on the three quadric types over `F_2`.

Genus-5 hyperelliptic models were enumerated too.  `P_C` comes from point
counts; a target "hits" when it divides `P_C`.  Every hit is
smoothness-checked; a miss needs no check, because the enumeration is
complete.  Both toy sizes have 2 primitive mod `n`, as 131 does
(`ord_3(2) = 2`, `ord_5(2) = 4`, `ord_131(2) = 130`), so `t^n − 1` has the same
irreducible block.  This is evidence about the mechanism at that size and no
more.

| genus | family | models (smooth/valid) | `A_3(E_0)` | `A_3(E_1)` | `W_3(E_1)` | `A_5`, `W_5` (both signs) |
|---:|---|---:|---:|---:|---:|---:|
| 2 | hyperelliptic (all curves) | 768 | **0** | **0** | — | — |
| 3 | hyperelliptic | 6,144 | 96 (× `E_1`) | 96 (× `E_0`) | 0 | 0 |
| 3 | plane quartic | 32,767 | 56 (× `T²−2T+2`) | 224 | **56 = Klein** | 0 |
| 4 | hyperelliptic | 49,152 | 0 | 0 | — | **0** |
| 4 | canonical `(2,3)` | 196,605 | 80 | 304 | — | **0** |
| 5 | hyperelliptic only | 393,216 | 896 | 896 | — | 0 |

Four readings:

1. **`A_n` is not a Jacobian at `n = 3` or `n = 5`, for either sign.**  The
   genus-`(n − 1)` searches (genus 2, genus 4) are complete and have no hit.
   This is the toy form of "is `A` itself a Jacobian", and it says no twice.
2. **At `n = 3`, curves of genus `n` carrying `A` exist for both signs.**  For
   every hyperelliptic hit the extra elliptic factor is *the other sign's* curve;
   the plane-quartic hits carry other cofactors.  For our sign:
   `y² + (x³ + x + 1)·y = x⁷ + x⁶ + 1` has `Jac ~ A_2(E_0) × E_1`.
3. **At `n = 5`, nothing at genus 4 and nothing hyperelliptic at genus 5.**
   Non-hyperelliptic genus 5 and genus ≥ 6 are not enumerated.
4. **The Klein quartic is the only plane quartic over `F_2` with
   `Jac ~ Res_{F_8/F_2}(E_1)`**: the 56 hits form one `GL_3(F_2)` orbit with
   stabiliser of order 3.

## 3. Route 1 — the modular curve

**The `n = 3` case is classical, and it says where to look.**  The Klein quartic:

- has Jacobian `Res_{F_8/F_2}(E_1)` over `F_2`;
- over `Q`, its point counts match the Weil restriction from `Q(ζ_7)^+` of the
  conductor-49 curve `y² + xy = x³ − x² − 2x − 1` (CM by `Q(√−7)`) at **all 93
  primes up to 499** — `#C(F_p) = p + 1 − 3a_p` when `p ≡ ±1 mod 7`, `p + 1`
  otherwise;
- is the modular curve `X(7) ≅ X_H(49)` with `±H = {d ≡ ±1 mod 7}` (classical;
  Elkies).  Exhaustive coset enumeration reproduces genus 3, 168 cosets and
  24 cusps.

**The construction at `n = 131`.**

- Take the newform `f` of a CM curve over `Q` with CM by `Q(√−7)` that reduces
  to `E_0`.  The conductor-49 curve has `a_2 = +1` (it reduces to `E_1`), so twist
  it by `χ_{−3}`, which has `χ_{−3}(2) = −1`.  That gives level 441 and
  `a_2 = −1`.
- Let `χ` be the even character of order 131 and conductor 263.  2 has order
  131 mod 263 up to sign, so 2 is inert in `Q(ζ_263)^+` with residue field
  `F_2^131`.
- `f ⊗ χ` has coefficient field `Q(ζ_131)` and nebentypus `χ²`, at level
  `441·263² = 30,503,529`.
- Its abelian variety has dimension 130.  By Eichler–Shimura its Frobenius at
  2 has eigenvalues `χ(2)^σ·{α, ᾱ}` with `α² + α + 2 = 0`: exactly the Weil
  numbers of `A`.  Numerically, `∏_ζ (T² + ζT + 2ζ²)` matches `P_A(T)` at
  `T = 1, 2, 3` to 1e-13.
- It is a factor of `J_H(N)` with `±H = {d ≡ ±1 mod 263}`, of index 131 over
  `X_0(N)`, and `X_H(N)` has good reduction at 2.

| `n` | sign | level `N` | genus of `X_H(N)` | how | `dim A` | smallest genus found by search (toy) |
|---:|:--:|---:|---:|---|---:|---:|
| 3 | `E_1` | 49 | **3** (Klein) | enumeration | 2 | 3 |
| 3 | `E_0` | 441 | 121 | enumeration | 2 | **3** |
| 5 | `E_1` | 5,929 | 2,841 | enumeration | 4 | — |
| 5 | `E_0` | 53,361 | 36,001 | formula | 4 | > 4 |
| 131 | `E_1` | 3,389,281 | 42,307,761 (`2^25.33`) | formula | 130 | — |
| 131 | **`E_0`** | **30,503,529** | **508,799,809 (`2^28.92`)** | formula | 130 | — |

The formula rows rest on one fact: `Γ_0(N)` has no elliptic points at these
levels, and every cusp splits in `X_H → X_0`, so the cover is étale and
`g = 131·(g_0 − 1) + 1`.  The cover is cut out by a character of conductor
`ℓ`, and `ℓ` divides `N/gcd(c, N/c)` for every cusp denominator `c`.  This was
checked by enumeration at levels 49, 121, 441, 539 and 5,929 (all split). The
unconditional range at `n = 131` is `[508,799,809, 509,348,929]`.

Two readings:

- **The construction is explicit, and it is not small.**  It replaces "no curve
  is known below genus `2^129`" with an explicit curve at `2^28.92`.  The
  correspondence is too — the Hecke projector onto `f ⊗ χ` and the modular
  parametrisation — though nobody has an algorithm that evaluates it at this
  level.
- **At `n = 3` it is optimal for one sign and far from optimal for the
  other.**  For `E_1` it is the Klein quartic, genus 3.  For `E_0` the extra
  level 9 the sign costs makes it genus 121, while search finds genus 3.  So
  `2^28.92` is an upper bound on the smallest genus, not an estimate of it.

## 4. Route 2 — cyclic covers

The Klein quartic is also a cyclic cover.  `(x:y:z) ↦ (y:z:x)` preserves it and
is defined over `F_2`.  It has two fixed points, a conjugate pair over `F_4`,
and the quotient has genus 1.  So the natural generalisation is a cyclic cover
of degree 131.

**Derived bound.**  Let `C → B` be geometrically cyclic of prime degree
`n = 131`, `B` of genus `g_B ≤ 1`, with `r` totally ramified branch points.
Then `g(C) = n(g_B − 1) + 1 + 65r` and `Jac(C) ~ Jac(B) × P`.  Since `A` is simple
of dimension 130, `A ⊆ P`.  Frobenius acts on the group by a unit of order
`d | 130`.

- *Class field theory over `F_{2^d}`:*
  - the group's 131-part comes from `Pic⁰(B)` only if Frobenius has an
    eigenvalue of order `d` on `Jac(B)[131]`.  For the five elliptic classes
    that requires a root of `x² − tx + 2 mod 131` of order `d`, which exists only
    for `t = 0` at `d ∈ {65, 130}`;
  - otherwise it comes from a branch place whose residue field contains
    `μ_131`, so `r ≥ 130/d`.
- *`μ_d`-stability:*
  - Frobenius permutes the 130 non-trivial characters in orbits of length
    `d`, so `P`'s Weil numbers are stable under `μ_d`;
  - containing `{ζα, ζᾱ}` they contain `μ_d·{ζα, ζᾱ}`, which has `260d`
    elements, so `r ≥ 2(d + 1 − g_B)`.

| `d` | 1 | 2 | 5 | **10** | 13 | 26 | 65 | 130 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `r ≥` (CFT) | 130 | 65 | 26 | 13 | 10 | 5 | 2 | 2 |
| `r ≥` (`μ_d`) | 2 | 4 | 10 | **20** | 26 | 52 | 130 | 260 |
| genus ≥ (base `E_0`) | 8451 | 4226 | 1691 | **1301** | 1691 | 3381 | 8451 | 16,901 |

Over all six bases of genus ≤ 1 — `P¹` and the five elliptic isogeny classes —
the minimum is **1300** (`P¹`, `d = 10`); every elliptic base gives 1301.  The
two-branch-point family, genus 131, is excluded for every twist.

**Toy validation.**  The bound gives genus ≥ 3 at `n = 3` (Klein attains it)
and ≥ 9 at `n = 5`.  Every Kummer cover `y^ℓ = f` of `E_0` branched at
`r ≤ 4` rational points was computed — 44 covers at `ℓ ∈ {3, 5}`:

- every one is a Weil polynomial with `E_0` splitting off;
- every Prym polynomial lies in `Z[T^d]`, which is the `μ_d`-stability the
  bound uses, observed directly;
- none contains `A`, including the three `ℓ = 3, r = 4` covers the bound would
  allow.

**Fermat quotients.**  `y^m = x^a(1−x)^b` carries a factor isogenous to a power
of the `√−7` CM curve only if its Stickelberger CM type is induced from
`Q(√−7)`.  At `m = 7·263 = 1841` all 3.39 million pairs give only the 12
level-7 triples, the Klein-type ones.  No factor of level 263 or 1841
qualifies, so the genus-920 family is excluded.

**What route 2 does not close:** non-cyclic group actions, non-Galois covers,
and covers of curves of genus ≥ 2.  The toy's genus-3 curves for `E_0` have the
*twist* `E_1` as their elliptic factor, and the bound covers covers of `E_1` as
well.  So a genus-131 curve with `Jac ~ A × E_1`, the naive analogue, would have
to come from outside the cyclic family.

## 5. Route 3 — direct search, priced

At genus 130:

- hyperelliptic models `y² + hy = f` number `2^(3g+5) = 2^395`;
- the moduli of curves has dimension 387, the hyperelliptic locus 259.

Each candidate can be tested in polynomial time (p-adic point counting
exists for hyperelliptic curves in characteristic 2); the count is the
obstruction.  The toy ladder in §2 is this search run to completion where it
can be.

## 6. Route 4 — no curve

Index calculus directly on `A`, i.e. decomposition attacks on
`E_0(F_{2^131})`, is the repository's main Koblitz line.  The scoreboard's
verdict for the challenge family: every route tried (decomposition targets,
Weil-descent SAT, residual walks, `m ≥ 4` algebra) is bounded above rho or
closed.  The `m = 4` algebraic audit reads a cost slope of 0.985–1.06 bits per
unit `n` against the 0.25 it would need (`docs/index-calculus-scoreboard.html`).
The curve matters because a Jacobian turns smoothness into polynomial
factorisation, which is where the `2^37` comes from; an abelian variety with
no curve has no such structure.

## 7. Correction to the hyperelliptic note (accounting)

`RESEARCH_ECC2K130_HYPERELLIPTIC.md` said the genus any construction reaches over
`F_2` is 1, `2^129` or `2^130`, that the window and "the constructions" miss by
`2^120.77`, and that "every construction is closed".  The trichotomy is a
theorem about GHS/Hess descent and stands.  The rest overreached:

- the Klein quartic is a non-GHS construction at `n = 3`;
- the modular curve carries `A` at genus `2^28.92`, so the gap from the window
  to an explicit construction is `2^20.69`, not `2^120.77`.

That note now says so.  No measured or modelled number changed; this is an
accounting correction.

## 8. What this does not settle

- **The window question is open.**  Neither a construction nor a proof of
  non-existence is given for genus 130…~300 outside GHS descent and cyclic
  covers of genus ≤ 1 curves.
- **The cyclic bound is a necessary condition.**  It assumes prime-cyclic
  degree 131 and a base of genus ≤ 1; it does not show that covers meeting it
  carry `A`.
- **The modular genus uses the split-cusp fact.**  It is validated at five
  levels, not proved here.  Without it the range is `[508,799,809, 509,348,929]`.
- **Extrapolated costs are `L(1/2, √2)`.**  Calibrated against the exact model
  from genus 130 to 1300, they overestimate by 4.1–10.9 bits.  At genus `2^28.92`
  the figure means "astronomically above rho" and nothing finer.
- **Toy fidelity is structural, not quantitative.**  2 is primitive mod 3 and
  mod 5 as mod 131.  The toys still say nothing about genus ~130 directly.
- **Genus 5 is incomplete at `n = 5`.**  Only hyperelliptic models were
  enumerated.

## 9. Reproduction

```
cargo run --release --bin ecc2k130_curve_construction -- \
  --out research/notes/ecc2k130/curve_construction_20261006/results/curve_construction.json \
  --legacy experiments/ecc2k130_hyperelliptic_cover_boundary.json
```

Deterministic (seed 20261006 for the GHS census); about 19 s on a 4-core
container.  Unit tests:
`cargo test --release --lib curve_construction` (18 tests).  An uncommitted
Python scratch prototype ran the Klein check and the `n = 3, 5` hyperelliptic
searches before the protocol was frozen.  It is disclosed in the protocol and
not cited.

## References

- P. Gaudry, F. Hess, N. Smart, *Constructive and destructive facets of Weil
  descent on elliptic curves*, J. Cryptology 15 (2002).
- F. Hess, *Generalising the GHS attack on the elliptic curve discrete
  logarithm problem*, LMS J. Comput. Math. 7 (2004).
- C. Diem, *On the discrete logarithm problem in elliptic curves*,
  Compositio Math. 147 (2011).
- A. Enge, P. Gaudry, *A general framework for subexponential discrete logarithm
  algorithms*, Acta Arith. 102 (2002).
- G. Shimura, *Introduction to the arithmetic theory of automorphic functions*
  (1971) — `A_f`, Eichler–Shimura.
- E. Kani, M. Rosen, *Idempotent relations and factors of Jacobians*,
  Math. Ann. 284 (1989).
- N. Koblitz, D. Rohrlich, *Simple factors in the Jacobian of a Fermat curve*,
  Canad. J. Math. 30 (1978).
- J. Denef, F. Vercauteren, *An extension of Kedlaya's algorithm to
  hyperelliptic curves in characteristic 2*, J. Cryptology 19 (2006).
- N. Elkies, *The Klein quartic in number theory*, in *The Eightfold Way*,
  MSRI Publ. 35 (1998).
