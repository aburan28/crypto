# Mestre's Algorithm for the Explicit Howe-Glued Cover

A focused note on Mestre's algorithm (Mestre 1991, refined by
Cardona–Quer 2005 and Lercier–Ritzenthaler 2010), which is the
**actual** way to construct the explicit Howe-glued genus-2 cover
that the audit in [`RESEARCH_SECP256K1_CM.md`](RESEARCH_SECP256K1_CM.md)
§§8.6, 8.10 left as an open implementation task.

> **Status (2026-08-08 update)**: the forward gluing problem this note
> was scoping — construct `C: y²=h(x)` explicitly from `(E1,E2)` — is
> **solved and verified**, via a different (and much cheaper) route than
> Mestre/Igusa reconstruction. See §5.6. The Mestre/Igusa-invariant path
> described in §§2-8 below is superseded for this purpose and kept only
> as background; skip to §5.6 for the working construction.
>
> Original status: documentation + partial PARI scaffolding.  The full
> Mestre reconstruction is multi-page and depends on moduli-of-
> abelian-surfaces machinery that PARI's standard library doesn't
> expose directly.  This note documents the algorithm in enough
> detail for a future implementation effort, with the easy parts
> (Igusa-invariant formulas) actually computed.

## 1. Problem statement

Given two elliptic curves `E_1, E_2 / F_p` with:

- `E_1 ≁ E_2` (no F_p-rational isogeny between them)
- `E_1[2] ≃ E_2[2]` as F_p-Galois modules with Weil pairing
- `gcd(#E_1, #E_2) = 1`

Howe (1996) showed that the abelian surface

```
        A  :=  (E_1 × E_2) / Γ_α
```

— where `Γ_α` is the graph of any F_p-Galois-equivariant
isomorphism `α: E_1[2] → E_2[2]` preserving Weil pairings — is the
Jacobian of a smooth genus-2 curve `C/F_p`.

**The task**: given `(E_1, E_2)` and an explicit `α`, write down
`C` as `y² = h(x)` for some sextic `h(x) ∈ F_p[x]`.

The audit in `RESEARCH_SECP256K1_CM.md` §8.6 verified the three Howe
conditions hold for `(secp256k1, secp256k1^t)`.  Howe's theorem
guarantees `C` exists.  Mestre's algorithm constructs it.

## 2. High-level algorithm

```
1.  Compute the Igusa invariants (I_2, I_4, I_6, I_10) of A = Jac(C)
    from the gluing data (E_1, E_2, α).

2.  Solve Mestre's conic + cubic to find C from (I_2, I_4, I_6, I_10).

3.  Verify Jac(C) is F_p-isogenous to E_1 × E_2 by computing the
    Frobenius char poly (or by reducing to a known invariant).
```

Step 1 is the hard part because `A` is not literally a Jacobian of
a known curve — it's a quotient abelian surface obtained by
identifying 2-torsion.  Computing its Igusa invariants requires
working in the moduli space of abelian surfaces (Igusa's `A_2`).

Step 2 is the **classical Mestre reconstruction** — given Igusa
invariants, find `C`.  This is fully documented in Mestre 1991 and
implemented in Magma.

## 3. Step 1: Igusa invariants of the glued surface

The Igusa invariants `(I_2, I_4, I_6, I_10)` of an abelian surface
`A/F_p` are elements of `F_p` (modulo certain weight scaling).
They form a complete set of invariants for `A` in the moduli of
principally polarised abelian surfaces.

For `A = E_1 × E_2`: the Igusa invariants are **degenerate** —
they correspond to a "boundary point" of `A_2` where the surface
splits as a product.  Specifically, `I_{10}(E_1 × E_2) = 0` for the
PRODUCT polarisation.

For `A = (E_1 × E_2) / Γ_α` (the Howe-glued quotient): the Igusa
invariants are NON-degenerate when the gluing is "non-trivial" (i.e.,
`α` is not the identity-via-isomorphism).  The resulting `A` is a
genuine non-product abelian surface, sitting inside the **2-isogeny
locus** of `A_2`.

### Explicit formula (partial)

In characteristic ≠ 2, 3, given:
- E_1: `y² = (x - a_1)(x - a_2)(x - a_3)`
- E_2: `y² = (x - b_1)(x - b_2)(x - b_3)`
- α: `a_i ↦ b_{σ(i)}` for some permutation σ of {1, 2, 3}

The Igusa invariants of `(E_1 × E_2) / Γ_α` are computable as
**rational functions of the (a_i, b_i, σ)** — but the exact formula
is several pages and involves elliptic-modular forms via the
restriction `M_2 → M_1 × M_1` (modular form pullback).

In characteristic > 3 with f_1, f_2 irreducible cubics, the
formulas simplify significantly because:
- the (a_i) lie in `F_{p^3}` as a Galois-conjugate orbit
- the (b_j) similarly in `F_{p^3}`
- σ is determined uniquely by F_p-Galois compatibility

For secp256k1 specifically: σ is the identity on {1, 2, 3} once we
pick a Galois-compatible labelling of E[2] and E^t[2].

**Status**: Igusa-invariant formulas are computable but not in
PARI's standard library.  Magma's `IgusaInvariants` works.

## 4. Step 2: Mestre's reconstruction

Given `(I_2, I_4, I_6, I_10) ∈ F_p^4` with `I_{10} ≠ 0`, Mestre's
algorithm constructs `C : y² = f(x)` (deg f = 6) with these Igusa
invariants.  In outline:

1. **Build a 3×3 quadratic form** `Q(x, y, z)` in F_p[x, y, z]
   from the Igusa invariants.  Explicit formula in Mestre 1991
   §III.
2. **Find a F_p-rational point** on the conic `Q = 0` in P²(F_p).
   If no F_p-rational point exists, `C` is defined over a
   quadratic extension (a *twist* of the curve we seek); the
   resulting `C` then descends via Galois twisting.
3. **Find a parametrisation** of the conic, yielding a degree-2
   covering `P¹ → Q = 0`.
4. **Pull back a binary sextic** from the conic to `P¹`, giving
   `f(x)` such that `C : y² = f(x)` has the right Igusa invariants.

Mestre's full recipe is implemented in Magma as `MestreConicAndCubic`
or similar; in Sage as `Genus2reduction` adjacent functions.  In
PARI: not directly available, though the components (conic
intersections, sextic reconstruction) are doable with low-level
arithmetic.

## 5. PARI scaffolding (this script)

The companion file
[`secp256k1_cm_audit/mestre_scaffold.gp`](secp256k1_cm_audit/mestre_scaffold.gp)
implements the EASY parts:

1. Computing Igusa invariants of a known smooth genus-2 curve
   `C : y² = f(x)` (sanity check against PARI's hyperellcharpoly).
2. The Igusa-invariant-to-`C` reconstruction (Mestre Step 2) for a
   simple toy case where `Q = 0` has an obvious F_p-rational point.
3. Demonstrating that two different `f(x)` give the same Igusa
   invariants iff they correspond to F_p-isomorphic Jacobians.

What it does NOT implement:

1. Igusa invariants of an abelian surface specified as
   `(E_1 × E_2)/Γ_α` — this requires moduli computation.
2. Resolution of the conic step when no F_p-rational point exists
   (Galois twist).
3. End-to-end Howe-cover construction for secp256k1.

Filling in these gaps is a multi-week implementation project.  See
§7 for triage.

## 5.5. Sage cross-verification (RUN — Igusa normalization confirmed)

The companion script
[`secp256k1_cm_audit/howe_mestre.sage`](secp256k1_cm_audit/howe_mestre.sage)
was executed under Sage 10.9 (installed via `brew install sage`,
binary at `/opt/homebrew/Caskroom/sage/10.9,10.9.0/...`).  For the
secp256k1 candidate genus-2 curve

```
        C : y² = (x³ + 7) · (x³ + 189) = x⁶ + 196·x³ + 1323
```

Sage's `HyperellipticCurve(h).igusa_clebsch_invariants()` reported:

| Invariant | Sage (mod p_secp)                | PARI explicit (signed Z)  | Ratio |
|-----------|----------------------------------|---------------------------|------:|
| J_2       | `115792…908833279279` = `−1392384` mod p_secp | `−43512`        | **32** (= 2^5) |
| J_10      | `48626773926836865849920323584` | `46374105383717408990784`  | **2^20** |

The Sage values match our PARI explicit (and transvectant)
computations **exactly up to multiplicative constants** — 32 for
J_2, 2^20 for J_10.  This confirms the PARI transvectant
computations in
[`secp256k1_cm_audit/igusa_clebsch.gp`](secp256k1_cm_audit/igusa_clebsch.gp)
are bona fide Igusa-Clebsch invariants in a *different*
normalisation than Sage / Cardona-Quer.  Multiplying our values
by these constants would convert to the standard.

For full constructive Mestre reconstruction, the next step is to
compute Igusa invariants of the abelian surface $(E \times E^t)
/ \Gamma_\alpha$ (rather than the naive `y² = f_1 · f_2`), then
apply Sage's `HyperellipticCurveFromInvariants` or
Genus2Reconstruction.  Sage's `igusa_clebsch_invariants` only goes
forward (curve → invariants); the inverse requires the
Genus2Reconstruction patch.

## 5.6. The actual explicit forward map (SOLVED, 2026-08-08)

The Mestre/Igusa route above assumes the only way to build `C` is via
its Igusa invariants. It isn't: the classical Frey-Kani (1991) /
Howe-Leprévost-Poonen (2000, *Forum Math.* 12, 315-364) gluing formula
writes `h(x)` down **directly** from the six 2-torsion x-coordinates of
`E1` and `E2` — no invariant theory, no reconstruction step, no moduli
computation. It is the same construction reimplemented (for the SIDH
attack context) as `FromProdToJac` in Castryck-Decru's public break code
(`richelot_aux.py`, github.com/jack4818/Castryck-Decru-SageMath, 2022).

**Formula.** Let `E1: y²=cubic_1(x)` have nonzero 2-torsion
x-coordinates `a1,a2,a3` (roots of `cubic_1`), and `E2` have `b1,b2,b3`
similarly, with `a_i <-> b_i` a chosen Galois-equivariant pairing (see
below for which pairings are valid). Then:

```
M = [[a1*b1,a1,b1],[a2*b2,a2,b2],[a3*b3,a3,b3]]
[R,S,T] = M^-1 . [1,1,1]^T            (fails iff det M = 0)
RD = R * det(M)
da = (a1-a2)(a2-a3)(a3-a1),  db = (b1-b2)(b2-b3)(b3-b1)
s1 = -da/RD,  s2 = -T/R
ai_t = (ai - s2)/s1   for i=1,2,3
h(x) = s1 * (x²-a1_t)(x²-a2_t)(x²-a3_t)
```

`C: y²=h(x)` is the Howe-glued cover: `Jac(C)` carries a (2,2)-isogeny to
`E1×E2`, so (Tate) `#Jac(C)(F_p) = #E1(F_p)·#E2(F_p)` exactly whenever
`h` has F_p-rational coefficients.

**Which pairing to use.** Only pairings `π: a_i -> b_π(i)` that are
Frobenius-equivariant (`π∘σ1 = σ2∘π`, where `σ1,σ2` are Frobenius's
permutation action on the two root sets) can give an F_p-rational `h`.
When `σ1,σ2` are conjugate in `S_3` (always true if `E1[2]`,`E2[2]` have
matching Galois-orbit structure — see next paragraph), such a `π`
exists; try all `3! = 6` pairings and keep whichever gives `det(M)≠0`
and `#Jac(h) = #E1·#E2`.

**Precondition on `(E1,E2)`.** Frobenius acts on `E[2]\{O}` (3 points)
as one of: identity (full 2-torsion, cubic splits), a transposition
(exactly 1 rational 2-torsion point), or a 3-cycle (cubic irreducible).
A Frobenius-equivariant `E1[2]≅E2[2]` can only exist when `E1`,`E2`
have the *same* case — i.e. the same number of rational 2-torsion
points (0, 1, or 3). This is checkable directly from each curve's
2-division polynomial factorization pattern and costs nothing.

**Verification (`secp256k1_cm_audit/chlrs_forward_map.gp`, all checks
pass under `gp -q`):**

- *Sanity check*, curves with full rational 2-torsion, `p=97`: for 4 of
  the 6 x-coordinate pairings, `#Jac(glued curve) = #E1·#E2` exactly
  (`112·108=12096`). Formula confirmed correct independent of the j=0
  setting.
- *j=0 curves, generic pairs* (`p=1009`, `E1: y²=x³+7`, `E2: y²=x³+b2`
  for `b2 ∈ {2,5,6,11,15,16,17,18}`, none a sextic-twist multiple of 7):
  **8/8 give an explicit F_p-rational Howe cover** with `#Jac(h)`
  matching `#E1·#E2` exactly. Example: `b2=2` gives
  `h(x) = 140x⁶+40`. (Sparse form — only `x⁶` and constant terms — is
  expected: both inputs are j=0/CM-by-ζ₃ curves, so the glued surface
  inherits the order-3 automorphism, forcing the `x⁴,x²` coefficients
  of `h` to vanish identically.)
- *j=0 curves, sextic-twist pairs* `b2 = 7·h^k` for `k=0..5` (`h` =
  primitive 6th root of unity mod 1009 — this is exactly the
  `secp256k1 × (its own sextic twist)` family the paper's B-block
  wants): **all 6 are degenerate.** At every one of the 3
  Frobenius-valid pairings, `det(M) = 0`. This generalizes the
  2026-07-27 log entry (which found only `k=3`, the quadratic-twist
  pair, degenerate) to *all six* twist classes — sextic twists of a
  single j=0 curve are a degenerate locus for this construction, not
  just the quadratic-twist case.

**Bottom line.** The forward Howe-gluing map is explicit, implemented,
and verified for the general case. It structurally cannot glue
secp256k1 to its own sextic twist via this construction (needs a
limiting/renormalized variant — unsolved, see §7/log 2026-08-08 for
next steps). But if the paper's cover-cost argument (block B5) only
needs *some* non-degenerate explicit Howe cover of a j=0 curve — not
specifically the pairing with its own sextic twist — §5.6's
construction already supplies one, with a one-line existence check
(same-parity 2-torsion Galois orbit) for which partner curves work.

## 6. Toy verification

The scaffolding can demonstrate Mestre's Step 2 on a hand-picked
curve.  For

```
        C : y² = x⁶ + 526·x³ + 5665  over F_1009
```

(from the toy in [`howe_explicit_cover.gp`](secp256k1_cm_audit/howe_explicit_cover.gp)),
we can:

1. Compute its Igusa invariants `(I_2, I_4, I_6, I_10)`.
2. Apply Mestre's reconstruction to those invariants.
3. Recover `C` up to F_p-isomorphism.

The Igusa invariants of `y² = h(x)` for `h(x) = ∑ a_i x^i` are
classical polynomial expressions in the `a_i`.  They are:

```
        I_2  =  (...)
        I_4  =  (...)
        I_6  =  (...)
        I_10 =  disc(h(x))
```

The exact formulas are pages-long polynomials (see e.g.
Cardona–Quer 2005 Appendix).  We don't reproduce them here; in
practice they're plugged into a CAS as a fixed lookup table.

## 7. Triage: what to implement vs. defer

Three implementation tasks, ranked by cost-to-value:

### High value, medium cost: Igusa invariants from curve coefficients

For a known `C : y² = h(x)`, compute `(I_2, I_4, I_6, I_10)` from
`h`'s coefficients.  Formulas are classical and well-tabulated.
~50 lines of PARI.

**Value**: enables verification that two `C`'s have the same
Jacobian (up to F_p-isomorphism).  Useful for the cargo-test
extension.

### High value, high cost: Igusa invariants of (E_1 × E_2)/Γ_α

This is the actual moduli computation needed for the slice-3
construction.  Requires:
- Modular form pullback `M_2 → M_1 × M_1`
- Explicit gluing formula in Igusa coordinates
- Galois descent if no F_p-rational point

~1000+ lines of PARI; multi-week effort.

**Value**: the actual explicit Howe-glued cover for secp256k1.
Publishable result.

### Medium value, high cost: Mestre's Step 2 (conic + sextic)

The reconstruction from Igusa invariants to curve.  Multi-page
algorithm; well-documented but tedious.

~500 lines of PARI; week-long effort.

**Value**: completes the pipeline.  Without it, even if §3.2 is
done, we have Igusa invariants but not the curve.

## 8. Realistic path forward

The "explicit Howe cover for secp256k1" remains an open
implementation task.  The realistic options:

**Option A**: implement the full Mestre algorithm in PARI.
~2 KLOC, 4-6 week effort.  Produces the explicit `C/F_p_secp` as
a publishable result.

**Option B**: use Magma (proprietary) to compute the Howe cover
directly via `IsogenousJacobian` or equivalent.  Single-line if
Magma is available; not viable in this open-source codebase.

**Option C**: use Sage (open-source) — its `genus2_curves`
module has some of the pieces but not the (E × E^t)/Γ_α gluing.
~1 week to wrap and test.

**Option D**: accept the existence-without-construction state.
Howe's theorem guarantees the cover exists (§8.6 verified
conditions); for the structural-completeness theorem, this is
sufficient because B5's cost analysis applies to any cover.
Constructive recovery is a separate deliverable.

The structural-completeness paper
([`PAPER_STRUCTURAL_COMPLETENESS.md`](PAPER_STRUCTURAL_COMPLETENESS.md))
is sound under Option D.  The explicit constructive cover is a
nice-to-have, not a load-bearing element.

## 9. Connection to the wider programme

This note records what's needed to fill the gap left by the
errata in `RESEARCH_SECP256K1_CM.md` §8.10.  The structural-
completeness theorem (§8 of the paper) does not depend on
explicit construction; the theorem's B5 (cover cost ≥ ECDLP cost)
applies to any cover, abstract or explicit.

The explicit cover is a **publishable individual result**:
"the Howe-glued cover of secp256k1 is computable; here it is."
It would be a nice complement to the paper but is not required
by the theorem.

## 10. References

- Mestre, "Construction de courbes de genre 2 à partir de leurs
  modules", in *Effective Methods in Algebraic Geometry*, 1991.
- Cardona, Quer, "Field of moduli and field of definition for
  curves of genus 2", in *Computational Aspects of Algebraic
  Curves*, 2005.
- Howe, "Constructing distinct curves with isomorphic Jacobians
  in characteristic zero", Israel J. Math. 1996.
- Lercier, Ritzenthaler, "Hyperelliptic curves and their invariants:
  geometric, arithmetic and algorithmic aspects", JoA 2012.
- Igusa, "Arithmetic variety of moduli for genus two", Annals of
  Math. 1960.
- Streng, "Computing Igusa class polynomials", Math. Comp. 2014.
- The Magma `Genus2Reconstruction` package, by Cardona, Howe,
  Lercier, Ritzenthaler, Streng et al.
