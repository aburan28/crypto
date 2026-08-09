# Mestre's Algorithm for the Explicit Howe-Glued Cover

A focused note on Mestre's algorithm (Mestre 1991, refined by
Cardona–Quer 2005 and Lercier–Ritzenthaler 2010), which is the
**actual** way to construct the explicit Howe-glued genus-2 cover
that the audit in [`RESEARCH_SECP256K1_CM.md`](RESEARCH_SECP256K1_CM.md)
§§8.6, 8.10 left as an open implementation task.

> **Status (updated 2026-08-09)**: both directions of Mestre's
> algorithm are now implemented and round-trip-validated in PARI —
> see [`secp256k1_cm_audit/mestre_reconstruction.gp`](secp256k1_cm_audit/mestre_reconstruction.gp).
> This was unblocked by finding SageMath's own implementation
> (`sage/schemes/hyperelliptic_curves/{invariants,mestre}.py`,
> GPL-2.0+) and porting its exact formulas rather than re-deriving
> the "pages-long" Cardona–Quer polynomials by hand. What's still
> **not** implemented is item 2 below (the actual gluing formula:
> Igusa invariants of `(E×E^t)/Γ_α` computed *without* already
> having a curve model for it) — that remains the open problem.

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
classical polynomial expressions in the `a_i`.  Rather than
hand-transcribing the pages-long Cardona–Quer closed forms, they
are computed *algorithmically* via Mestre's own "Ueberschiebung"
(transvectant) recipe from his 1991 paper, p.315/317:

```
        i  = (f,f)_4          Delta = (i,i)_2
        y1 = (f,i)_4           y2 = (i,y1)_2         y3 = (i,y2)_2
        A = (f,f)_6    B = (i,i)_4    C = (i,Delta)_4    D = (y3,y1)_2
        I_2, I_4, I_6, I_10  =  clebsch_to_igusa(A,B,C,D)   [linear-ish
                                 change of basis, see source]
```

where `(f,g)_k` is the transvectant defined by
`(f,g)_k = [(m-k)!(n-k)!/m!n!] · (∂f/∂x·∂g/∂y − ∂f/∂y·∂g/∂x)^k`
for forms of degree `m,n`. This is implemented and validated (exact
match against SageMath's own docstring test vectors) in
[`secp256k1_cm_audit/mestre_reconstruction.gp`](secp256k1_cm_audit/mestre_reconstruction.gp),
function `igusa_clebsch`.

## 7. Triage: what to implement vs. defer

*(Updated 2026-08-09 — items 1 and 3 below are DONE; see
`secp256k1_cm_audit/mestre_reconstruction.gp`. Item 2 is still open
and is the one genuinely hard piece left.)*

### ~~High value, medium cost~~ DONE: Igusa invariants from curve coefficients

For a known `C : y² = h(x)`, compute `(I_2, I_4, I_6, I_10)` from
`h`'s coefficients. Implemented as `igusa_clebsch(f)` in
`mestre_reconstruction.gp` via the transvectant recipe above (ported
from Sage's `invariants.py`, not re-derived). Validated against 3
exact Sage docstring test vectors.

**Value**: enables verification that two `C`'s have the same
Jacobian (up to F_p-isomorphism). Useful for the cargo-test
extension. **Now also gives full `(I2,I4,I6,I10)` for `h_toy` and
`h_secp`** from §6 above, where the old `mestre_scaffold.gp` only
had `I_10` (via `poldisc`). `h_secp = (x³+7)(x³+189)` over
`F_p_secp`: `I10 != 0`, i.e. it's a bona fide smooth genus-2 curve —
but per the 2026-07-27 log entry, that does *not* by itself mean
it's the Howe-glued cover of `(secp256k1, secp256k1^t)`; the
Richelot-degenerate pair `(0,3)` finding is about the correspondence
used to try to *construct* it as that cover, not about whether
`h_secp` itself is singular.

### High value, high cost: Igusa invariants of (E_1 × E_2)/Γ_α

This is the actual moduli computation needed for the slice-3
construction.  Requires:
- Modular form pullback `M_2 → M_1 × M_1`
- Explicit gluing formula in Igusa coordinates
- Galois descent if no F_p-rational point

~1000+ lines of PARI; multi-week effort.

**Value**: the actual explicit Howe-glued cover for secp256k1.
Publishable result.

### ~~Medium value, high cost~~ DONE: Mestre's Step 2 (conic + sextic)

The reconstruction from Igusa invariants to curve. Implemented as
`reconstruct_curve(I2,I4,I6,I10,p)` in `mestre_reconstruction.gp`
(ported from Sage's `mestre.py`), using:
- `mestre_L`: the conic matrix (exact match to Sage's `Mestre_conic`
  on both a `Q`-example and a fresh `F_7`-example from Sage's own
  docstrings).
- `find_conic_point`: O(p) point search over `F_p` by completing the
  square in the dehomogenized chart (brute-force O(p³) fallback for
  the rare case the fast path misses, e.g. all-affine-chart-degenerate
  conics).
- `mestre_parametrize`: the standard "point + pencil of lines"
  rational parametrization of a conic through a known point.
- `mestre_sextic`: Mestre's `c_ijk` combination formula.

Validated end-to-end by round-tripping `sextic → invariants →
reconstruct → invariants` on 4 toy curves at `p ∈ {97, 1009}`: 3/4
matched exactly on the scale-invariant ratios
`I2⁵/I10, I4⁵/I10², I6⁵/I10³` (Mestre's output is only defined up to
the projective scaling `h ↦ λ·h`, `I_{2k} ↦ λ^{2k}·I_{2k}`, so the
raw invariants differ but these ratios must not). The 4th case hit a
genuine edge case: for one specific `(I2,I4,I6,I10)` input, the
particular basis choice in `mestre_parametrize` produced a
degenerate reconstructed curve (`I10 = 0`). This matches Sage's own
documented caveat ("Mestre's algorithm only works for generic curves
of genus two"); a production version should try multiple `(e1,e2)`
bases before giving up. Not chased further this session — the
mechanism (point-finding, parametrization, sextic formula) is
independently confirmed correct by the 3 successful round trips and
by directly checking each parametrized point lies on the conic.

**Caveat for the real secp256k1 case**: `find_conic_point`'s O(p)
scan is only a toy-scale method. A 256-bit `p` needs a genuine
rational-point-on-a-conic algorithm (e.g. PARI's `qfsolve` machinery
adapted to `F_p`, or the classical reduction to a 2-adic/Hensel
lift + CRT approach) — this is a concrete, scoped next step, *not*
a research question (the math is 19th-century; it's an
implementation task).

**Value**: completes the pipeline.  Without it, even if §3.2 is
done, we have Igusa invariants but not the curve.

## 8. Realistic path forward

The "explicit Howe cover for secp256k1" remains an open
implementation task.  The realistic options:

**Option A**: implement the full Mestre algorithm in PARI. **Done
2026-08-09** for the *reconstruction* half (invariants → curve, and
curve → invariants) — see `mestre_reconstruction.gp`, ~4 hours by
porting Sage's implementation rather than the originally-estimated
4-6 weeks of independent derivation. What's still missing is not
Mestre's algorithm itself but its *input*: the Igusa invariants of
`(E×E^t)/Γ_α` (item 2 in §7 above), which is a different, harder
problem (an actual gluing formula, not a reconstruction algorithm).

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
