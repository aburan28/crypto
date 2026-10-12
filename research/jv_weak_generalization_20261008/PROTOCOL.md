# Protocol: generalizing the Joux–Vitse weak class — which classes over F_{q³} admit a genus-3 cover at all

Frozen 2026-10-08, before any instrument is built.  **Stage diagnostic,
toy sizes.**  `S`, end-to-end cost and speedup are **unset**.  Nothing
here concerns a prime-field or deployed curve: the question lives over an
extension field of degree 3.  Status: **PENDING**.

## 1. Derivation, part A — the 2-torsion covers over F_{q³} are classified

Let `E : y² = f(x) = (x − e₀)(x − e₁)(x − e₂)` over `F = F_{q³}` with full
rational 2-torsion, `σ` the `q`-power Frobenius.  The Kummer-type descent
(GHS in odd characteristic, Diem's form) builds the compositum of the
three conjugate function fields `F(x, √(σ^k f))`, a `(Z/2)^m` cover of the
`x`-line with `m = dim_{F₂} span{[σ^k f]}` in `F(x)^× / squares`.  Each
place of the `x`-line at which some element of that span has odd valuation
ramifies, with `2^{m−1}` places above it of index 2, so Riemann–Hurwitz
gives

    g = 1 − 2^m + 2^{m−2} · R,

`R` the number of such places (finite roots plus `∞` when some element has
odd degree).  The roots of the conjugates are the Frobenius orbits of the
`eᵢ`, of size 1 (`eᵢ ∈ F_q`) or 3.  Every configuration:

| configuration of `{e₀, e₁, e₂}` | `m` | `R` | `g` |
|:--|--:|--:|--:|
| all three in `F_q` | 1 | 4 | 1 (the curve is a base change) |
| one `ρ ∈ F_q`, the other two in **one** orbit `{α, σα}` | 3 | 5 | **3** (the Joux–Vitse class) |
| one in `F_q`, the other two in two orbits | 3 | 8 | 9 |
| all three in one orbit `{α, σα, σ²α}` | ≤ 2 | 4 | 1 (a twist of a base change) |
| two in one orbit, the third in another | 3 | 7 | 7 |
| three distinct orbits (generic) | 3 | 10 | 13 |

The Joux–Vitse condition `N_{F/F_q}` of a cross-ratio `= 1` is exactly
"some model has the second configuration" (ledger §13), and `g = 3` occurs
for it alone.  So **among 2-torsion-type covers over `F_{q³}`, the weak
class is complete**: there is no other configuration with `g ≤ 3`, and
`g ≤ 3` is what a genus-`g` index calculus at `q^{2−2/g}` needs to beat rho
at `q^{3/2}`.  Curves with one or no rational 2-torsion point give `g ∈
{1, 9, 13}` by the same count.  This part is derivation, checked below
on the census's own curves (A-1).

## 2. Derivation, part B — covers that are not of 2-torsion type

The Weil restriction `A = Res_{F/F_q}(E)` is a 3-dimensional principally
polarised abelian variety over `F_q` with characteristic polynomial
`x⁶ − t x³ + q³`, `t` the trace of `E` over `F`.  A genus-3 curve `C/F_q`
has `Jac(C)` isogenous to `A` over `F_q` exactly when its `L`-polynomial is
`1 − tT³ + q³T⁶`, i.e. when

    #C(F_q) = q + 1,   #C(F_{q²}) = q² + 1,   #C(F_{q³}) = q³ + 1 − 3t.

Such a `C` with a map to `E` is a genus-3 cover whether or not `C` is
hyperelliptic; the plane-quartic case is the one Diem's index calculus
for non-hyperelliptic genus 3 handles, and the one the (ℓ,ℓ,ℓ)-isogeny
route of Tian (J. Cryptology 2023) makes explicit for some curves.  The
cover ledger records that [JV12]'s sieve needs the hyperelliptic form.
So the **largest class a genus-3 route could ever reach** is the set of
traces `t` realised by some genus-3 curve over `F_q` with `a₁ = a₂ = 0`,
and the **new** part is the traces realised only by plane quartics.

## 3. Instrument (Rust, standalone, `genus3.rs`)

`q` a prime with `q ≡ 1 (mod 12)`, `F_{q³} = F_q[s]/(s³ − c)`, so the
weak test and the census run as in `../jv_weak_class_law_20261007/census.rs`
with `q` prime in place of `p²`.

1. **Traces of `E/F_{q³}`**: every Legendre `λ` gives the full-2-torsion
   classes with their weak flag (the JV set `W`); the set `T₂` of their
   traces by BSGS; separately the set `T` of all traces in the Hasse
   interval realised by some curve (a sample of 4,000 random curves).
2. **Hyperelliptic genus 3 over `F_q`**: `y² = h(x)`, `h` squarefree of
   degree 7 or 8, sampled uniformly (`N_h` curves); `#C(F_q)` and
   `#C(F_{q²})` by direct counting; survivors with the signature get
   `#C(F_{q³})` and `t = (q³ + 1 − #C(F_{q³}))/3`.  Set `H` of traces.
3. **Plane quartics over `F_q`**: random ternary quartics, nonsingular
   (checked by the gradient system having no common projective zero),
   `N_c` curves; the same counts; set `Q` of traces.
4. Report `|W|`, `|T₂|`, `|T|`, `|H|`, `|Q|`, `|H ∖ W|`, `|Q ∖ W|`,
   `|Q ∖ (W ∪ H)|`, with the depth `v₂(f)` and `D_K mod 8` of every trace
   in `Q ∖ W`.

Sizes: `q ∈ {13, 37, 61}`; `N_h = N_c = 200,000` at `q = 13`, `50,000` at
`q = 37`, `20,000` at `q = 61`.  Seed `20261008`.

## 4. Predictions (pass/fail)

- **A-1 (classification check).**  On every full-2-torsion curve of the
  census, the orbit configuration's `g` from the table equals the
  Riemann–Hurwitz count computed from the actual span and ramification;
  and `g = 3` iff the JV test passes.
- **A-2 (hyperelliptic traces).**  `H ⊆ W`: every trace of a
  hyperelliptic genus-3 curve with the signature is a JV-weak trace.
  *Falsified by one trace in `H ∖ W`*, which would be a new explicit
  hyperelliptic class.
- **A-3 (new classes exist).**  `Q ∖ W` is non-empty at every size, and
  `|Q ∖ W| / |Q| ≥ 0.3`: plane quartics reach traces the weak class does
  not.
- **A-4 (depth-1 is excluded here too).**  No trace in `Q` has `v₂(f) = 1`.
  *This is the sharpest test of the depth-1 law: it asks whether the
  obstruction is a property of genus-3 Jacobians isogenous to `Res(E)`,
  not of the hyperelliptic form.*
- **A-5 (coverage).**  `|(W ∪ Q)| / |T|` is reported; no threshold.

## 5. Decision rule

A-2 passing confirms part A on the hyperelliptic side.  A-3 passing
names the new class: curves over `F_{q³}` whose Weil restriction is
isogenous to a plane-quartic Jacobian, reachable in principle by Diem's
non-hyperelliptic genus-3 index calculus, with explicitness the open
problem (Tian's isogenies).  A-4 passing or failing decides whether the
depth-1 exclusion is a theorem about `Res(E)` or about the hyperelliptic
form, which fixes what a proof must prove.  Class: **boundary** throughout;
no attack cost is claimed.

## 6. Inadmissible

Treating a trace in `Q` as an explicit attack; sampling quartics
non-uniformly to find signatures; omitting the `F_{q³}` count; reading
anything here as a statement about a prime-field curve.

[JV12] A. Joux, V. Vitse, EUROCRYPT 2012.  Diem, *An index calculus
algorithm for plane curves of small degree*, ANTS 2006.  S. Tian, *Cover
attacks for elliptic curves over cubic extension fields*, J. Cryptology
2023.
