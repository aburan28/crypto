# Isogenies, conductor gaps and ECDLP hardness differences in the ECC2K-130 class

**Date:** 2026-10-07
**Status:** computed structure (exact, reproducible with
`volcano_ecc2k130.py` below) plus research angles with falsifiers. No DLP
measured. Literature attributions are from memory and should be checked.
**Companions:** `PLAN_IC_ACCOUNTING_FIXES_20261007.md`,
`RESEARCH_IC_NOVEL_DIRECTIONS_20261007.md`.

## −1. Prior work in this workspace (found after the computation below was done)

`/Volumes/SSD990/ecdlp-hardness-work/report/REPORT.md` and
`/Volumes/SSD990/ecdlp-hardness-work/ic-conductor-leads/LEADS.md` (both
2026-10-07) derive the same class structure independently, build and verify
all 262 conductor-263 curves at `n = 131` (PARI `polclass` and Vélu from
`E0[263]` over `F_{q²}`, agreeing on all j-invariants), give rho hardness per
level (crater 60.81 bits, floor 64.33 native and 60.81 effective, `p` levels
`[60.81, 64.33]` conditional on a nonexistent characteristic-2 Kani-type
pipeline; Galbraith ePrint 2024/924), estimate the cost of constructing one
`p`-level curve at `2^70–2^71.3` `F_q`-mults, and measure toy index calculus
across levels with an exact oracle: the only level-dependent IC resource is
`τ` (ordered τ-slots, `log₂(m!)` bits). The numbers in §1–§2 agree with
theirs where they overlap. The angles in §3 are restated against that work in
`EXPERIMENTS_ISOGENY_CONDUCTOR_GAP_20261007.md`, which marks V2 as done
there and V7 as closed there, and targets what remains: the production
solvers of this repository on non-crater curves, the landed rungs
`n = 41, 61, 83`, the 263 coincidence, the horizontal cycle structure, and
the characteristic-2 Kani question.

## 0. The theory being probed

Jao–Miller–Venkatesan (ASIACRYPT 2005) prove, under GRH, that the ECDLP is
random self-reducible among curves with the **same endomorphism ring**: the
horizontal isogeny graph is an expander, so a polynomial number of
small-degree isogeny steps moves an instance to a random curve on the same
volcano level. Vertical steps change the endomorphism ring; a step at a prime
`ell` costs polynomial in `ell` (modular polynomial, Elkies/Kohel), or
`Õ(√ell)` with √élu (Bernstein–De Feo–Leroux–Smith 2020) when a kernel point
is at hand. JMV's open case: two curves in the same isogeny class whose
endomorphism-ring conductors differ by a **large prime** have no known
reduction in either direction. Their hardness could, as far as theory knows,
be drastically different. The one concrete historical example of hardness
varying inside an isogeny class is Galbraith–Hess–Smart (2002): over
`F_{2^155}` a small fraction of each class is GHS-weak and horizontal walks
find them. That example needs composite `n`; Menezes–Qu (2001) close it for
prime `n`.

The "conductor gap" below means the index of `Z[π]` in `End(E)`: for the
challenge curve `End = O_K` and the gap is the full conductor `f` of `Z[π]`.

## 1. The ECC2K-130 isogeny class, computed

`K_0: y² + xy = x³ + 1` over `F_{2^131}`, `τ² + τ + 2 = 0`, `End(K_0) = Z[τ] = O_K`,
`K = Q(√−7)`, `h(O_K) = 1`. With `π = τ^131`:

```
t = tr(π) = −22283658519494248867
#E = 2^131 + 1 − t = 4 · 680564733841876926932320129493409985129 = 4r
t² − 4q = −7 · f²,  f = cond(Z[π] ⊂ O_K) = 38531015900842053623  (66 bits)
f = 263 · 146505763881528721          (both prime)
```

Volcano profile (K_0 is alone on the surface since `h(O_K) = 1`):

| prime `ell` | depth `v_ell(f)` | `(−7/ell)` | horizontal edges at surface | curves one level down | `π mod ell` | field degree of kernel points |
|---|---:|---:|---:|---:|---:|---:|
| 263 | 1 | +1 (split) | 2 (the endomorphisms of norm 263) | 262 | −1 | **2** |
| 146505763881528721 ≈ 2^57.3 | 1 | −1 (inert) | 0 | 146505763881528722 | 65861677189822723 | **12208813656794060 ≈ 2^53.4** |

Total `F_q`-isomorphism classes in the isogeny class: `38531015900842054149 ≈ 2^65`.

Consequences, all exact:

1. **Two hardness components.** The class splits into `{K_0}` ∪ `{262 curves at
   the 263-level}` (reachable: kernel points live in `F_{q²}`, so a 263-isogeny
   is a Vélu computation on 131 kernel multiples) and `{≈ 3.85·10^19 curves}`
   at or below the `ell_big`-level, separated from `K_0` by one isogeny of prime
   degree `≈ 2^57` whose kernel points live in an extension of degree `≈ 2^53`.
   No known method computes that isogeny: the modular polynomial has degree
   `2^57`, the kernel field is unreachable, and √élu needs the kernel. Within
   the big component, JMV applies horizontally (2 splits in `O_K` and `2 ∤ f`,
   so 2-isogenies are all horizontal there). So the ECC2K-130 class realizes
   JMV's open case exactly: **99.99999999999999% of the curves of order `4r`
   over this field are not known to be DLP-equivalent to the challenge curve.**
2. **Those curves cannot even be written down.** Constructing a curve with
   `End = Z + ell·O_K` by CM needs the ring class polynomial of degree
   `ell + 1 ≈ 2^57`; random sampling hits order `4r` with probability
   `≈ 2^−66`. The only members of the class anyone can exhibit are `K_0` and
   the 262 curves at the 263-level. The "drastic difference" is therefore
   unobservable at `n = 131` in either direction.
3. **Horizontal moves from `K_0` are endomorphisms.** `h(O_K) = 1`, so every
   horizontal isogeny from the surface is an element of `O_K` (τ, τ̄, the two
   norm-263 elements, …). JMV random self-reducibility is vacuous for `K_0`;
   a GHS-style horizontal walk cannot leave the curve.
4. **Every known weakness is a class or field invariant except `τ`.** GHS is
   dead for all curves over `F_{2^131}` (prime degree, Menezes–Qu); the
   embedding degree of `r` and the anomalous condition depend only on `(q, r)`;
   the Frobenius-stable-subspace obstruction for index calculus is a property
   of `F_{2^131}` as an `F_2[τ]`-module, not of the curve. The one property
   that varies across the class is the degree-2 endomorphism `τ`, which only
   `K_0` has (`τ ∉ Z + f'·O_K` for `f' > 1`). So isogenies can only make the
   ECDLP **harder** relative to `K_0`: rho loses the `√(2n)` signed-Frobenius
   quotient; index calculus loses the `2n`-fold orbit folding in relations and
   the `(2n)²` in linear algebra. `K_0` is the easiest curve in its class by
   every known measure.
5. **The 263-level in detail.** The 2-power Frobenius `τ` has eigenvalues
   `123, 139 mod 263` on `E[263]`, ratio of order 131, so the 262 descending
   kernels form **two `F_2`-Galois orbits of 131 curves each**; none of those
   curves is defined over `F_2`. `π ≡ −1 mod 263` means `E[263] ⊂ E(F_{2^262})`
   is fully rational over the quadratic extension, `263 | 2^262 − 1`, and the
   Tate pairing on `E(F_{2^262})[263]` is non-degenerate, which is exactly the
   Ionica–Joux "pairing the volcano" setting. The prime-order subgroup `⟨G⟩` is
   trivial in `E/263E`, so this pairing says nothing about the DLP of interest.
   Curiosities: `263 = 2n + 1` is also the type-II ONB prime, and
   `r ≡ 1 mod 263` and `r ≡ 1 mod 131`.

## 2. The toy ladder has the same shape at tractable sizes

Same computation for the `a = 0` rungs the repos use:

| `n` | `f = cond(Z[π])` | factorization | kernel field degree per prime | curves in class |
|---|---:|---|---|---:|
| 23 | 967 | 967 | 21 | ≈ 2^9 |
| 41 | 703889 | 409 · 1721 | 408, 215 | ≈ 2^19 |
| 53 | 68476319 | prime | 646003 | ≈ 2^26 |
| 61 | 1146312001 | 1951 · 587551 | 975, 117510 | ≈ 2^30 |
| 71 | 31797598073 | prime | ≈ 3.2·10^10 | ≈ 2^34 |
| 73 | 23369856751 | prime | ≈ 1.2·10^10 | ≈ 2^34 |
| 83 | 347450761417 | 6473 · 53676929 | 6472, ≈ 1.3·10^7 | ≈ 2^38 |
| 97 | 264777735625649 | 751943 · 352124743 | 375971, ≈ 5.9·10^7 | ≈ 2^47 |

So `n = 23, 41, 61, 83` have a level that is constructible (small prime, small
kernel field) **and** solvable by rho on both sides, which is what makes
hardness-difference experiments possible at all; `n = 71, 73` reproduce the
`n = 131` situation (one unbridgeable large prime). The twist `K_1` has the
same `f` (it depends on `t²`), hence the same volcano shape with a different
group order.

## 3. Angles

Each angle names the computation or experiment, the pre-registered
expectation, and what a deviation would mean.

### V1. Record the volcano profile in every candidate manifest

The `cryptanalysis/AGENTS.md` manifest already has `endomorphism` fields
(conductor, `v_ell` per prime, proof status). Fill them for every rung from
the table above and for `n = 131`, with `f`, its factorization, and the
reachability verdict per prime (kernel field degree, isogeny cost class).
*Expectation:* purely descriptive, but it makes "ISO1 transport" proposals
check against the reachability column before anything is run.

### V2. Constructive equivalence on the 263-component at `n = 131`

Compute one descending 263-isogeny `φ: K_0 → E_1`: pick a random point of
`E(F_{2^262})`, multiply by `#E(F_{2^262})/263²`-type cofactors to land in
`E[263]`, pick a non-eigen kernel, apply Vélu over `F_{2^262}` (131 kernel
multiples), and reduce the resulting curve to `F_{2^131}`. Verify `#E_1 = 4r`
by transporting `G` and checking `[r]φ(G) = O` and `φ(G) ≠ O`. Record both
Galois orbits (two inequivalent `E_1` up to `F_2`-conjugation).
*Expectation:* succeeds in minutes. Deliverable: the only two non-Koblitz
curves of order `4r` anyone can exhibit, with their DLP constructively
equivalent to the challenge via one isogeny evaluation on `(G, Q)`.

### V3. Hardness-difference panel on the toy ladder (the measurable version of JMV's question)

At `n = 23` (ell = 967), `n = 41` (409 and 1721) and `n = 61` (1951), construct
`E_1` one level down as in V2, transport `(G, Q)`, and run both pipelines on
`K_0` and on `E_1` with matched parameters: rho with negation only versus
signed-Frobenius, and index calculus with the same factor-base size `B`.
*Pre-registered expectations:* rho operation ratio `E_1 / K_0 ≈ √(2n)`
(6.8, 9.1, 11.0); IC relations needed `× 2n`, LA `× (2n)²`, **per-solve PDP
cost unchanged** (Semaev's `S₃` for `y² + xy = x³ + a₂x² + b` is
`(x₁x₂ + x₁x₃ + x₂x₃)² + x₁x₂x₃ + b`: same five terms, `a₂` absent, so only the
constant `b` leaves `F_2`); first-fall and solving degree identical within
noise. *A deviation* in the PDP solving degree between the two curves would be
the first measured sign that algebraic hardness depends on the curve model
beyond the Frobenius structure, which is exactly what IDEA-434 asks.

### V4. Hardness within the unbridgeable component, measured where it can be

At `n = 41` the 1721-level (kernel degree 215) and at `n = 61` the 587551-level
(kernel degree 117510) are reachable with effort. Compare rho and IC on a curve
at the 409-level versus the 1721-level versus the bottom (both primes) at
`n = 41`. *Expectation:* no difference beyond the loss of `τ` (which all three
share); any difference between levels would contradict JMV-style equivalence
inside the class and must be replicated before being believed.

### V5. Decide, for a given curve of order `4r`, which component it is in

An attacker handed an arbitrary curve `E'` with `#E' = 4r` should first test
whether it is `K_0` or one of the 262 reachable curves (compute `j(E')`,
compare against the two Galois orbits from V2; or try to ascend at 263 via the
rational ascending kernel). If not, `E'` is in the big component and the
Koblitz speedups are unavailable; Bisson–Sutherland endomorphism-ring
computation certifies the level. *Expectation:* a random curve of order `4r`
is in the big component with probability `1 − 2^−57`. This is the honest
statement of "same order does not imply same hardness" for this challenge.

### V6. Why the `ell_big` gap is not just expensive but structurally closed

Write down, for the descending isogeny at `ell_big ≈ 2^57` from `K_0`: the
modular-polynomial route (`Φ_ell` has `≈ 2^114` coefficients), the kernel
route (field degree `≈ 2^53`), the CM route (ring class polynomial of degree
`2^57`), and the Elkies/root-finding route (needs a root of `Φ_ell(X, j)`).
All four are `≫ 2^61`, the rho cost. *Research question worth stating
plainly:* is there any representation of an `F_q`-rational prime-degree
isogeny from a CM surface curve whose cost is polylogarithmic in `ell` when
`π` acts as a scalar on `E[ell]`? Nothing known; a positive answer would
collapse the two components and would be a result about isogenies, not about
ECC2K-130.

### V7. Vertical relation collection (negative, record it)

Pulling relations back through `φ̂φ = [263]` from the 263-level yields nothing
new: relations among `φ`-images are relations on `K_0`, and relations on `E_1`
among non-image points introduce unknown logs. The 2-power torsion is
unchanged by odd-degree isogenies (`E_1(F_q) ≅ Z/4 × Z/r`), so the FHJRV
symmetry group has the same size on every curve in the class. *Expectation:*
no gain; one paragraph in the ledger closes "ISO1 as an IC lever" for this
class.

### V8. The 263 coincidence

`263 = 2n + 1`, `263 | f`, `263 | r − 1`, `μ_263 ⊂ F_{2^131}^×`, and
`E[263] ⊂ E(F_{2^262})`. Check whether this is forced (it is at least
half-likely: `τ/τ̄` has order dividing 131 in `F_263^×` with probability about
1/2) or accidental, and whether the Weil pairing `E[263] × E[263] → μ_263`
interacts with the cyclotomic-sparse factor base of the directions note.
*Expectation:* accidental and inert for the DLP (the pairing never sees
`⟨G⟩`), but it costs an afternoon and it closes a door that will otherwise be
reopened.

### V9. Cost table for vertical steps across the ladder

For each `(n, ell)` in §2 give the measured cost of one descent by the kernel
route (`F_{q^d}` arithmetic with `d = ord(π mod ell)`, √élu or Vélu) and the
modular-polynomial route, and the cost of one DLP transfer (evaluating the
isogeny on two points). *Expectation:* kernel route feasible for
`d ≤ 10^3`, modular route for `ell ≤ 10^4`; this table is the operational
meaning of "conductor gap" and belongs next to V1.

### V10. Equidistribution check inside the reachable component

JMV's expander statement is about horizontal graphs at a fixed level. The
263-level at `n = 131` has 262 curves connected horizontally by 2-isogenies
(all horizontal there, since `2 ∤ f`). Compute the 2-isogeny cycle structure
on those 262 curves (it is the action of the ideal above 2 in
`Cl(Z + 263·O_K)`, a group of order `262`). *Expectation:* one or two cycles;
purely structural, but it is the smallest complete JMV component anyone can
draw for a cryptographic-size curve.

### V11. Pre-register what would count as a hardness difference

A measured ratio of operation counts between two curves in the same class
that is not explained by (a) the `√(2n)` or `2n` Frobenius factors, (b) the
factor-base size, or (c) implementation speed, with both arms measured under
the accounting fixes in the plan note. Anything else is noise. This angle
exists so that V3 and V4 cannot be read as "found a difference" when they
found the Frobenius factor.

## 4. Reproduction

```python
# volcano_ecc2k130.py  (sympy)
# Lucas sequences for tau, taubar with tau+taubar = -1, tau*taubar = 2:
# V_{k+1} = -V_k - 2 V_{k-1}, U_{k+1} = -U_k - 2 U_{k-1};  t = V_n,  f = |U_n|
# check t*t - 4*2**n == -7*f*f ;  factorint(f) ; chi(ell) = legendre(-7, ell) (chi(2)=+1)
# curves one level down at ell: ell - chi(ell);  kernel field degree: n_order(t/2 mod ell, ell)
```

The scripts used for the numbers above are checked in at
`research/isogeny_conductor_gap_20261007/` (`volcano_ecc2k130.py`,
`volcano2.py`, `volcano_output.txt`). The first script has a cosmetic
formatting bug in its final print loop; the second supersedes it and produced
the tables in §1 and §2.
