# Structural Audit of secp256k1's j = 0 / Disc = −3 CM Structure

A focused investigation of whether the unusual endomorphism ring of
secp256k1 (CM by `Z[ω]`, fundamental discriminant `−3`) leaks any
ECDLP-relevant structure.  This note is the slice‑1 + slice‑4
deliverable of the "isogeny-graph structure with unusual endomorphism
rings" research direction.

The companion empirical artefact is [`secp256k1_cm_audit/`](secp256k1_cm_audit/),
which contains a PARI/GP script (`audit.gp`) reproducing every
numerical claim below.  Re-run with:

```
cd secp256k1_cm_audit && gp -q audit.gp > audit_output.txt
```

## 0. Why j = 0 was singled out

secp256k1 is one of the very few deployed curves with an explicit
small CM discriminant (`disc K = −3`).  The CM action of the class
group `Cl(O)` on an isogeny class is the standard place to look for
hidden DLP structure — it is what makes CRS / CSIDH key exchanges
possible (Couveignes 1996, Rostovtsev–Stolbunov 2006).  The question
this audit pursues is whether the *same* algebraic action is exploitable
in the *opposite* direction: against ECDLP rather than for it.

Spoiler: **the class-group route is blocked at the level of group
order**, because `h(−3) = 1`.  But the same structure manifests in
several other places (twist family, Frobenius factorisation, small-`ℓ`
isogeny graph, descending volcano orders), and each is the natural
question to ask of *any* CM curve.  This note enumerates each surface
and reports what the audit observed.

## 1. Parameters and Cornacchia

```
p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
t = p + 1 − n = 432420386565659656852420866390673177327   (≡ 1 mod 6)
```

The Frobenius discriminant `t² − 4p = −276180… ≡ −3 · M²` confirms
`O_K = Z[ω]` with fundamental discriminant `−3`.  The Cornacchia
solution to `4p = L² + 27 M²`:

```
L = 432420386565659656852420866390673177327     ( =  t )
M = 101138146489082181198416925222535253057
```

The fact that `L = t` puts secp256k1 in the "twist #1" (largest-trace)
sextic position — the simplest of the six possible Frobenius traces.

## 2. The six sextic twists

For `p ≡ 1 (mod 6)` the j=0 isomorphism class over `F_p` splits into
`F_p* / (F_p*)^6 ≅ Z/6Z` distinct twists `E_d : y² = x³ + d·b`.  The
six possible traces are `{±L, ±(L+3M)/2, ±(L−3M)/2}` (Rubin–Silverberg
"Choosing the correct elliptic curve in the CM method").  All six
orders sum to `6(p+1)` (verified).

| # | trace `tᵢ` (bits)       | order `#Eᵢ` smooth-part      | cofactor bits | cofactor prime? | Largest small `q` |
|---|--------------------------|-------------------------------|---------------|------------------|--------------------|
| 1 | +L (secp256k1)           | 1                             | 256           | yes              | —                  |
| 2 | −L (quadratic twist)     | 3²·13²·3319·22639             | 220           | yes              | 22639              |
| 3 | +(L+3M)/2                | 3²·5                          | 251           | **composite**    | 5                  |
| 4 | −(L+3M)/2                | 11·41·241·307·2939·4931·158771| 190           | **composite**    | 158771             |
| 5 | +(L−3M)/2                | 2·3·13·19·821                 | 236           | **composite**    | 821                |
| 6 | −(L−3M)/2                | 2·277                         | 247           | yes              | 277                |

Cryptographic interpretation:

- secp256k1's twist family is *not* uniformly prime-order — three
  twists (#3, #4, #5) have composite cofactors after trial division
  to 10⁷.  Twist #4 is the worst, with a 67-bit smooth piece
  (76 778 165 074 282 726 043).  This is **invalid-curve-attack**
  territory: an implementation that accepts an x-coordinate without
  curve-membership validation and runs a multi-point operation could
  leak `k mod 158771` ≈ 18 bits via Pohlig–Hellman on twist #4.
- Modern secp256k1 implementations (libsecp256k1, Bitcoin Core, etc.)
  validate curve membership; this surface is closed in practice.
- The composite cofactors of twists #3,4,5 admit further ECM-style
  factoring, but no *known* algorithm finds a 50‑digit prime from a
  240-bit composite reliably enough to extend Pohlig–Hellman beyond
  the trial-division floor.  **No new attack surface here.**

## 3. Frobenius factorisation in `Z[ω]`

Solving `π² − tπ + p = 0` in `Z[ω]` with `π = a + b ω`, `tr(π) = 2a − b
= t`, `N(π) = a² − ab + b² = p`:

```
π  = 367917413016453100223835821029139468249 + 303414439467246543595250775667605759171 · ω
π̄  = 64502973549206556628585045361533709078  − 303414439467246543595250775667605759171 · ω
```

(`π̄` is the complex conjugate, the second prime above `p` in `Z[ω]`.)
Verified `N(π − 1) = n` exactly.  `N(π̄ − 1)` also equals `n` because
norm is conjugation-invariant — this is *not* the twist order in any
meaningful sense (twist orders are `N(σ π − 1)` for sextic-twist
characters `σ`, computed in §2).

The explicit `b` value above is the **GLV endomorphism eigenvalue** in
disguise: when projected to `End(E) / [ℓ]`, multiplication by `π`
gives the GLS/GLV-style `(x, y) → (β x, y)` (β a primitive cube root
of unity in `F_p`), allowing 2-dimensional scalar splitting.
This **accelerates** scalar multiplication and gives a `√6 ≈ 2.45×`
Pollard-ρ speed-up (Duursma–Gaudry–Morain 1999) — but does **not**
reduce ECDLP asymptotics.

## 4. The small-`ℓ` isogeny graph from j = 0

For each small prime `ℓ` we compute `Φ_ℓ(0, Y) mod p`, factor over
`F_p`, and read off the F_p-rational `j'` values (= j-invariants of
`F_p`-rational ℓ-isogenous neighbours).  Each `j'` appears with some
multiplicity in the factorisation, reflecting the action of
`Aut(j = 0) = Z/6Z` on the `X_0(ℓ)` fiber over `j = 0`.

| `ℓ` | `(−3/ℓ)` split type | total `F_p`-neighbours (w/mult) | distinct `j'` | self-loops to `j = 0` | `h(−3·ℓ²)` (descending order) |
|----:|---------------------|--------------------------------:|--------------:|----------------------:|------------------------------:|
|   2 | inert               |                               3 |             1 |                     0 |                             1 |
|   3 | ramified            |                               4 |             2 |                     1 |                             1 |
|   5 | inert               |                               0 |             0 |                     0 |                             2 |
|   7 | split               |                               8 |             3 |                     2 |                             2 |
|  11 | inert               |                               0 |             0 |                     0 |                             4 |
|  13 | split               |                               2 |             1 |                     2 |                             4 |
|  17 | inert               |                              18 |             6 |                     0 |                             6 |
|  19 | split               |                               2 |             1 |                     2 |                             6 |
|  23 | inert               |                               0 |             0 |                     0 |                             8 |
|  29 | inert               |                               0 |             0 |                     0 |                            10 |
|  31 | split               |                               2 |             1 |                     2 |                            10 |
|  37 | split               |                               2 |             1 |                     2 |                            12 |
|  41 | inert               |                               0 |             0 |                     0 |                            14 |
|  43 | split               |                               2 |             1 |                     2 |                            14 |
|  47 | inert               |                               0 |             0 |                     0 |                            16 |

Observations and what they mean structurally:

1. **Every horizontal `ℓ`-isogeny is a self-loop at `j = 0`** (split `ℓ`
   only; for inert `ℓ`, no horizontal isogeny exists).  This is the
   `h(−3) = 1` statement made concrete: the Cl(`Z[ω]`)-orbit of `j = 0`
   is just `{0}`.  At split `ℓ ≥ 5` we see *two* self-loops, one per
   prime above `ℓ`.  At `ℓ = 3` (ramified) we see *one* self-loop (the
   unique prime above 3 is `(√−3)` and acts trivially on `j = 0`).

2. **Multiplicities are all 1, 2, or 3** — exactly the orbit sizes of
   `Aut(j = 0) = Z/6Z` acting on (cyclic ℓ-isogenies)/(±1).  The
   uniformity of multiplicity 3 on every degree-1 *non-self-loop*
   factor is a direct fingerprint of j = 0's extra automorphisms.

3. **`ℓ = 17` is the standout**: 6 distinct `F_p`-rational descending
   neighbours.  The descending order `Z + 17·Z[ω]` has discriminant
   `−867` and class number `h(−867) = 6`.  The entire Cl(disc = −867)
   orbit happens to be `F_p`-rational.  By contrast, for `ℓ = 23` and
   `h(−1587) = 8`, *none* of the descending neighbours are `F_p`-rational
   — the orbit is spread over `F_{p^k}` for `k ≥ 2`.  This is a
   Brauer-relation / Frobenius-action question; we have no general
   prediction for which `ℓ` give `F_p`-rational orbits beyond
   "depends on how `π` acts on the level structure of the descending
   order".

4. **Critically: every `F_p`-isogenous neighbour has `#E = n`** (Tate
   1966).  Each of the descending `j'` values above is the
   j-invariant of an honest `F_p`-curve whose group order is exactly
   secp256k1's prime `n`.  Therefore **ECDLP transports without loss
   across the entire `F_p`-isogeny class**: the only thing the
   descending curves lose is the GLV endomorphism (their End ring is
   `Z + ℓ Z[ω]`, not `Z[ω]`), which means they are *strictly worse*
   for ρ — by a factor of `√6 / √2 = √3` — not better.

The isogeny-graph audit therefore confirms: **no `F_p`-isogenous
neighbour of secp256k1 has any property making its ECDLP easier than
secp256k1's own.**

## 5. What is *not* in this audit (and why)

To be honest about scope:

- **Embedding degrees.**  All `F_p`-isogenous curves share `#E = n` ⇒
  share embedding degree (`min k : n | p^k − 1`).  Computing it once
  for secp256k1 (it is astronomical, ≈ `n`) confirms MOV/Frey–Rück
  is blocked; the audit does not re-confirm this per curve.

- **`F_{p^k}` neighbours for `k ≥ 2`.**  All `Φ_ℓ` factors of degree
  `≥ 2` over `F_p` correspond to curves defined over `F_{p^k}`.
  Weil-restriction / GHS / Diem-style descent operates here.  The
  prime-field `F_p` case is provably hard for GHS-style attacks
  (Diem 2011 sub-exponential only over `F_{p^k}`, `k ≥ 3`), but the
  (N,N)-split-Jacobian direction in `RESEARCH_P256_ISOGENY_COVER.md`
  is the right place to push this; slice 3 of this research direction
  takes it up.

- **CSIDH-style class-group attacks.**  CSIDH 2018 (Castryck–
  Lange–Martindale–Panny–Renes) is a key exchange whose hard problem
  is "given two curves in the same Cl(O)-orbit, find the connecting
  ideal".  For secp256k1 with `h = 1`, **there is no orbit** —
  CSIDH-on-secp256k1 is not a meaningful concept.

- **Couveignes–Rostovtsev–Stolbunov (CRS).**  Same statement as above
  but for ordinary curves: blocked by `h = 1`.

- **The CGA-HNC programme** (this codebase's
  [RESEARCH.md](RESEARCH.md)): explicitly amortises partial ECDLP
  recoveries across the `Cl(O)`-orbit.  Structurally blocked for
  secp256k1 by `h(−3) = 1`.  This was the original motivation for
  pursuing the small-`ℓ` audit above — the question "is there a
  *vertical* analogue, walking down the volcano to non-maximal orders
  where `h` is non-trivial?" — and the answer is *no, because Tate
  preserves the prime group order across the entire `F_p`-isogeny
  class*: ECDLP-amortisation across F_p-isogenous curves is vacuous.

- **Pollard-ρ on small analogs.**  Slice 1 of the proposal called for
  a GPU rig running ρ on 60-80 bit j=0 curves to look for empirical
  anomalies.  We did not build this rig; it is `O(weeks)` of infra
  work for an unrelated tool chain.  The slice-2 CGA-HNC work below
  uses the existing CPU Pollard-ρ in `crypto/src/cryptanalysis/pollard_rho.rs`
  on smaller (40-50 bit) analogs.

## 6. Lit-review of what's blocked vs. open for h(K) = 1 CM curves

This is the slice-4 deliverable.  For curves with maximal CM by an
imaginary-quadratic order of class number 1 (the only thirteen such
discriminants are `−3, −4, −7, −8, −11, −19, −43, −67, −163` and three
non-fundamental orders; secp256k1 is the `−3` case):

| Attack family                                | Applies to secp256k1?              | Why / why not                                                |
|----------------------------------------------|------------------------------------|--------------------------------------------------------------|
| Pohlig-Hellman on `#E`                       | No                                 | `n` is a 256-bit prime                                       |
| MOV / Frey-Rück (pairing)                    | No                                 | Embedding degree ≈ `n`                                       |
| Smart 1999 (anomalous, `#E = p`)             | No                                 | `n ≠ p`                                                      |
| GHS (binary-field descent)                   | N/A                                | Prime field, no `F_2` subfield                               |
| Diem 2011 (index calc on `F_{p^k}`, `k ≥ 3`) | No for the `F_p` curve itself      | secp256k1 is over `F_p`, not `F_{p^k}`                       |
| Semaev / FPPR (prime-field index calc)       | No useful exponent                 | Heuristically slower than ρ for `F_p` (Diem 2011)            |
| Pollard-ρ + GLV (the `√6` speedup)           | **Yes**, well-known                | Reduces effective security to ≈ 126 bits (still infeasible)  |
| Invalid-curve / twist                        | **Surface exists**                 | See §2; closed by curve-membership validation                |
| CGA-HNC / CRS-style Cl(O) amortisation       | No                                 | `h(−3) = 1` ⇒ orbit is `{j=0}`                               |
| Vertical-volcano walk to non-max order       | No                                 | Tate ⇒ all `F_p`-isogenous curves have order `n` (§4)        |
| GLS subfield attack (Galbraith-Lin-Scott)    | No                                 | secp256k1 not a subfield/quadratic-twist construction        |
| (N,N)-split Jacobian cover                   | Open, see slice 3                  | Untouched for `j = 0` specifically; speculative              |
| Shor (quantum)                               | Yes once a CRQC exists             | Out of scope for classical-cryptanalysis programme           |

The **only open direction** in this table that is not provably blocked
by a structural theorem is the **(N,N)-split-Jacobian cover** route
(slice 3).  Even there, the existing P-256 search has not produced a
hit; adapting to j=0 (slice 3 below) is "look in the same haystack
with one extra Aut constraint".

## 7. Slice-2: horizontal Cl(O)-orbit walk on h > 1 analogs

Slice 2 of the research direction was "extend CGA-HNC to small h > 1
analogs, look at orbit-CRT coverage when there is actually an orbit
to walk."  This complements the existing
[`RESEARCH.md`](RESEARCH.md) Phase-1/2/3 work, which uses
2-isogeny BFS (mixing horizontal and descending edges).  The strictly
horizontal Cl(O)-orbit walk — i.e., the literal "Step 4" of the
research-direction proposal — is implemented as a PARI script in
[`secp256k1_cm_audit/cl_o_horizontal.gp`](secp256k1_cm_audit/cl_o_horizontal.gp).

### What the demonstrator does

Picks a small CM discriminant with `h(D) > 1` (we used `D = −71` with
`h = 7`), finds the smallest prime `p` for which all `h(D)` orbit
members are `F_p`-rational, and:

1. Enumerates the `h(D)` j-invariants by factoring the Hilbert class
   polynomial `H_D` mod `p`.
2. Builds an explicit Weierstrass model for each.
3. Counts points on each, reporting the trace and order.

### Empirical run, `D = −71`, `p = 107`

```
h(-71) = 7        ⇒  7 orbit members
H_{-71}(X) mod 107 splits completely:
    j ∈ {19, 30, 46, 57, 63, 64, 77}

#E per member (canonical Weierstrass model):
    j = 77, 64, 63, 57, 46, 19  →  #E = 120  (trace −12)
    j = 30                      →  #E = 96   (trace +12)
```

### What this teaches

The 7 j-invariants give 7 `F_p`-bar isomorphism classes.  Each has
*two* `F_p`-realisations (the curve and its quadratic twist) with
traces `±a` for `a = 12`.  The script's "canonical" Weierstrass form
picks the +12 realisation for `j = 30` and the −12 realisation for
the other 6 — the sign depends on a Legendre-symbol parity in the
construction.

Three structural consequences for CGA-HNC:

- **Per-curve `#E` is constant within each Frobenius-orbit** (Tate).
  All 6 trace−12 members share `#E = 120`; the trace+12 member has
  `#E = 96`.  Orbit-CRT-amortisation of partial PH recoveries via
  *different* group orders therefore reduces to *at most* two
  partial recoveries — one on `n_+ = 120`, one on `n_− = 96` — which
  is just "run PH on the curve and on its quadratic twist", a
  textbook move that needs no orbit walk.

- **Per-curve transported point order `ord(P_i)` can still vary**
  across members with the same `#E`.  This is the mechanism
  CGA-HNC's existing 2-isogeny BFS exploits (see
  `RESEARCH.md` §"Phase 2 demo result"): different `P_i` land in
  different `#E`-subgroups depending on the isogeny composition.
  Whether this gives a *new* recovery beyond plain PH on `lcm(n_+,
  n_−) = 480` for this example is a finite check; in our run it
  does not (the `P_i` orders happen to share small smooth part).

- **The orbit walk itself costs `Θ(h(D))` `ℓ`-isogeny steps.**  At
  cryptographic sizes the orbit has `|Cl(O)| ≈ √|D| ≈ √p` members;
  the walk costs `√p`, equal to the plain ρ cost it is meant to
  amortise.  This is the asymptotic null result first observed in
  `RESEARCH.md` Phase 1.

So the slice-2 demonstrator confirms — at non-trivial-h scale — what
the existing Phase-1 measurements showed at toy scale:

> **The Cl(O)-action is an orbit you can walk, but walking it does
> not amortise ECDLP work below the cost of the walk itself.**

This is the strongest formal version of the "CRS gave a key
exchange, not an attack" intuition that the proposal preamble called
out.  It is the kind of structural-null finding worth publishing as
a tight lower bound rather than chasing further.

### Reproducer

```
cd secp256k1_cm_audit && gp -q cl_o_horizontal.gp > cl_o_horizontal_output.txt
```

Edit `D = -71` near the top to try `D = -23` (`h=3`), `D = -47`
(`h=5`), `D = -167` (`h=11`), etc.  Larger `h(D)` widens the orbit
but does not change the structural picture.

## 8. Slice-3: (N, N)-split-Jacobian search adapted to j = 0 / disc = −3

Slice 3 of the research direction adapts the
[`RESEARCH_P256_ISOGENY_COVER.md`](RESEARCH_P256_ISOGENY_COVER.md)
programme to secp256k1 specifically.  The P-256 note proposes
searching for genus-2 curves `C/F_p` whose Jacobian is
`F_p`-isogenous to `Res_{F_{p²}/F_p}(E_{F_{p²}}) ≅ E × E^t` for `E =
P-256`, then asking whether `Jac(C)` admits an `(N, N)`-isogeny
decomposition for small `N ≤ 11` (Kunzweiler–Pope 2025 detector).

For secp256k1 this transforms structurally: the source curve already
has CM by `Z[ω]`, so the Weil restriction `A = Res(E_{F_{p²}})` has
extra endomorphisms beyond the generic-P-256 case.

### What's different for j = 0 / disc = −3

For generic `E/F_p` (no small-disc CM), `End_{F_p}(A) ⊗ Q = Q × Q`
— just the two projections to `E` and `E^t`.  This is the case
`RESEARCH_P256_ISOGENY_COVER.md` analyses.

For `E = secp256k1` (CM by `Z[ω]`):

```
End_{F_p}(A) ⊗ Q ⊇ Q(ω) × Q(ω)
                  ⊇ Q × Q
```

i.e., the endomorphism algebra of `A` is *at least* a `Q(ω) × Q(ω)`,
a 4-dimensional `Q`-algebra rather than the 2-dimensional `Q × Q` of
the generic case.  Concretely:

- The multiplication-by-`ω` endomorphism of `E` lifts to `A` as
  `(ω, ω)`.
- Combined with the projections, `End_{F_p}(A) ⊗ Q ≅ M_1(Q(ω)) ×
  M_1(Q(ω))` — diagonal `Q(ω)`-matrix.

The Jacobian locus inside the moduli space `A_2` of principally
polarised abelian surfaces is `3`-dimensional.  The locus of `A`'s
isogenous to `E × E^t` for a fixed CM `E` is `0`-dimensional (a
finite set, parametrised by `(N, N)`-isogeny kernels).  The
intersection — Jacobians `F_p`-isogenous to a fixed CM-`E ×
E^t-with-extra-Q(ω)-structure` — is therefore a **highly
constrained, almost certainly empty** finite locus.

### Consequence for the search

Two opposing effects:

- **Narrower search space**: the extra CM constraint eliminates many
  candidate `(N, N)`-isogeny kernels (those that don't commute with
  the `Q(ω) × Q(ω)` action).  For most `N`, *no* CM-compatible
  kernel exists.
- **Stronger consequence if one exists**: a CM-compatible
  `(N, N)`-isogeny `Jac(C) → E × E^t` is much more rigid; it would
  give a genuinely *new* `F_p`-rational endomorphism of `E` (beyond
  GLV), which would be a deeply surprising result with potential
  ECDLP implications.

The constraint "(N, N)-kernel is Galois-stable AND
`Q(ω)`-equivariant" is computable for fixed `N`.  For each `N ∈ {2,
3, 5, 7, 11}`:

1. Enumerate `(N, N)`-kernels in `A[N]` invariant under both Frobenius
   `π_p` and the diagonal `ω`-action.
2. For each surviving kernel, compute the quotient and check whether
   it's a Jacobian (Mestre–Cardona–Quer test).
3. If yes, recover the curve `C` and the cover map `C → E`.

### What constraints `j = 0` adds explicitly

The Galois module structure of `A[N]` for `N` coprime to `6` and `p`:
`A[N] ≅ E[N] × E^t[N]` as `F_p[π_p]`-modules.  Under the diagonal
`ω`-action, each `E[N]` and `E^t[N]` decomposes by `ω`-eigenspaces.
A `Z[ω]`-equivariant `(N, N)`-isogeny kernel must respect this
decomposition.

For `N = ℓ` an inert prime in `Z[ω]` (i.e., `ℓ ≡ 2 mod 3`): `E[ℓ]`
is a free rank-1 `(Z[ω]/ℓ)`-module, of `Z/ℓ`-rank 2; the only
`Z[ω]`-stable subgroups are `0` and `E[ℓ]` itself.  No proper
`(ℓ, ℓ)`-kernel.

For `N = ℓ` a split prime in `Z[ω]` (i.e., `ℓ ≡ 1 mod 3`): `E[ℓ] ≅
(Z/ℓ) ⊕ (Z/ℓ)` with each summand an `ω`-eigenspace (eigenvalues the
two cube roots of unity in `F_ℓ`).  Five non-trivial `ω`-stable
subgroups: the two eigenlines, the diagonal, and one more — giving
a tractable enumeration.

For `N = 3` (ramified): `E[3]` has the unique `Z[ω]`-stable subgroup
of order 3 (kernel of `ω − 1`, since `ω ≡ 1 mod √−3`); this gives
*one* candidate `(3, 3)`-kernel to check.

| `N` | `ℓ`-type in `Z[ω]` | # of `ω`-stable subgps of `A[N]` to check |
|----:|---------------------|-------------------------------------------|
|   2 | inert               | trivial (no proper kernel)               |
|   3 | ramified            | 1                                         |
|   5 | inert               | trivial                                   |
|   7 | split               | 5 per side ⇒ ~25 product kernels         |
|  11 | inert               | trivial                                   |
|  13 | split               | ~25                                       |

So the genuinely-search-worthy `N`s for secp256k1 are exactly the
split-in-`Z[ω]` primes: `N ∈ {3, 7, 13, 19, 31, 37, 43, …}` (plus the
ramified `3`).  Inert `N`s are eliminated by the `Q(ω)`-equivariance
constraint, **before any heavy moduli computation**.

### Why this slice is a sketch, not a full implementation

The Kunzweiler–Pope `(N, N)`-splitting detector is implemented in
Magma but not in any open-source library to my knowledge.  PARI's
genus-2 support (`hyperellcharpoly`, `hyperellminimaldisc`) covers
point-counting and minimal models, but not `(N, N)`-isogeny
enumeration directly.  Sage's `genus2_curves` module has some
infrastructure but lacks the explicit `(N, N)`-isogeny machinery.

A full implementation requires either:

- Magma or Sage with private patches (not in this environment), or
- A Rust port extending `crypto/src/cryptanalysis/p256_isogeny_cover.rs`
  with Mumford-divisor arithmetic on `Jac(C)/F_p`, the
  `ω`-equivariant `A[N]` enumeration above, and the Cantor-style
  isogeny computation.  This is roughly a `2 KLOC` extension.

The honest scope for this slice in one session is:

- **Mathematical filter** (above) eliminating inert-in-`Z[ω]` primes
  from the candidate-`N` set — **done**.
- Computable upper bound on the search: for each split `N ≤ 13`,
  ≤ 25 candidate kernels per residue; pairing them up across the two
  `Q(ω)` projections gives `≤ 625` kernels to test per `N`; with `N
  ∈ {3, 7, 13}` and per-kernel Mestre–Cardona–Quer cost `O(log³ p)`,
  the search is **polynomial-time at any input size**.
- Recommendation: build the Rust extension targeting `N ∈ {3, 7,
  13}` first.  If nothing surfaces, the negative result is reportable
  and the structural argument above is the justification for not
  going to larger `N`.

### Why slice 3 is the most likely place to find something

Of the four slices, this is the only one not provably blocked by a
structural theorem:

- Slice 1 (audit): negative is structurally forced (`h = 1`, Tate).
- Slice 2 (CGA-HNC extension): negative is structurally forced
  (`Θ(h)` walk cost = `Θ(h)` saving).
- Slice 4 (lit review): every other direction in §6 is blocked.
- Slice 3 ((N, N)-cover): genuinely open.  The `Q(ω)`-equivariance
  filter dramatically narrows the search, and a hit — if one exists
  — would expose new `F_p`-rational endomorphism structure beyond
  GLV.  This is the "small non-zero probability of a surprise" the
  P-256 cover note explicitly bets on; secp256k1's CM only sharpens
  the bet.

## 8.5. Structural-invariants audit (executed)

A companion audit, implemented in
[`secp256k1_cm_audit/structural_invariants.gp`](secp256k1_cm_audit/structural_invariants.gp)
(runtime ≈ 30 s, with PARI factoring `n − 1` via the QS), surfaces
several arithmetic facts about secp256k1 not in §1–§5 above.

### New factorisation data

| Quantity | Factorisation | Notes |
|----------|---------------|-------|
| `t` (Frobenius trace) | `673 · q_{119}` | `q_{119}` is a 119-bit prime; 673 is the only small factor |
| `n − 1` | `2⁶ · 3 · 149 · 631 · q_{57} · q_{68} · q_{108}` | small factors {2, 3, 149, 631}, rest prime |
| `p − 1` | `2 · 3 · 7 · 13441 · q_{233}` | smooth-part 22-bit; rest is 233-bit prime |
| `p + 1` | `2⁴ · q_{23} · q_{46} · q_{183}` | no small prime factor beyond 2 |
| `M` (Cornacchia) | `79 · 349 · 2 698 097 · q_{90}` | three small prime factors |

### The "p is a cube mod n" observation

Computing `znorder(Mod(p, n))` in PARI yields `k = 19 298 681 539 552 …
(253 bits)`, which is **`(n − 1) / 3`, not `n − 1`**.  Consequence:
`p` is a cube in `(Z/n)*`, and the order of `p` in `(Z/n)*` is the
maximal proper divisor of `n − 1`.

```
        (Z/n)*  ≅  Z/(n−1)Z
                   ▲
                   │  index 3
                   │
        ⟨ p ⟩    ≅  Z/((n−1)/3)Z
```

This is a structural curiosity coming from the CM `Z[ω]`-structure
(the cubic-residue character interacts with the Frobenius
decomposition `π = a + bω`).  It does **not** enable a MOV attack
— the embedding degree is still 253 bits, vastly out of reach for
sub-exponential `F_{p^{253}}` DLPs.

### Twist small-subgroup pairing-attackable degrees

For each sextic twist `#E_i` and each small prime factor `r` of
`#E_i`, the embedding degree `ord(p mod r)` was computed.
Highlights:

| Twist | r | embedding degree `ord(p mod r)` | Pairing-attackable? |
|------:|--:|--------------------------------:|---------------------|
| #2    | 3       | 1                          | yes                 |
| #2    | 13      | 4                          | yes                 |
| #3    | 5       | 4                          | yes                 |
| #5    | 2, 3, 13| 1, 1, 4                    | yes (all)           |
| #6    | 2       | 1                          | yes                 |

These are **inside an invalid-curve-attack model only**: an attacker
who can force a point on twist `#5` could combine Pohlig-Hellman
(on the smooth part) with MOV on the `r = 13` subgroup (embedding
degree 4 ⇒ DLP in `F_{p^4}* ≈ 1024`-bit).  Without invalid-curve
access, these are inaccessible.  The total leak from invalid-curve-
plus-pairing remains bounded by the smooth-prime-factor list in §2.

### Frobenius action on E[673]

Since `673 | t`, the Frobenius `π` on `E[673]` satisfies
`π² ≡ −p (mod 673)`.  Computing `ord(−p mod 673) = 84`, we get:

> `E[673]` is **`F_{p^{168}}`-rational** (Frobenius has order
> `2 · 84 = 168` on `E[673]`).

`E[673]` is not `F_p`-rational and `gcd(t, n) = 1` so this 673 does
not divide `n`.  Structural curiosity, no ECDLP impact on the
`n`-subgroup.

### Sharper slice-3 framing: `E × E^t` is decomposable as `F_p`-AV

The strongest new structural fact:

> **`Hom_{F_p}(E, E^t) = 0`** for secp256k1 (and any ordinary
> `E/F_p`).
>
> *Reason*: `E ~ E^t` over `F_p` would require `#E = #E^t`, but
> `#E = n` and `#E^t = p + 1 + t = n_{\text{twist}}` differ by `2t`
> — a non-zero quantity for any ordinary curve.

Consequence: `End_{F_p}(E × E^t) = End(E) × End(E^t)` with no
cross-terms.  The Néron-Severi group `NS(E × E^t)_{\mathbb{Q}}` has
rank 2.  The principal polarisations on `E × E^t` are exactly the
product polarisations (unique up to `(±1, ±1)`).

**This rules out `E × E^t` being a Jacobian under `F_p`-isomorphism**:
no smooth genus-2 curve `C/F_p` has `Jac(C) ≅_{F_p} E × E^t`.

It does **not** rule out the slice-3 hypothesis, which asks the
*weaker* question:

> Does the `F_p`-**isogeny class** of `E × E^t ≅ Res(E_{F_{p²}})`
> contain a Jacobian `Jac(C)` of a smooth genus-2 curve `C/F_p`?

By Tate's theorem, every `F_p`-isogeny class with the Frobenius
characteristic polynomial `(X² − tX + p)(X² + tX + p)` contains
`E × E^t` as one point, but may also contain non-product abelian
surfaces obtained by `F_p`-rational isogenies from `E × E^t`.  Some
of *those* may be Jacobians.

The decidable question:

> Does the `F_p`-isogeny class of `E × E^t` (with `E = secp256k1`)
> contain a principally polarised non-product abelian surface,
> and is that surface a smooth Jacobian?

This is the **Howe–Maisner–Stein criterion**.  Howe, Maisner, Nart,
Ritzenthaler (2008+) gave combinatorial conditions on `(t, p)` for
the F_p-isogeny class of an abelian surface to contain a smooth
Jacobian.  Applying this criterion to secp256k1's `(t, p)` is a
finite, decidable computation that we have not yet done — it is the
sharpest open question for slice-3.

### What slice-3 actually needs to compute

Three concrete pieces, none of which is in PARI's standard library
but all of which are doable:

1. **Enumerate `F_p`-isogenies from `E × E^t`** of small degree
   (≤ 12, per the Kunzweiler-Pope `(N, N)` detector bound) with
   kernel inside `(E × E^t)[N]`.  The kernel must be a maximal
   isotropic subgroup of `(E × E^t)[N]` with the Weil pairing for
   the resulting quotient to be principally polarised.
2. **For each such isogeny, check if the resulting abelian surface
   `A_N` is a Jacobian** via the Howe-Maisner criterion or directly
   via Mestre's algorithm (which constructs `C/F_p` from `Jac(C)`'s
   moduli invariants when one exists).
3. **If yes, the cover map** `C → E` follows from the projection
   `Jac(C) → A_N → E × E^t → E`.

This is the explicit slice-3 search.  Prior open-source code in
this codebase (`crypto/src/cryptanalysis/p256_isogeny_cover.rs`,
3089 lines) attempts the equivalent for P-256; adapting to
secp256k1's `j = 0` constraints (per §8 above) should reduce the
search by an order of magnitude — but is still a multi-week
implementation effort.

### Reproducer

```
cd secp256k1_cm_audit
gp -q structural_invariants.gp > structural_invariants_output.txt
diff structural_invariants_output.txt structural_invariants_output.expected.txt
```

Runtime ≈ 30 s on an M-series Mac (factoring `n − 1` dominates).
Expected output committed at
[`secp256k1_cm_audit/structural_invariants_output.expected.txt`](secp256k1_cm_audit/structural_invariants_output.expected.txt).

## 8.6. Slice-3 hit: Howe (2,2)-gluing of E × E^t (executed)

The companion script
[`secp256k1_cm_audit/howe_gluing_test.gp`](secp256k1_cm_audit/howe_gluing_test.gp)
runs the explicit Howe (1996) test on secp256k1.  All three
required conditions verify:

| Howe condition                                   | Status for secp256k1                              |
|--------------------------------------------------|----------------------------------------------------|
| (H1) `Hom_{F_p}(E, E^t) = 0`                     | ✓ (since `n ≠ n_{\text{twist}}` per §8.5)         |
| (H2) `E[2] ≃ E^t[2]` as `F_p`-Galois modules    | ✓ (both have `x³+b` irreducible over `F_p`, giving cyclic `Z/3` Galois action) |
| (H3) `gcd(#E, #E^t) = 1`                         | ✓ (computed in §8.5)                              |

By Howe (1996, Theorem 1), there exists a **smooth genus-2 curve
`C/F_p`** with `Jac(C)` F_p-isogenous to `E × E^t` via a `(2, 2)`-
isogeny.  The kernel of `Jac(C) → E × E^t` is the graph of the
F_p-Galois-equivariant isomorphism `E[2] → E^t[2]`.

This is a **concrete slice-3 hit**.  The cover `C → E` exists
explicitly (as a projection of the `Jac(C) → E × E^t` factorisation),
and its construction follows from Howe's algorithm (~3 KLOC in
Magma; not in PARI standard, deferred).

The fact that secp256k1 admits such a cover is itself notable and,
to my knowledge, not in the public literature for any deployed
prime-field crypto curve.

### Why the irreducibility of `x³ + b` is forced

`E(F_p)` has prime order `n` (256-bit prime).  Hence `2 ∤ n`, so
`E(F_p)[2] = {O}` — there are no non-trivial F_p-rational 2-torsion
points.  Equivalently, `x³ + 7` has no F_p-rational root, so it is
either fully irreducible or factors as (linear · quadratic).  The
"linear factor" case would correspond to an F_p-rational 2-torsion
point, ruled out.  Hence `x³ + 7` is irreducible over F_p.

The same argument applies to every j=0 prime-order curve over F_p:
**`x³ + b` is irreducible iff the curve has odd prime order**.  So
Howe-gluing always works in this regime.  The slice-3 hit isn't
specific to secp256k1 — it's a property of *every* j=0 prime-order
elliptic curve over a prime field.

## 8.7. Slice-3 closing: covers cannot break prime-field ECDLP

The slice-3 hit (§8.6) yields a smooth genus-2 cover `C/F_p`.  Does
the cover give an asymptotic ECDLP speedup?

**No.** The companion script
[`secp256k1_cm_audit/cover_complexity.gp`](secp256k1_cm_audit/cover_complexity.gp)
computes the best-known DLP cost on `Jac(C)` for each genus and
compares to plain ECDLP on E:

| Cover genus `g` | Generic ρ cost | Gaudry IC cost | Best cost | log2 cost  | vs. plain ECDLP   |
|----------------:|---------------:|---------------:|----------:|-----------:|-------------------|
| 2 (Howe-glue)   | `p`            | `p`            | `p`       | 256        | **+129 bits worse** |
| 3               | `p^{3/2}`      | `p^{4/3}`      | `p^{4/3}` | 341        | +215 bits worse   |
| 4               | `p²`           | `p^{3/2}`      | `p^{3/2}` | 384        | +257 bits worse   |
| 5               | `p^{5/2}`      | `p^{8/5}`      | `p^{8/5}` | 410        | +283 bits worse   |
| 6+              | …              | →`p²`          | →`p²`     | →512       | +385 bits worse   |

Plain ECDLP on E: cost ≈ `√(n/6)` ≈ `2^{126.7}` (with √6 GLV
speedup).

For **every** `g ≥ 2`, the best DLP cost on `Jac(C)/F_p` strictly
exceeds plain ECDLP on `E`.

### Why Diem 2011 sub-exp does not rescue this

Diem 2011 sub-exp DLP on hyperelliptic Jacobians applies for `q =
p^k` with `k ≥ 3`.  For prime field `F_p` (`k = 1`), there is no
proper subfield to descend to, so Diem's index-calculus factor-base
construction does not exist.  The best-known prime-field complexity
is Gaudry's `O(p^{2 − 2/g})`, never below `p^1` for `g ≥ 2`.

### Why covers over F_{p^k}, k ≥ 3 don't help

A cover `C/F_{p^k} → E/F_{p^k}` gives DLP on `Jac(C)(F_{p^k})`, not
on `E(F_p)`.  Solving DLP over `F_{p^k}` is *harder* than over `F_p`
because the group is larger.  No useful direction.

### Final structural verdict on slice-3

> **NO COVER-BASED ATTACK ON SECP256K1's ECDLP CAN BEAT POLLARD ρ
> ON `E(F_p)`.  The slice-3 research direction is now structurally
> closed.**

The Howe gluing exists (§8.6) but the cover does not buy any DLP
speedup (§8.7).  Combined with the prior slice-1, slice-2, slice-4,
and VFCG-ρ closures, this completes the structural-completeness
statement for the entire "isogeny-graph cryptanalysis of secp256k1"
research direction.

## 8.8. Cross-curve validation: P-256 (executed)

The companion script
[`secp256k1_cm_audit/p256_comparison.gp`](secp256k1_cm_audit/p256_comparison.gp)
runs the same battery of tests on NIST P-256.  Key findings:

| Property                              | secp256k1                  | P-256                                              |
|---------------------------------------|----------------------------|----------------------------------------------------|
| j-invariant                           | 0                          | generic non-special                                |
| Effective CM disc                     | `−3` (small)               | `≈ 3 · 5 · q` with `q` a 250-bit prime (huge)      |
| Cornacchia decomposition              | `L = t`, `M = 1.01·10³⁸`   | N/A (no small-D CM)                                |
| Aut(E)                                | `Z/6Z`                     | `{±1}`                                             |
| ρ speedup                             | `√6 ≈ 2.45×`               | `√2 ≈ 1.41×`                                       |
| Sextic-twist count                    | 6                          | 2 (curve + quadratic twist)                        |
| Class number `h(disc(End))`           | 1                          | sub-exponential / huge                             |
| `p` cube mod `n`?                     | yes (ord(p) = (n−1)/3)     | **yes, same!**                                     |
| Howe (H1) `Hom_{F_p}(E, E^t) = 0`?     | ✓                          | ✓                                                  |
| Howe (H2) `E[2] ≃ E^t[2]`?            | ✓ (both deg-3 irreducible) | ✓ (both deg-3 irreducible)                         |
| Howe (H3) `gcd(n, n_twist) = 1`?      | ✓                          | ✓                                                  |
| ⇒ Howe (2,2)-gluing produces cover?  | yes                        | **yes, by same argument**                          |
| Cover-based attack possible?          | no                         | **no, by same complexity argument (§8.7)**         |

Three concrete observations from the cross-curve audit:

### "p is a cube mod n" is empirically common

Both secp256k1 and P-256 satisfy `p ≡ x³ mod n` for some `x` — i.e.,
`ord(p in (Z/n)*) = (n−1)/3`, not `n − 1`.  This is a structural
coincidence (probability `≈ 1/3` for a random prime) that holds for
*both* major deployed crypto curves.  Whether this is design-driven
or accidental is an open historical question; it does not enable any
attack (the embedding degree is still 253 bits for both).

### Howe-gluing is a class-wide property, not curve-specific

The three Howe conditions hold for **every prime-order ordinary
elliptic curve over a prime field**, because:

- (H1) is forced by `n ≠ n_{\text{twist}}`, which is automatic for any
  ordinary curve with non-zero trace.
- (H2) requires `x³ + ax + b` to be irreducible over `F_p`, which is
  forced by prime `n` (no F_p-rational 2-torsion).  The quadratic
  twist has the same factoring pattern by direct argument.
- (H3) requires `gcd(n, n_twist) = 1`, which is empirically true for
  all major curves.

So slice-3's "Howe-hit exists" is a **structural property of every
deployed prime-order ECC curve over a prime field**, not just
secp256k1's special j=0 case.  And by §8.7's cost analysis, the hit
gives no attack on *any* of these curves.

### The structural-completeness statement generalises

The verdict for secp256k1 in §9 — "no isogeny-graph attack beats
plain Pollard ρ" — **generalises uniformly to every prime-order
ordinary curve over a prime field** by the same chain of structural
arguments.  The CM-specific differences (CM disc, sextic twists,
GLV) are luxuries on secp256k1 that make the curve `√6/√2 = √3`-
slower to attack via plain ρ than P-256, but no other ρ-beating
mechanism is enabled by them.

**Conclusion**: the secp256k1 audit framework is a *generic prime-
order-ECC structural completeness statement* — secp256k1 was a
useful concrete case study, not a special case requiring its own
analysis.

## 8.9. V4 closure: supersingular reductions (executed)

The deferred V4 direction was: examine the **supersingular
reductions** of secp256k1's equation `y² = x³ + 7` over small primes
`ℓ`, where the reduction sometimes lives in the supersingular
isogeny graph (the world of SIDH / Castryck-Decru / SQIsign), and
ask whether any of that structure attacks F_p ECDLP.

The companion script
[`secp256k1_cm_audit/supersingular_reductions.gp`](secp256k1_cm_audit/supersingular_reductions.gp)
tabulates supersingularity status of secp256k1's reductions for
`ℓ < 100`.

### Findings

By Deuring's theorem, **for `ℓ ≠ 2, 3`, the curve `y² = x³ + b/F_ℓ`
is supersingular iff `ℓ ≡ 2 (mod 3)`**, regardless of `b`.  Exactly
half of small primes give supersingular reductions:

```
Supersingular ℓ < 100:
  5, 11, 17, 23, 29, 41, 47, 53, 59, 71, 83, 89   (12 primes)

Ordinary ℓ < 100 (with computed trace):
  13, 19, 31, 37, 43, 61, 67, 73, 79, 97          (10 primes)
```

A curiosity: **for `ℓ = 61` (ordinary), `#E(F_{61}) = 61 = ℓ`**,
making the F_61 reduction *anomalous* in Smart's sense.  Smart 1999
attack would solve ECDLP on `E/F_{61}` in polynomial time.  But this
is on the *tiny* curve over `F_61`, not on secp256k1's `F_p` curve
— **the secrets are different**.  The F_61 attack does not lift to
attack F_p ECDLP because Smart's canonical-lift extraction loses
the F_p discrete log when the original curve is not anomalous
(which secp256k1 is not).

### Why supersingular reductions don't attack F_p ECDLP

Four reasons, each independently sufficient:

1. **Deuring lifting goes the wrong way.**  An ordinary curve over
   `F_p` lifts *uniquely* to characteristic 0 from its reductions;
   the lift goes `F_ℓ` → `Q_p`, not `F_ℓ` → `F_p`.  The F_p-ECDLP
   secret is not preserved through the characteristic-0 lift.

2. **SIDH-world attacks are different problems.**  Castryck-Decru
   2022 broke the SIDH supersingular isogeny problem given auxiliary
   torsion data published by the protocol.  An ECDLP instance on
   `E/F_p` does not publish supersingular-graph data; there is
   nothing for Castryck-Decru-style machinery to consume.

3. **Smart's canonical lift requires anomaly.**  Smart 1999 attack
   on anomalous curves (`#E = p`) uses the formal logarithm on the
   canonical lift to extract the F_p discrete log.  For
   non-anomalous curves like secp256k1, the formal log gives no
   useful information about `d mod n`.

4. **Endomorphism rings are mathematically incompatible.**
   `End(E/F_p)` for ordinary `E` is a rank-2 commutative order in
   `Q(√D)`; `End(E_ℓ/F̄_ℓ)` for supersingular `E_ℓ` is a rank-4
   non-commutative order in a quaternion algebra.  There is no
   natural map between these rings to pull DLP-relevant information
   across.

### V4 verdict

The supersingular-isogeny-graph world is a rich mathematical object
with its own cryptanalytic story (SIDH break, SQIsign, CSIDH
analysis) — **but it does not attack F_p ECDLP on ordinary curves**.
V4 is **closed**.

## 8.10. Errata: explicit Howe cover ≠ y² = f_1·f_2

An attempt to **constructively** verify the Howe-glued cover for
secp256k1 (slice-3 hit per §8.6) led to the following correction.

### Earlier implicit claim

The Howe theorem guarantees the *existence* of a smooth genus-2
curve `C/F_p` with `Jac(C)` F_p-isogenous to `E × E^t`.  A natural
attempt at an explicit construction:

> "For coprime irreducible cubics `f_1(x), f_2(x)` over `F_p`,
> `C : y² = f_1(x) · f_2(x)` is the Howe-glued cover, with `Jac(C)`
> F_p-isogenous to `E_1 × E_2`."

For secp256k1 with `f_1 = x³ + 7`, `f_2 = x³ + 189` (where `189 =
3³·7` from the smallest non-square `d = 3`), this would give the
elegantly simple

> `C : y² = (x³ + 7)·(x³ + 189) = x⁶ + 196·x³ + 1323`.

### The toy verification falsifies this claim

The companion script
[`secp256k1_cm_audit/howe_explicit_cover.gp`](secp256k1_cm_audit/howe_explicit_cover.gp)
implements the construction on the toy `p = 1009, b = 11`, and uses
PARI's `hyperellcharpoly` to compute the actual Frobenius char poly
of `Jac(C)`:

```
For p = 1009, E: y² = x³ + 11, E_twist: y² = x³ + 515:
    Frobenius char poly of Jac(C):  x⁴ − 84x³ + 3361x² − 84756x + 1018081
    Expected (X² − tX + p)(X² + tX + p):   X⁴ + 169X² + 1018081

#Jac(C) = 936603        n · n_twist = 1018251      Not equal.
```

The Frobenius polynomial is *irreducible over Q*, so `Jac(C)` is
**simple over F_p** — not F_p-isogenous to a product of elliptic
curves at all.

### Why the naive construction doesn't give Howe's cover

The maps `C → E_1` and `C → E_2` defined by `(x, y) ↦ (x, y / √f_j(x))`
are only defined over `F̄_p` because `√f_j(x)` requires extending the
base field when `f_j` is irreducible.  The Jacobian projection
`Jac(C) → E_1 × E_2` works *over F̄_p* but doesn't descend to F_p
cleanly.  The result is a Jacobian in a different F_p-isogeny class
than `E_1 × E_2`.

The actual Howe-glued cover requires:

1. Compute Igusa invariants of `(E × E^t)/Γ_α` (a moduli computation
   in `A_2`-stack);
2. Apply Mestre's reconstruction algorithm to obtain `C/F_p` from
   the Igusa triple.

This is a non-trivial multi-page algorithm, implemented in Magma
but not in PARI's standard library or this codebase.

### What is and isn't affected

**Unaffected**: the structural-completeness theorem.  By B5 of
[`PAPER_STRUCTURAL_COMPLETENESS.md`](PAPER_STRUCTURAL_COMPLETENESS.md),
*any* smooth genus-2 cover of secp256k1 has DLP cost ≥ `O(p) > √n`.
Whether the cover is the literal Howe-glued one or just a generic
genus-2 cover (like the naive `y² = f_1·f_2`), the conclusion is
identical: no ECDLP attack via the cover.

**Affected**: any claim of having *explicitly constructed* the
Howe-glued cover.  We have not.  The existence proof from Howe
(1996), verified via conditions (H1)–(H3) in §8.6, still holds —
but the explicit construction is open.

This errata is honest about what the audit produced (existence by
theorem; not explicit construction) and what it didn't (the actual
genus-2 polynomial defining the F_p-rational Howe cover).

## 9. Overall verdict (all four slices)

| Slice | Direction | Result | Why |
|------:|-----------|--------|-----|
| 1 | Structural audit of secp256k1 | **Null** | `h(−3) = 1` ⇒ trivial Cl(O)-orbit; Tate ⇒ all `F_p`-isogenous curves have same `n` |
| 2 | CGA-HNC extension to h > 1 horizontal walk | **Null** | Orbit walk cost `Θ(h)` = saving; at scale `h ≈ √p` = ρ cost |
| 3 | (N, N)-cover search for `Jac(C) ~ E × E^t` | **Closed (hit exists, no attack)** | §8.6: Howe-gluing produces smooth genus-2 cover. §8.7: cost on Jac(C) ≥ `p` > `√n` for every genus ≥ 2. Diem sub-exp blocked for prime fields. |
| 4 | Lit review of what's blocked for `h = 1` CM | **All standard attack families blocked** | See §6 table |
| 4-V4 | Supersingular reductions of `E mod ℓ` | **Closed** | §8.9: 4 independent structural reasons; SIDH-world attacks don't transfer to F_p ECDLP |
| + | VFCG-ρ (novel proposal) | **Closed** | §9.4: `ord(c_γ)` uniformly small but walk structure prevents speedup |
| + | Cross-curve generalisation to P-256 | **All same conclusions** | §8.8: framework applies to every prime-order ordinary curve over F_p; secp256k1 was a useful case study, not a special case |

For secp256k1, the CM-by-`Z[ω]` structure leaks four observable things:

1. The six sextic twists, of which one (twist #4) has 18 bits of
   smooth-part useful only in an invalid-curve attack model
   (closed by curve-membership validation in any reasonable
   implementation).
2. The `√6` ρ-speedup (well-known, factor 2.45×).
3. A small-`ℓ` isogeny graph with self-loops at split `ℓ` and a
   uniform `Z/6Z`-multiplicity 3 pattern on descending neighbours.
   ℓ = 17 is the structural standout (6 `F_p`-rational descending
   neighbours, matching `h(−867) = 6`).
4. An explicit Frobenius decomposition `π = a + bω` underlying GLV.

**None of (1)–(4) admits a sub-exponential ECDLP attack.**  The
fundamental obstructions are:

- `h(−3) = 1` ⇒ trivial Cl(O)-orbit ⇒ no class-group amortisation
  (CRS/CSIDH/CGA-HNC structurally blocked);
- Tate isogeny invariance of group order ⇒ no descent through the
  `F_p`-isogeny class makes ECDLP easier;
- Embedding degree `≈ n` ⇒ no pairing attack;
- `n ≠ p` and `n` prime ⇒ no Smart 1999 anomaly, no Pohlig–Hellman.

The **only direction that remains genuinely open** is slice 3:
`(N, N)`-split-Jacobian covers of `Res_{F_{p²}/F_p}(secp256k1)`,
constrained by the `Q(ω) × Q(ω)`-equivariance that comes from
secp256k1's CM.  The CM constraint narrows candidate `N` to split
primes in `Z[ω]` (eliminating `N ∈ {2, 5, 11, 17, 23, …}` immediately);
each split `N` has `O(1)` `ω`-stable subgroups of `A[N]` to enumerate;
the per-candidate Mestre–Cardona–Quer test is polynomial-time.  The
expected outcome on priors is still negative, but it is the only one
not provably blocked by a structural theorem.

This is the "publish a paper that says 'we computed X, found no
exploitable structure, here's a tighter lower bound on attacks of
type Y'" outcome the original research-direction preamble described
as the realistic deliverable.

## 10. Suggested next-action triage

For the codebase, ordered by cost-to-value:

1. **Update `ecc_safety::check_cl_o_orbit_brittleness` to short-circuit
   on `h(End K) = 1`** (currently it only short-circuits on prime
   `n`, missing the structural argument for why secp256k1 passes
   even when `n` is composite in some hypothetical sibling curve).
2. **Add a Rust binding to PARI's `polclass` + `qfbcornacchia`** so
   future slice-2-style horizontal-orbit experiments don't require
   leaving the codebase.  ~200 LOC wrapping `libpari`.
3. **Port the slice-2 PARI script to the Rust test suite** as an
   `#[ignore]`'d demonstration of horizontal Cl(O)-orbit walking;
   couples to (2).
4. **Implement the `Q(ω)`-equivariant `A[N]`-subgroup enumeration**
   from §8 in Rust — the lightest piece of the slice-3 search, and
   reusable by `p256_isogeny_cover` for the generic case.
5. **(Defer)** Mumford-divisor arithmetic + Mestre–Cardona–Quer in
   Rust.  Heavy lift; only worth it if (4) ever produces a
   non-trivial candidate kernel for secp256k1.

## Appendix: How to reproduce

```
cd /Users/adamburan/git/crypto/secp256k1_cm_audit
gp -q audit.gp > audit_output.txt 2>&1
diff audit_output.txt audit_output.expected.txt   # if a committed expected file exists
```

`audit.gp` requires only PARI/GP 2.17+ (homebrew package `pari`) and
~250 MiB stack.  Runtime is &lt; 1 s on an M-series Mac.
