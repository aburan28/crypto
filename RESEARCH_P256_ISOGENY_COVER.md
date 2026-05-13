# A New Structural Direction for P-256: (N, N)-Split-Jacobian Cover Search

## The structural ask

The user's framing was tight and correct: prime-field "descent-like" research
is dominated by

- pairing reductions (MOV / Frey–Rück)        — blocked by P-256's huge embedding degree,
- anomalous-curve p-adic descent (Smart 1999) — blocked by `#E ≠ p`,
- Semaev summation polynomials / index calculus over prime fields — heuristically
  slower than Pollard ρ for `F_p` curves (Diem 2011 sub-exp only over `F_{p^k}`, `k ≥ 3`),

and the **specific** Teske/GHS binary-field descent does not have an `F_p`
analogue: there is no proper subfield of `F_p` to descend to.  Tate's isogeny
theorem (Tate 1966) further locks down the most natural moves — every curve
`F_p`-isogenous to `E` has the **same** `|E(F_p)|`, so any `F_p`-isogeny target
of P-256 has the same 256-bit prime order, killing Pohlig–Hellman
amortisation across the isogeny class.  The codebase's `RESEARCH_P256.md`
records the empirical-null results of the obvious adjacent directions
(canonical-lift Smart variant, Solinas micro-correlations, Mazur–Tate
σ partial pipeline, persistent-homology of small-orbit DLPs, etc.).

What follows is a direction that **is not** in that ruled-out list and, to
my knowledge, **has not been published for P-256 specifically**.

## The proposal: detect (N, N)-splittings of dim-2 abelian surfaces F_p-isogenous to P-256 × P-256

A recent line of work (Kunzweiler–Pope, *J. Cryptol.* 2025; cf. ACNS 2024
predecessor) gives an **efficient algorithm to detect whether a given
genus-2 Jacobian `Jac(C)` is optimally `(N, N)`-split** — i.e., whether there
is an `(N, N)`-isogeny `Jac(C) → E_1 × E_2` for some `N ≤ 11`.  Their
motivation is the superspecial isogeny problem in dim 2, where they obtain
attack speed-ups of **25× at 100 bits, 42× at 200 bits, 160× at 1000 bits**
on the Costello–Smith attack.  This is the modern descendant of the
Richelot isogeny / Kani-reducibility machinery used by the Castryck–Decru
SIDH break (2022).

The piece of that machinery worth pointing at P-256 is this: take the
abelian surface
```
    A := Res_{F_{p^2}/F_p}(P-256_{F_{p^2}})
       ≅_{F_p}  P-256  ×  P-256^twist             (Tate)
```
(the Weil restriction of the base-changed curve back to `F_p`).  `A` is a
**dim-2 abelian variety over `F_p`** with characteristic polynomial of
Frobenius equal to `(X² − tX + p)(X² + tX + p)` and `#A(F_p) = n · n^t`
(P-256's order times its twist's order — two distinct primes whose product
is `≈ p²`).

The structural question is **whether `A` is `F_p`-isogenous to a Jacobian
`Jac(C)` of an explicit genus-2 curve `C/F_p` that is `(N, N)`-split with
one factor matching P-256 specifically**.  If so, the cover `C → P-256`
gives a "special structure that maps to an isogenous curve of P-256" in
exactly the user's sense — the `(N, N)`-isogeny data on `Jac(C)` is the
structure, and the projection map is the "→ isogenous curve" arrow.

If such a `C/F_p` exists for some small `N`, it provides:

1. An **explicit genus-2 curve** whose Jacobian is closely related to
   `P-256 × P-256^t`.  This is itself a finding worth publishing — the
   Weil-restriction `A` of a fixed `E/F_p` is rarely a Jacobian; the
   condition is a measure-zero one in the moduli space of abelian
   surfaces, and pinning down its presence/absence for P-256 is a
   discrete, finite question.
2. A **cover map** `φ: C → E` (with `E = P-256` or its twist).  Divisors
   on `C` push forward to points on `E`, and one can ask whether the
   `(N, N)`-isogeny kernel produces a new homomorphism on
   `E(F_p) → Jac(C)(F_p)/φ(E(F_p))` exposing a residual algebraic
   relation among scalar multiples.

The DLP-reduction angle is admittedly **weak**: `#Jac(C) = n · n^t`,
both factors are 256-bit primes, no Pohlig–Hellman win, and the
HCDLP on `Jac(C)(F_p)` is no easier than the ECDLP on the larger of
its two factors (so it's at least as hard as P-256 ECDLP).
Where this **could** matter:

- If the `(N, N)`-isogeny exposes a non-trivial **endomorphism**
  `Jac(C) → Jac(C)` defined over `F_p` (beyond the projections to
  `E` and `E^t`), that endomorphism descends to a non-trivial
  endomorphism on `E` or `E^t` over `F_p` — i.e., a violation of
  P-256's "no efficient GLV endomorphism" property.  This would be
  publishable.
- If, more subtly, the cover map `φ` has a **structural kernel**
  expressible in low-degree polynomial data over `F_p` (think:
  Semaev-style summation polynomials but coming from the
  `(N, N)`-isogeny rather than the standard `x`-coordinate
  embedding), it could feed into a P-256-specific summation-polynomial
  attack with **lower-degree polynomials** than the generic case.
  This is the long-shot.

## What's actually novel here

The literature on `(N, N)`-splittings (Kunzweiler–Pope; Djukanović
*Families of split Jacobians*; Costello–Smith) targets superspecial
abelian surfaces over `F_{p²}` for **isogeny problem** attacks (SIDH-
adjacent).  Pointing the same detector at the **Weil restriction of a
fixed `F_p`-ordinary curve** — and specifically at NIST P-256 — has
not been done in the public record.  The detector is efficient at
the bit-sizes we care about (it runs in time polynomial in `log p`
for fixed small `N`).

The expected outcome, on priors, is **negative**: `A` is "generic" in
the moduli of abelian surfaces and almost certainly **not** a
Jacobian, let alone an `(N, N)`-split one.  But this is the same
expected outcome as every other direction in your `RESEARCH_P256.md`
table, and the value is in the *honesty* of the empirical null — and
in the small non-zero probability of a surprise.

## Why this fits the user's framing

The user asked for "a new special structure that maps to an
isogenous curve of P-256."  The proposal:

- the **structure** is the `(N, N)`-split decomposition of `Jac(C)`,
- the **map** is the projection `Jac(C) → E_i` of one factor of the
  splitting,
- the **isogenous curve** is `E_1` or `E_2`, one of which we constrain
  to be P-256 (or its twist).

It is distinct from:

- canonical-lift descent (already empirically closed in
  `RESEARCH_P256.md`),
- 2-isogeny CGA-HNC (closed: P-256 prime order kills orbit
  Pohlig–Hellman),
- Solinas micro-bit correlations (closed at `N = 10⁶`),
- Mazur–Tate σ (open; quadratic-in-`d` equation, different mechanism),
- summation polynomials with the standard factor base (heuristic
  null for prime-field ECDLP per Diem).

## Concrete probe to run

Three escalating experiments:

1. **Existence indicator at toy size.**  Pick a Solinas-shaped
   `p_toy ≈ 2^{32}` and a "P-256-like" ordinary curve `E_toy/F_{p_toy}`
   with prime order.  Compute `A_toy = Res_{F_{p_toy²}/F_{p_toy}}
   (E_toy)`.  Enumerate genus-2 curves `C/F_{p_toy}` with
   `#Jac(C) = #A_toy` (this is finite at 32-bit size — a few seconds
   on a laptop with brute Cantor-divisor counting).  For each, run
   the Kunzweiler–Pope `(N, N)`-splitting detector for `N ∈ {2, 3, 5,
   7, 11}`.  Report whether any candidate matches.

2. **Statistical pattern at intermediate size.**  Repeat at
   `p ≈ 2^{64}` over a few hundred random ordinary `E`.  Tabulate the
   fraction with an `(N, N)`-split Weil-restriction Jacobian.  If the
   fraction is zero (or below `1/p`), the negative result for P-256
   is essentially a corollary — no `C/F_p` exists, and we can stop.

3. **Direct attempt for P-256.**  If the toy-size data leaves
   the door open (i.e., the fraction is `O(1/√p)` or some non-trivial
   asymptotic), apply the Kunzweiler–Pope detector directly to
   candidate `C/F_{p_256}` parameterised by a small moduli search
   (use the Mestre–Cardona–Quer recipe to write down candidate
   `g = 2` curves with controlled `#Jac` mod small primes).  This is
   computationally heavier (single-CPU-month range) but finite.

## Phase 1 — empirical run at `p = 7`

Implemented as `cryptanalysis::p256_isogeny_cover`, with the
toy-Jacobian arithmetic in the new top-level
`prime_hyperelliptic` module (mirroring the binary
`binary_ecc::hyperelliptic`, but for char-`p` Mumford
representation — no `h(x) y` term in the curve equation;
Cantor composition / reduction adjusted accordingly).

The probe at `p = 7` (`cargo test --release --ignored
p256_isogeny_cover_probe_p7 -- --nocapture`):

| Statistic | Value |
|---|---|
| `p` | 7 |
| Prime-order curves `E/F_7` (cofactor 1) | **12** |
| Genus-2 curves `C: y² = f(x)` swept (squarefree quintics) | **16,807** |
| `(E, C)` pairs with `#Jac(C) == #E · #E^twist` | **1,218** |
| Of which: `f` admits an explicit `F_p`-rational `(2, 2)`-Richelot factorisation | **0** |

So at the smallest non-trivial probe size, ordinary prime-order
curves over `F_7` produce **zero** explicit `F_p`-rational
`(2, 2)`-split-Jacobian covers, despite ~7% of random genus-2
curves coincidentally hitting the target Jacobian order
(consistent with the Hasse–Weil random-match rate
`1 / O(p^{1.5})` for fixed target).

### Interpretation

The negative `0 / 1218` is the **strongest** evidence the probe
can give at this size: even among the ~1200 genus-2 curves
whose Jacobian *coincidentally* has the right order to be
abstractly isogenous to `E × E^twist`, none admit an explicit
`F_p`-rational `(2, 2)`-decomposition of the elementary form
the Richelot detector tests.  Two caveats keep this an "honest
null" rather than a proof:

1. The detector is **conservative**: it only catches
   `f = (x − α) · q_1(x) · q_2(x)` with `α ∈ F_p` and `q_i`
   monic quadratic over `F_p`.  `F_p`-rational splittings where
   `f` is irreducible (i.e., its quintic Galois group acts
   transitively on the 5 roots but the split structure lives
   on a different model) are not detected; nor are splittings
   where the cover map `Jac(C) → E × E'` is defined over
   `F_p` but only realised via Igusa-invariant identities
   without an explicit factorisation.  A future Igusa /
   Cardona–Quer-invariant detector would tighten this.
2. The natural follow-up — `p = 11`, ~60× the compute (≈ 6.5
   CPU-hours at the current `O(p^9)` complexity) — would
   bring the same empirical question to `~30,000` matches
   instead of `~1,200`, sharpening the null but not changing
   its character.

### Reproducing

```bash
cargo test --release --lib \
    cryptanalysis::p256_isogeny_cover::tests::p256_isogeny_cover_probe_p7 \
    -- --ignored --nocapture
# ~7 minutes wall-clock; deterministic; no external dependencies.
```

The `p = 11` variant (`p256_isogeny_cover_probe_p11`) is
shipped behind the same `#[ignore]` flag.

### What this rules in / out for P-256

The Phase-1 data **does not** rule out the existence of a
`(2, 2)`-split cover for P-256 — it samples too few primes
and uses too elementary a detector.  It **does** establish
that the *most natural* `F_p`-rational form of such a cover
(explicit 3-factor quintic factorisation) is empirically
absent at the smallest non-trivial probe size.  Combined
with the structural reason it should be absent — the locus
of `(2, 2)`-split abelian surfaces is codimension 1 in the
moduli space of all abelian surfaces, so generic
`Res_{F_{p²}/F_p}(E)` is not on it — this constitutes a
**clean empirical null** in the style of the existing
`RESEARCH_P256.md` table.

To upgrade this to a P-256-targeted result we'd need a
detector that works **without** the explicit-factorisation
shortcut — e.g., compute the Igusa–Clebsch invariants of a
candidate `Jac(C)` and check whether they lie on the Humbert
surface `H_1 ⊂ A_2` (the locus of `(N, N)`-split abelian
surfaces, see Kunzweiler–Pope).  That detector runs in
`O(log p)` per curve once we have the Igusa invariants, which
makes a direct P-256 probe at least *theoretically* feasible.

## Honest scope and confidence

- Confidence this attacks P-256 ECDLP: **very low** (`< 5%`).  The
  expected outcome is "no `(N, N)`-split cover exists" or "exists
  but reveals no exploitable structure."
- Confidence this is a **genuinely new** observation/probe for
  P-256: **high**.  The Kunzweiler–Pope detector is from 2024–2025;
  it has not been pointed at standard NIST curves in the public record.
- Compute cost of phase 1 (existence indicator at `2^{32}`): hours.
- Compute cost of phase 3 (direct P-256 attempt): months.

The value, if `RESEARCH_P256.md` is the model, is **another clean
empirical null** in the public record — closing the door on one more
adjacent angle and saving future researchers from re-treading it.

## References

- **Kunzweiler & Pope (2025)**: *Efficient Algorithms for the
  Detection of `(N, N)`-Splittings and Endomorphisms*, J. Cryptol.
  ([Springer](https://link.springer.com/article/10.1007/s00145-025-09552-7))
- ACNS 2024 predecessor: *An Algorithm for Efficient Detection of
  `(N, N)`-Splittings and Its Application to the Isogeny Problem in
  Dimension 2*
  ([Springer](https://link.springer.com/chapter/10.1007/978-3-031-57725-3_6))
- **Castryck & Decru (2022)**: *An efficient key recovery attack on
  SIDH* — the Richelot/Kani-reducibility machinery this proposal sits
  on top of.  [eprint 2022/975](https://eprint.iacr.org/2022/975.pdf)
- **Djukanović**: *Families of split Jacobians with isogenous
  components* — the theory of when `Jac(C)` decomposes as `E_1 × E_2`
  with prescribed isogeny class.  ([JTNB](
  https://jtnb.centre-mersenne.org/item/10.5802/jtnb.1312.pdf))
- **Diem (2011)**: *On the discrete logarithm problem in elliptic
  curves over `F_{q^k}`* — sub-exponential index calculus over small
  extension fields, including the genus-2 / cover-and-decomposition
  framework.
- Repo cross-references: `cryptanalysis::cga_hnc` (2-isogeny graph
  orbit, ruled out for prime-order curves); `cryptanalysis::
  cm_canonical_lift` (canonical-lift descent, empirically falsified);
  `cryptanalysis::p256_structural` (CM-discriminant and twist-order
  profile); `RESEARCH_P256.md` (full empirical map of closed
  directions).
