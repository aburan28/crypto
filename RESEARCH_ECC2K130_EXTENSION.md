# Raising ECC2K-130 to an extension field

**Experiment:** `scripts/ecc2k130_extension_field_boundary.py`
**Frozen artefact:** `experiments/ecc2k130_extension_field_boundary.json`
**Related:** `RESEARCH_KOBLITZ_INDEX_CALCULUS.md` (the invariant-factor-base
pipeline this note bounds), `RESEARCH_QUASI_SUBFIELD.md` (the same barrier
approached from the factor-base side), `research/ecdlp_autolab/phase6_weil_restriction.md`
(the prime-field version of the question), `RESEARCH_ISOGENY_CLASS_SEARCH.md`
(the other "change the object, not the algorithm" lever).

The question.  `E(F_2^131) ≤ E(F_2^(131e))` for every `e ≥ 1`, so an attacker
is free to work in any of those larger fields, and the larger ones have
structure the challenge field conspicuously lacks — composite extension
degree, proper subfields, Weil restrictions of small dimension.  Does moving
up buy an attack cheaper than Pollard rho on the original group?

**Bottom line.  No, and not because a search came up empty.**  All five
boundaries below are derived, and between them they leave no `e` any room.
The reason is a single arithmetic fact — **131 is prime** — which forces a
dichotomy: in `F_2^(131e)`, every subfield and every Frobenius-stable
subspace is either *too large to enumerate* or *too small to generate*, with
nothing in between, at every `e`.

- **A.** Base change moves neither the target group nor its automorphisms, so
  the reference boundary is `2^60.81` at every `e` (§1).
- **B.** Every Frobenius-stable `F_2`-subspace of `F_2^(131e)` has dimension
  `≤ e` and lies inside `F_2^e`, or has dimension `≥ 130`.  Proved in §2 and
  checked exhaustively for `e = 1 … 16, 131, 262`.  There is no third case.
- **C.** The large horn (`131 | d`) carries a factor base of `2^131` points;
  its published complexity bottoms out at `2^131` operations, `2^70.19×` the
  rho reference, and *rises* with `e` (§3).
- **D.** The small horn is not expensive, it is **empty**: the factor base is
  contained in the subgroup `E(F_2^e)` of order `≈ 2^e`, so no target in `⟨G⟩`
  has a decomposition at any cost (§4).
- **E.** The one route that really does trade a curve group for an easier
  field — MOV/Frey–Rück transfer — needs the embedding degree, which is
  `> 10^7`: the target field has more than `1.31 × 10^9` bits (§5).

What *is* worth taking from the idea is one register down: not a bigger
**field**, a bigger **ring**.  `F_2^131` embeds in `F_2[x]/(x^263 − 1)`, where
Frobenius is a cyclic shift.  That is the type-II optimal normal basis the
`ecc2k130/` client already runs on (§6) — an **engineering** gain by §3 of
`AGENTS.md` — already banked, and it does not touch the exponent.

## 0. The boundary, stated before anything is measured

`K_0 : y² + xy = x³ + 1` over `F_2`, used over `F_2^131`; trace `t = −1`.

```
#E(F_2^131) = 4 · r,   r = 680564733841876926932320129493409985129   (prime, 2^129.0000)
```

The reference is Pollard rho on `⟨G⟩` with the `⟨−1⟩ × ⟨π⟩` speed-up, the
same accounting the rest of the repository uses:

| quantity | value |
|---|---|
| automorphisms on `⟨G⟩` | `2 · 131 = 262` |
| rho, plain | `2^64.8257` |
| **rho, reference** | **`2^60.8090`**, `S = 0.0774` |

Matching `experiments/isogeny_class_search.json` (`log2_rho_reference`
`60.8090`) to the digit.  Every figure below is a ratio to that number.

**Falsification target.**  This thread is a success if some `e ≥ 2` admits a
complete attack — relation collection, linear algebra and descent, all phases
priced per §5 of `AGENTS.md` — costing under `2^60.81` group operations.
It is abandoned if the two horns can be shown to exhaust the options, which
is what §2 does.

## 1. Boundary A — what base change does not move

`⟨G⟩` is a cyclic group of order `r`.  It is the *same* cyclic group inside
`E(F_2^131)` and inside `E(F_2^(131e))`: base change adds points, not
structure on the ones already there, and the discrete logarithm is asked in
`⟨G⟩` either way.  So rho's cost is untouched — the reference is the same
`2^60.81` at every `e`.

Nor does the rho *speed-up* improve.  Its size is the order of the group of
automorphisms acting on `⟨G⟩`, and:

- `t = −1` is odd, so `E` is **ordinary**, and for an ordinary curve
  `End_{F_q}(E) = End_{F̄_q}(E)` (Waterhouse).  Every endomorphism is already
  defined over `F_2`; a larger field supplies no new one.
- `π` acts on `⟨G⟩` as `λ` with `λ² + λ + 2 ≡ 0 (mod r)` and `λ^131 = 1`.
  That holds in every ambient field, so the orbit length stays `131` and the
  usable automorphism group stays `262`.

This is the part of the question that has a one-line answer: **the thing you
would be attacking, and the best generic attack on it, are invariant.**  Only
the *non-generic* attacks can possibly notice the extension, which is what
B–E are about.

## 2. Boundary B — the dichotomy, and why 131 being prime is the whole story

Every factor base this repository's pipeline can build (GGMP §4, and the
quasi-subfield generalisation) is the point set above a Frobenius-stable
`F_2`-subspace `V ⊆ F_2^N`.  Those are exactly the kernels `ker g(σ)` for
divisors `g` of `t^N − 1` over `F_2`, with `dim V = deg g`.

**Lemma.**  Let `N = 131 e`.  Every Frobenius-stable `F_2`-subspace `V` of
`F_2^N` satisfies either `V ⊆ F_2^e` (hence `dim V ≤ e`) or `dim V ≥ 130`.

*Proof.*  Write `e = 2^v e'` with `e'` odd, so `t^N − 1 = (t^(131e') − 1)^(2^v)`
in characteristic 2.  The irreducible factors of `t^(131e') − 1` are indexed by
the divisors `d | 131e'`, a factor for `d` having degree `ord_d(2)`.  Split
those divisors by whether `131` divides them.  If `131 ∤ d` then `d | e'`,
since `131` is prime, and the factor is one of `t^(e') − 1`.  If `131 | d`
then reduction modulo `131` gives `ord_131(2) | ord_d(2)`, and
**2 is a primitive root mod 131**, so `ord_d(2) ≥ ord_131(2) = 130`.  Either
way every irreducible factor of `t^N − 1` divides `t^e − 1` or has degree
`≥ 130`.  A divisor `g` with
`deg g < 130` is therefore built only from factors of `t^e − 1`, with
multiplicity at most `2^v`, so `g | (t^(e') − 1)^(2^v) = t^e − 1`, and
`ker g(σ) ⊆ ker(σ^e − 1) = F_2^e`. ∎

Checked by exhaustive factor-degree enumeration for `e = 1 … 16`, plus
`e = 131` and `e = 262` — the family where `131` and `e` share a factor and
the split in the proof needs the divisibility step rather than a coprime one.
The assertion is in the script, not in the prose.

| `e` | `N = 131e` | factor degrees of `t^N − 1` | max dim inside `F_2^e` | min dim outside |
|---:|---:|---|---:|---:|
| 1 | 131 | `1, 130` | 1 | 130 |
| 2 | 262 | `1², 130²` | 2 | 130 |
| 3 | 393 | `1, 2, 130³` | 3 | 130 |
| 4 | 524 | `1⁴, 130⁴` | 4 | 130 |
| 5 | 655 | `1, 4, 130, 260²` | 5 | 130 |
| 6 | 786 | `1², 2², 130⁶` | 6 | 130 |
| 7 | 917 | `1, 3², 130, 390²` | 7 | 130 |
| 8 | 1048 | `1⁸, 130⁸` | 8 | 130 |

The same dichotomy holds for **subfields**, and there it needs no lemma:
`F_2^d ⊆ F_2^(131e)` means `d | 131e`, and `131` prime means `131 | d` or
`131 | n` for `n = N/d`.  Either the base field is at least `2^131`, or the
decomposition needs at least `131` summands.

At `e = 1` this already answers the unextended question, and it is the same
barrier `RESEARCH_QUASI_SUBFIELD.md` measures from the other side: at
`n = 131` the only subspace dimensions available are `0, 1, 130, 131`, so the
only invariant factor bases are `E(F_2)` (4 points) and one of size `2^130`.
The census in that note found no non-subfield quasi-subfield polynomial at any
prime `n ≤ 600` with `n₀ < 64`, `131` included — measured agreement with what
the divisor criterion predicts here.

**`131` being prime is not one obstacle among several.  It is the obstacle,
and the extension does not dilute it: `131 | 131e` for every `e`.**

## 3. Boundary C — the large horn costs more than the answer

Take the horn with `131 | d`: base field `F_2^d`, `q = 2^d ≥ 2^131`,
extension degree `n = 131e/d`.  This is the textbook Gaudry/Diem setting and
it *works* — it is just priced against the wrong thing.  Its published cost is
`Õ(q^(2−2/n))` for fixed `n ≥ 2`, with a constant exponential in `n`, and
`2 − 2/n ≥ 1` for every `n ≥ 2`, so the cost is at least `q ≥ 2^131`.

| `e` | `N` | `d` | `n` | factor base | `q^(2−2/n)` | ratio to rho |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 262 | 131 | 2 | `2^131` | `2^131.00` | `2^70.19` |
| 3 | 393 | 131 | 3 | `2^131` | `2^174.67` | `2^113.86` |
| 4 | 524 | 131 | 4 | `2^131` | `2^196.50` | `2^135.69` |
| 5 | 655 | 131 | 5 | `2^131` | `2^209.60` | `2^148.79` |
| 8 | 1048 | 131 | 8 | `2^131` | `2^229.25` | `2^168.44` |

Best case `e = 2`, at `2^70.19×` the reference — and the column is monotone
increasing, so no larger `e` recovers anything.  The mechanism is worth saying
plainly, because it is the general answer to the question in the title:

> **Index calculus is priced by the ambient field; rho is priced by the
> subgroup.**  Raising the field raises the first and leaves the second fixed.

ECC2K-130 is the worst possible instance for that trade: `r ≈ 2^129` already
fills its `2^131` field, so there is no slack to exploit and every step up
costs `131` bits of factor base to buy nothing.

## 4. Boundary D — the small horn has no relations, at any price

Now the horn with `V ⊆ F_2^e`.  A factor base here is not merely small, it is
**closed**:

- for an abscissa `x ∈ F_2^e` the two ordinates solve `y² + xy = x³ + b`,
  i.e. `z² + z = (x³ + b)/x²` with `y = xz` — both roots lie in `F_2^e` or
  both in `F_2^(2e)`;
- `2e ∤ 131e` (131 is odd), so `F_2^(2e) ⊄ F_2^N`, and the second case
  contributes no points of `E(F_2^N)` at all;
- hence `{P ∈ E(F_2^N) : x(P) ∈ V} ⊆ E(F_2^e)`, a **subgroup** of order
  `≈ 2^e`.

Sums of factor-base points never leave `E(F_2^e)`.  The target lies in `⟨G⟩`,
of order `2^129`, and `G ∈ E(F_2^e)` would need `λ^e = 1`, i.e. `131 | e`,
which is exactly the large horn again.  So for every `e` in the small horn:

| `e` | factor-base subgroup | order | target reachable |
|---:|---|---:|:--:|
| 1 | `E(F_2)` | 4 | no |
| 2 | `E(F_4)` | 8 | no |
| 3 | `E(F_8)` | 4 | no |
| 4 | `E(F_16)` | 16 | no |
| 5 | `E(F_32)` | 44 | no |
| 6 | `E(F_64)` | 56 | no |

Not a bad cost — **no relation exists**.  This is the subfield-curve
degeneracy that makes GGMP use invariant subspaces rather than subfields in
the first place; §2 says the extension never supplies an invariant subspace
that escapes it without landing in horn C.

A factor base that is *not* Frobenius-stable forfeits the `n` / `n²` GGMP
savings; [`RESEARCH_ECC2K130_DECOMPOSITION.md`](RESEARCH_ECC2K130_DECOMPOSITION.md)
prices that family at `e = 1` in full and finds `m·2^131`, `2^71.77×` rho, with
the `n²` saving granted for free still leaving `2^57.71×`.  For a
fixed summand count `m` the subspace dimension must be `≈ N/m = 131e/m`, so
`|F|` grows exponentially in `e` while the target stays at `2^129`.  At `e = 1`
that is exactly the pipeline of `RESEARCH_KOBLITZ_INDEX_CALCULUS.md`, whose
measured verdict at the sizes reached is `528×` to `4,000×` rho.  Extending the
field does not add an option there; it inflates `N` in the one place the cost
is exponential in `N`.

## 5. Boundary E — the transfer route, priced

The genuinely different way to "raise to an extension field" is not to move
the curve but to leave it: MOV/Frey–Rück carries the ECDLP into `F_2^(131k)^*`
where `k` is the embedding degree, and small-characteristic DLP is the one
place quasi-polynomial algorithms live.  That is a real mechanism — it is how
supersingular binary curves died — and it is the strongest form of the
question.

It needs `k = ord_r(2^131)` to be small.  Exhaustion finds no `k ≤ 10^7`:

```
embedding degree > 10^7   ⇒   target field > 1.31 × 10^9 bits
```

An element of that field does not fit in 160 MB, and the quasi-polynomial
cost `N^O(log N)` at `log₂ N ≈ 30` is past `2^900` under the most generous
reading of the constant.  Ordinary Koblitz curves have embedding degree of
size comparable to `r`; this is the expected answer and the bound above is
what this thread contributes to it, not a surprise.

## 6. Where the idea does pay — a ring, not a field

The instinct behind the question is sound at the level of *representation*.
Search for primes `p` with `ord_p(2) = 131`: there is exactly one below 4000,
`p = 263 = 2·131 + 1`, and

```
F_2[x]/(x^263 − 1) ≅ F_2 × F_2^131 × F_2^131
```

so `F_2^131` sits in a 263-bit redundant ring in which Frobenius is a cyclic
shift and squaring is a permutation.  That is the type-II optimal normal basis
— `ecc2k130/README.md` already names it ("the permuted type-II optimal normal
basis of Bailey et al."), and `ecc2k130/FROBENIUS-NETWORK.md` is the packed
backend applying those shifts as fixed bit-permutation networks.

The contrast in the same README is the cleanest possible evidence that this is
the mechanism and not a coincidence: ECC2K-95 over `F_2^97` has
`2·97 + 1 = 195` composite, no type-II ONB, and therefore runs on a
polynomial-basis backend instead.

By §3 of `AGENTS.md` this is **engineering**: `S` falls, the ratio to the
floor does not move.  It is worth roughly a constant factor in the iteration
cost and it is already banked.

## 7. The counterfactual, which is what makes the question a good one

Extension structure is not harmless in general — it is how binary-curve
ECDLPs actually fall.  Had the challenge been posed one degree lower, over
`F_2^130`, the same group size would have been broken by descent:

| write `F_2^130` as | `q` | `n` | `q^(2−2/n)` | rho on that group |
|---|---:|---:|---:|---:|
| `(F_2^65)²` | `2^65` | 2 | `2^65.00` | `2^60.31` |
| `(F_2^26)⁵` | `2^26` | 5 | **`2^41.60`** | `2^60.31` |
| `(F_2^13)^10` | `2^13` | 10 | `2^23.40` | `2^60.31` |

The `n = 5` row needs only `S_6`, which `RESEARCH_HIGHER_SEMAEV.md` already
computes, and beats rho by nineteen bits.  Certicom's choice of prime
extension degrees — 131, 163, 191, 239 — is precisely the defence, and it is
a defence the attacker cannot undo by moving up, because `131 | 131e` for
every `e`.

## 8. What this does not settle

- **Characteristic 2 and this curve.**  The lemma in §2 uses `ord_131(2) = 130`
  and the parity of `131`; it is not a statement about ECDLP in general.
- **The published complexities are taken at face value.**  §3 quotes
  `Õ(q^(2−2/n))` rather than re-deriving it; the thread's own contribution
  there is only that `2 − 2/n ≥ 1`, which needs no constant.
- **Non-invariant factor bases are bounded only by the repository's measured
  numbers**, not by a derivation.  §4's argument that extending `N` inflates
  them is a monotonicity claim, not a lower bound on the oracle cost.
- **Quasi-subfield polynomials that are not Frobenius-stable** remain outside
  §2's lemma.  `RESEARCH_QUASI_SUBFIELD.md` §3 gives measured agreement
  between the divisor criterion and exhaustive census on every cell brute
  force reaches, and §5 gives an exponent argument (`0.9487` against the
  generic `0.5`) that closes the framework regardless — but neither is a proof
  at `n = 131`.
- **The embedding degree is bounded, not computed.**  `k > 10^7` is an
  exhaustion; `k | r − 1` and factoring `r − 1` was not attempted.

## References

- P. Gaudry, *Index calculus for abelian varieties of small dimension and the
  elliptic curve discrete logarithm problem*, J. Symbolic Comput. 44 (2009).
- C. Diem, *On the discrete logarithm problem in elliptic curves*,
  Compositio Math. 147 (2011).
- P. Gaudry, F. Hess, N. Smart, *Constructive and destructive facets of Weil
  descent on elliptic curves*, J. Cryptology 15 (2002).
- S. Galbraith, R. Granger, S.-P. Merz, C. Petit, *On index calculus algorithms
  for subfield curves*, SAC 2020 (ePrint 2020/1315).
- M.-D. Huang, M. Kosters, C. Petit, S. L. Yeo, Y. Yun, *Quasi-subfield
  polynomials and the ECDLP*, J. Math. Cryptol. 14 (2020).
- D. Bailey et al., *Breaking ECC2K-130*, ePrint 2009/541.
- W. Waterhouse, *Abelian varieties over finite fields*, Ann. Sci. ENS 2 (1969).
