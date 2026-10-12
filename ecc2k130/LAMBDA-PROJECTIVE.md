# λ-projective coordinates and the invariant-hash question

Question: the lever table in [ITERATION-FUNCTION.md](ITERATION-FUNCTION.md) §1
dismisses projective coordinates in one line ("projective forms need the affine
`x` for the hash anyway"). The proposed direction is **λ-projective coordinates
plus a projective-invariant hash and distinguished-point predicate**: drop the
per-step inversion entirely by adding in λ-projective form (Oliveira, López,
Aranha, Rodríguez-Henríquez, CHES 2013: mixed addition 8M + 2S, the cheapest
inversion-free addition on a binary curve) and drive the walk with a branch
selector and a DP predicate that are functions of the point, not of the
representative `(X, L, Z)`. Does that beat the affine walk with Montgomery's
batched inversion?

Answer: **no, on the operation count alone, before the hash is priced.** At
the kernel's batch of 16 the affine step costs 4.81 products, one squaring
and one sixteenth of an inversion, 38.1 carry-less lane-instructions per
update; λ-projective mixed addition costs 8 products and two squarings,
60.5, with the invariant hash *given for free*. Its ceiling on the RTX PRO
6000 is 11.8 – 12.5 B/s at 100% of the carry-less unit, below the 14.41 B/s
the shipping affine walk measures. On the binding pipe affine wins from batch
3 upward (at batch 2 the affine walk is ALU-bound by the inversion's own
4,192 slots, and only there does the projective ceiling come out ahead), and
the kernel runs batch 16 because batch 8 and 32 both measured worse
([BATCH-TUNING.md](BATCH-TUNING.md)). The hash question is therefore moot for
the rate, and it is also closed on its own terms: the cheapest known
scale-invariant of `[X : Z]` is `X/Z`, one inversion, which is the cost
projective coordinates set out to remove; every alternative priced in §3 costs
more than that inversion; and the invariant *is* required, because a selector
on the representative parts two trails of the same point within one step
(§3.3, measured). Class: **accounting** (a row priced, nothing built, no gain
claimed), and the lever row in ITERATION-FUNCTION.md now points here.

Everything here is reproduced by
[`codegen/lambda_projective.py`](codegen/lambda_projective.py) on the
generator's own model of `GF(2^m)` and `y² + xy = x³ + 1`, and the parts that
are cheap enough for a unit test are in
[`codegen/testlambdaprojective.py`](codegen/testlambdaprojective.py)
(`make check-cli`).

## 1. Boundary

**Unit.** Carry-less (`CLMAD`) and ALU lane-instructions per scalar update,
the unit of ITERATION-FUNCTION.md §1, with the per-routine prices of
[THROUGHPUT-20B.md](THROUGHPUT-20B.md) §3 (static SASS, shipping preset,
`sm_120`):

| routine | ALU slots | clmad | source |
|---|---:|---:|---|
| field product (`mulPolynomial131`, Karatsuba + reduction) | 161 | 6 | THROUGHPUT-20B §3 |
| squaring in the polynomial basis (`spread32p` + reduce) | 83 | 6.25 | THROUGHPUT-20B §3 |
| inversion (Itoh–Tsujii: 8 products at 360 with basis conversions, five Frobenius networks) | 4,192 | 48 | THROUGHPUT-20B §3, ×16 |

Pipe rates on the 188-SM RTX PRO 6000, measured: 1.65 carry-less and 62.1
ALU lanes per SM-clock, 2.30 – 2.43 GHz (ITERATION-FUNCTION §1).

**Reference.** The affine walk the kernel ships: `λ = e/d` (1M), `λ(x + x')`
(1M), `λ²` (1S), Montgomery's trick at `3(B − 1)/B` products and `1/B` of an
inversion per update. At `B = 16` this is 4.81 products (the 77 per 16 slots
the kernel issues), 1,120 ALU and **38.1 clmad** per update, which is
ITERATION-FUNCTION's "arithmetic floor" row (≈ 1,090 – 1,180 ALU, 38.5 clmad)
re-priced routine by routine. It is the same for the shipping walk `R + σʲ(R)`
and the table walk `R + ε σᵏ(T_h)`: `σ` is a free permutation in the normal
basis and the addend is affine either way.

**Floor for the direction.** The direction replaces "`3M + I/B`" (the
batched inversion) with an inversion-free addition. It can only win if

```
cost(inversion-free addition) + cost(invariant selector) < cost(affine addition with batched inversion)
```

Both sides are formula counts: the left is 8M + 2S for the best known
inversion-free binary-curve addition (λ-projective mixed; López–Dahab mixed
is 8M + 5S or 9M + 4S), the right is `5M + 1S − 3M/B + I/B`. This is stated
before anything is priced, and the table below is the comparison.

**Falsification target.** The direction is alive iff a scale-invariant
selector `h(X, L, Z)` and predicate exist whose cost `c_h` satisfies
`60.5 + c_h < 38.1` clmad per update at `B = 16`, i.e. `c_h < 0`. It is
revived only in a regime where the batch is forced to `B ≤ 2`, where the
affine walk is ALU-bound by the inversion and its ceiling (9.8 – 10.3 B/s at
`B = 2`) drops under the projective one (11.8 – 12.5); there the bar on the
invariant is `c_h ≤ 13` clmad and `≤ 1,300` ALU per update, about two
products, and §3.2 finds nothing under one inversion (48 clmad, 4,192 ALU).
BATCH-TUNING.md measured `B = 8` at 12.6 against `B = 16` at 13.2 B/s and
`B = 32` at 8.7, so the kernel is not in that regime and moving toward it
costs more than it saves.
Inadmissible: counting the squarings as free (they are 6.25 clmad each in the
polynomial basis the products need; moving them to the ALU costs ~100 ALU
each and the ALU is the tighter pipe), pricing the mixed formula for the
shipping walk (whose addend `σʲ(R)` is projective too: 11M + 2S), or a hash
that is invariant on a fraction `1 − p` of representatives with `p` above
`2^−25` (ITERATION-FUNCTION §4.1: trails that meet as points must stay
together for ~2^25 steps).

## 2. The single table

Formula counts verified against affine addition on 300 random pairs at
`m = 23` and `41` and 30 at `m = 131`, with every projective input given a
random `Z` (§4). Prices from §1; the ceiling is updates per second at 100% of
the binding pipe over the clock range, and every row is **priced from a
verified formula, not built**. The one measured row is the reference.

| variant | M | S | I | ALU / update | clmad / update | ceiling B/s at 100% | binding pipe | measured | class |
|---|---:|---:|---:|---:|---:|---:|---|---|---|
| affine, Montgomery batch 1 | 2.00 | 1 | 1 | 4,597 | 66.2 | 5.8 – 6.2 | ALU | | |
| affine, Montgomery batch 2 | 3.50 | 1 | 1/2 | 2,742 | 51.2 | 9.8 – 10.3 | ALU | | |
| affine, Montgomery batch 3 | 4.00 | 1 | 1/3 | 2,124 | 46.2 | 12.6 – 13.4 | ALU | | |
| affine, Montgomery batch 4 | 4.25 | 1 | 1/4 | 1,815 | 43.8 | 14.8 – 15.6 | ALU | | |
| affine, Montgomery batch 8 | 4.62 | 1 | 1/8 | 1,352 | 40.0 | 17.8 – 18.8 | clmad | 12.6 B/s screen (BATCH-TUNING) | |
| **affine, Montgomery batch 16 (shipping)** | **4.81** | **1** | **1/16** | **1,120** | **38.1** | **18.7 – 19.8** | clmad | **14.41 B/s** (ITERATION-FUNCTION §6.3; 2,324 ALU with the walk's non-arithmetic work) | **reference** |
| affine, Montgomery batch 32 | 4.91 | 1 | 1/32 | 1,004 | 37.2 | 19.2 – 20.3 | clmad | 8.7 B/s (BATCH-TUNING: state traffic) | |
| affine, batch → ∞ | 5.00 | 1 | 0 | 888 | 36.3 | 19.7 – 20.8 | clmad | | limit |
| λ-projective mixed (table-walk addend), **invariant hash free** | 8 | 2 | 0 | 1,454 | **60.5** | **11.8 – 12.5** | clmad | not built | below reference |
| λ-projective full (shipping-walk addend `σʲ(R)`), invariant hash free | 11 | 2 | 0 | 1,937 | 78.5 | 9.1 – 9.6 | clmad | not built | below reference |
| λ-projective mixed + invariant hash by one inversion of `Z` per update | 8 | 2 | 1 | 5,646 | 108.5 | 4.8 – 5.0 | ALU | not built | below reference |
| λ-projective mixed + invariant hash by batched inversion of `Z` (16) | 10.81 | 2 | 1/16 | 2,169 | 80.4 | 8.9 – 9.4 | clmad | not built | below reference; and once `1/Z` is in hand the affine step is 3M + 1S cheaper |

Reading the table. The ratio that matters is the λ-projective row against
the reference on the binding pipe: **60.5 / 38.1 = 1.59×** the carry-less
work, with the hash free, and 1,454 / 1,120 = 1.30× the arithmetic ALU. The
ceiling of the best projective row, 12.5 B/s at 100% of the carry-less unit,
is under the 14.41 B/s the affine kernel measures. On the carry-less column
affine is ahead at every batch from 2 (51.2 against 60.5); on the ceiling,
which takes the binding pipe, the crossover is between `B = 2` (affine
ALU-bound at 9.8 – 10.3) and `B = 3` (12.6 – 13.4 against 11.8 – 12.5):
Montgomery's trick pays for itself against the best inversion-free addition
as soon as three walks share an inversion. The two hash rows show that the natural way to get the invariant,
computing `x = X/Z`, costs exactly what affine coordinates cost, and once
`1/Z` is available the affine formulas need 3M + 1S fewer than the projective
ones; "projective plus an invariant obtained by inversion" is affine with
extra steps.

Class of every projective row: **priced, below the reference**; nothing was
built because the priced ceiling is below the measured baseline, which is the
condition ITERATION-FUNCTION §1 uses for its own lever rows. Per `AGENTS.md`
§3 the round is **accounting**: the numbers on the lever row moved from a
sentence to a table, the algorithm did not.

## 3. The invariant hash on its own terms

The rate argument of §2 gives the projective walk the hash for free and it
still loses. This section records why the hash is not free either, so that
the direction is not reopened from that side.

### 3.1 What the walk needs from the hash

Write `f` for the step on representatives: `f(X, L, Z) = (X, L, Z) + Q(h)`,
where the addend `Q(h)` is selected by `h = h(X, L, Z)` (the branch `j` or the
table index and sign). The addition formulas are homogeneous: `f(cX, cL, cZ)`
is a representative of the same point as `f(X, L, Z)` *iff* `h(cX, cL, cZ) =
h(X, L, Z)`. So the walk descends to a function on points exactly when `h` is
scale-invariant, and the same for the distinguished-point predicate, which
must fire on the point and not on the representative. (Both also have to be
invariant under `σ` and negation, as today. `σ` acts coordinate-wise on any
representative and `−(X, L, Z) = (X, L + Z, Z)`, so the same normal-basis
tricks apply once the affine coordinate is in hand; they do not remove the
scaling problem.)

The invariance has to be exact. ITERATION-FUNCTION §4.1 already made the
argument for the negation map: two trails that meet as points must stay
together for the ~2^25 steps to the next distinguished point, so a selector
that disagrees on a fraction `p` of representative pairs splits them after
`1/p` steps and needs `p ≪ 2^−25`. Every representative-based quantity
(weight of `X`, trace of `X`, any bits of `X`, `L` or `Z`) has `p` near
`1 − 1/H` for an `H`-way selector.

### 3.2 What a scale-invariant costs

A scale-invariant function of `(X, L, Z)` is a function of the point
`(X/Z, L/Z)`. The candidates, priced in products (`M`) on `GF(2^131)` with
free squarings in the normal basis:

| candidate | invariant? | cost | verdict |
|---|---|---:|---|
| `HW(X/Z)`, `Tr(X/Z)`, any function of `x` | yes | one inversion: Itoh–Tsujii 8M + 130 free squarings (`ref.h`), or `3M + I/B` batched | the affine walk's cost, with the affine formulas then 3M + 1S cheaper (§2) |
| `HW(X)`, `Tr(X)`, bits of `X`, `L`, `Z` | **no** | free | splits trails at step 1 with probability ~0.8 (§3.3) |
| `X = 0`, `Z = 0`, `X = Z` and other 2-minor tests `X Z' = X' Z` | yes | 0 – 2M | one bit, constant on the walk's points (no point has `x = 0` or `Z = 0`); a 2-minor compares two points, it does not hash one |
| `(X/Z)^k` for small `k` via `X^k · Z^{−k}` | yes | ≥ 1 inversion | contains the inversion |
| the character `(X/Z)^{(2^131 − 1)/263} ∈ μ_263` as `X^e / Z^e` | yes | two exponentiations to `e = (2^131 − 1)/263` (≥ 8M each, 123-bit exponent) plus an inverse in `μ_263` (2M) | ≥ 18M, more than two inversions, for a 263-way hash |
| membership `X/Z ∈ V` for an `F_2`-subspace `V` (`X ∈ Z · V`) | yes | a basis of `Z · V`: `dim V` products, or `X/Z` itself | ≥ 1 inversion for any `V` of dimension ≥ 9 |
| `Tr(X · Z^{2^k})`, `X · Z^{2^k}` and other bilinear forms | **no** (degree `1 + 2^k ≠ 0`) | 1M | homogeneous of nonzero degree; scales by `c^{1+2^k}` |
| any polynomial in `X, L, Z` whose monomials have total degree `d ≢ 0 mod 2^131 − 1` | **no** | | scales by `c^d`; a scale-invariant polynomial function needs every monomial at degree `≡ 0 mod (q − 1)`, and the first non-constant ones, `X^k Z^{q−1−k} = (X/Z)^k`, are exponentiations of the inversion's size (`1/Z = Z^{q−2}`) |

The last row is the general statement. Every function on `F_q³` is a
polynomial, so "no polynomial invariant exists" is false; the true statement
is about degree. Group a polynomial `g(X, L, Z)` by total degree,
`g = Σ_d g_d`; `g(cX, cL, cZ) = Σ_d c^d g_d`, and this equals `g` for all
`c ∈ F_q*` iff `g_d = 0` as a function whenever `d ≢ 0 mod (q − 1)`. The
surviving homogeneous parts have degree `≥ q − 1 = 2^131 − 1`, and the
simplest, `X^k Z^{q−1−k}`, is `(X/Z)^k` on `Z ≠ 0`: evaluating it is an
exponentiation with a 131-bit exponent, the size of the inversion itself
(`1/Z = Z^{q−2}`). So any non-constant invariant on the walk's points is a
function of `x` and `λ`, and computing one from the representative means, in
effect, computing `X/Z` or an exponentiation of the same size. This is a **reference boundary** (the cheapest known, not a
theorem that no cheaper evaluation exists for some specific `h`): the table
above is the list of the ways one would try to evade it, each priced, none
below one inversion.

### 3.3 Measured: a representative-based selector does split merged trails

`lambda_projective.py` walks the shipping iteration `R ← R + σʲ(R)`,
`j = 3 + ((HW/2) mod 8)`, on λ-projective representatives, from two
representatives `(x, λ, 1)` and `(cx, cλ, c)` of the same random point with a
random `c`, on the test curves of `codegen/curves.py`:

| `m` | selector | trails on the same point after 64 steps | first divergence | predicate `HW(·) ≤ w` on the representative |
|---:|---|---|---|---|
| 23 | `HW(X/Z)` (invariant, by inversion) | 300 / 300 | never | |
| 23 | `HW(X)` (representative) | 0 / 300 | step 1: 82.3%, by step 2: 95.7%, mean 1.23 | `w = 8`, true rate 0.103: agrees 0.810, fires spuriously 0.107, misses 0.083 |
| 41 | `HW(X/Z)` | 300 / 300 | never | |
| 41 | `HW(X)` | 0 / 300 | step 1: 86.3%, by step 2: 99.0%, mean 1.15 | `w = 16`, true rate 0.110: agrees 0.787, fires spuriously 0.110, misses 0.103 |

With the invariant selector the twist identity holds exactly: the two trails
are representatives of the same point at every step (the `Z`s differ, the
points do not). With the representative selector they separate at the first
step in four cases of five, as the `1 − 1/8` estimate says, and the
distinguished-point predicate on the representative fires on non-distinguished
points at its own rate `θ` and misses distinguished ones at `θ`: a walk driven
that way is not a function on the class set, has no collision structure, and
reports points the host re-walk would reject. This is the measurement behind
the "invariance must be exact" sentence of §3.1; on `GF(2^131)` with
`θ = 2^−25` nothing changes except that the spurious reports would be rarer
and the trails would still part at step 1.

## 4. Formula verification

The λ-projective formulas priced in §2, as implemented and checked
(`checkFormulas`), with every projective input scaled by a random nonzero
`Z` and the result compared to the generator's affine `Curve.add` after
converting back:

```
mixed  (X_P, L_P, Z_P) + (x_Q, λ_Q):                              8M + 2S
  A = L_P + λ_Q Z_P            B = (X_P + x_Q Z_P)²
  X = (A X_P)(A x_Q Z_P)       Z = A B Z_P
  L = (A x_Q Z_P + B)² + A B (L_P + Z_P)

full   (X_P, L_P, Z_P) + (X_Q, L_Q, Z_Q):                        11M + 2S
  A = L_P Z_Q + L_Q Z_P        B = (X_P Z_Q + X_Q Z_P)²
  X = (A X_P Z_Q)(A X_Q Z_P)   Z = (A B Z_Q) Z_P
  L = (A X_Q Z_P + B)² + (A B Z_Q)(L_P + Z_P)
```

| `m` | random additions checked | mixed | full | `R + σʲ(R)` with a projective `σʲ(R)` | negation `(X, L + Z, Z)` |
|---:|---:|---|---|---|---|
| 23 | 300 | 8M + 2S, all agree | 11M + 2S, all agree | all agree | all agree |
| 41 | 300 | 8M + 2S, all agree | 11M + 2S, all agree | all agree | all agree |
| 131 | 30 | 8M + 2S, all agree | 11M + 2S, all agree | all agree | all agree |

The same runs check the λ-affine addition reduced to one inversion
(`6M + 2S + 1I`, §6) and the López–Dahab mixed addition (`8M + 5S`), each
against the same affine `Curve.add`, with the same result at all three sizes.

The counts are taken by a counting wrapper on the field on every call and
asserted constant across inputs; a first draft of the full formula came out at
12M because `A B Z_Q` was formed twice, and the shared product is what brings
it to the published 11M. The mixed formula assumes `x_P ≠ x_Q`, which the walk
guarantees except on the `2^−131` event that the affine walk's division by
`x + x'` already has.

## 5. What this leaves open

Nothing in the direction as stated. The two things that could reopen a
projective row are outside it:

- **A regime with `B ≤ 2`.** If state traffic ever forced the batch below 3
  (it favours 16 over 32 today for that reason, and 16 over 8), the
  projective row with a free hash would win on the ceiling and the
  invariant-hash question would become the binding one; §3.2 says the hash
  would then cost an inversion, and projective-plus-inversion (108.5 clmad,
  5,646 ALU) is under `B = 1` affine (66.2, 4,597) on both pipes. The row
  loses at every batch once the hash is paid.
- **A cheaper invariant than `X/Z`.** §3.2 is a reference boundary. At
  `B = 16` no invariant revives the row, since the bar is `38.1 − 60.5 < 0`
  clmad: even a free one loses. Only in a hypothetical `B = 2` regime would
  the bar be positive (about two products, §1), and §3.2 lists nothing under
  one inversion. So the hash question decides nothing at any batch the
  kernel has run, and would decide against the direction at the one batch
  where it matters.

The lever row in ITERATION-FUNCTION.md §1 stands, with this note as its
price. No measured number moved; the campaign status page carries the table
of §6 as a static section that cites this note, and the site build pins its
figures to it.

## 6. Every point representation, priced

The question generalises: is there *any* point representation on this curve
whose addition, plus the class-invariant selector it needs, costs less per
step than Weierstrass affine with Montgomery's trick? The rows below extend
the table of §2 to the other representations. Verified rows were checked by
`checkFormulas` against the generator's affine addition at `m = 23, 41, 131`
(300, 300 and 30 random pairs, random `Z` on every projective input) with
their operations counted; the last row is from the literature and is marked
so. Prices as in §1; the ceiling is at 100% of the binding pipe.

| coordinates | addition | M | S | I | clmad / update | ALU / update | ceiling B/s | invariant selector | status |
|---|---|---:|---:|---:|---:|---:|---:|---|---|
| **Weierstrass affine `(x, y)`, Montgomery batch 16** | `1I + 2M + 1S`, batched | 4.81 | 1 | 1/16 | **38.1** | 1,120 | 18.7 – 19.8 | `HW(x)`, free in the normal basis | **shipping, measured 14.41 B/s** |
| Weierstrass affine, batch 8 | same | 4.62 | 1 | 1/8 | 40.0 | 1,352 | 17.8 – 18.8 | same | measured 12.6 B/s screen (BATCH-TUNING) |
| Weierstrass affine, batch 32 | same | 4.91 | 1 | 1/32 | 37.2 | 1,004 | 19.2 – 20.3 | same | measured 8.7 B/s (state traffic) |
| λ-affine `(x, λ)`, batch 16 | `1I + 6M + 2S` reduced to one inversion, batched | 8.81 | 2 | 1/16 | 68.4 | 1,847 | 10.4 – 11.0 | same `x` | verified formula, below reference |
| λ-projective `(X, L, Z)`, mixed addend | `8M + 2S` | 8 | 2 | 0 | 60.5 | 1,454 | 11.8 – 12.5 | none without `X/Z`; given free here | verified, below reference |
| λ-projective, projective addend `σʲ(R)` | `11M + 2S` | 11 | 2 | 0 | 78.5 | 1,937 | 9.1 – 9.6 | given free | verified, below reference |
| López–Dahab `(X, Y, Z)`, `x = X/Z`, `y = Y/Z²`, mixed addend | `8M + 5S` (the `a Z²` term is a constant multiply on a Koblitz curve) | 8 | 5 | 0 | 79.2 | 1,703 | 9.0 – 9.5 | given free | verified, below reference |
| λ-projective mixed, invariant by one inversion of `Z` | `8M + 2S + 1I` | 8 | 2 | 1 | 108.5 | 5,646 | 4.8 – 5.0 | `HW(X/Z)` | below reference on both pipes |
| binary Edwards, Hessian, Huff, Jacobian | published inversion-free additions from about `10M` to `16M + 1S + 4D`; every affine form has more than one denominator per addition, so it costs more products than Weierstrass to fold into one inversion | | | | | | | negation is not a coordinate map on Edwards (`(x, y) ↦ (y, x)`), so the invariant needs symmetric functions | published counts, not re-derived here, not priced |

Reading it. Two structural facts make the column monotone. Every projective
model exists to remove the inversion, and Montgomery's trick already
amortises the inversion to three products plus a sixteenth of the chain, so
an inversion-free formula has to beat about five products per step and none
is below eight. Every affine model other than Weierstrass has more than one
denominator per addition (λ-affine has `B` and `AB`), so folding them into
one inversion costs the products the model was meant to save. The
Itoh–Tsujii chain is at the addition-chain lower bound for exponent 130
(`⌊log₂ 130⌋ + HW(130) − 1 = 8` products), so the inversion cannot get
cheaper either. What is left on the representation axis is the field and
state representation, which have their own notes (polynomial against normal
basis for the products, packed and compact state for the traffic that
decided batch 16 over 32, the denominator cache).

Two ideas that look like representation changes and are not: on this curve
`τ² + τ + 2 = 0`, so `R + σʲR = [1 + τʲ]R` and `x([1 + τʲ]R)` is a rational
function of `x(R)` of degree `1 + 2ʲ + Vⱼ`, which is 14 already at `j = 3`,
more than an addition; and x-only (Kummer-line) arithmetic needs
`x(P − Q)`, which a random-walk addend does not supply, while a fixed addend
makes the walk a permutation rather than a random function.
