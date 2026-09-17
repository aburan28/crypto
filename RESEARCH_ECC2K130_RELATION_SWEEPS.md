# Homogeneous relation sweeps on ECC2K-130, and what rank is worth

**Question.** Can index calculus built on homogeneous point relations beat
optimised parallel Pollard rho on ECC2K-130?

**Answer.** No, and the reason is a counting argument rather than a failure
to optimise. The sweeps in this note are a confirmation of that argument and
a calibration of the instrument; they are not where the answer comes from.

The instance is the Certicom challenge curve `y^2 + xy = x^3 + 1` over
`F_2^131`, cofactor 4, prime subgroup order

    r = 680564733841876926932320129493409985129,   log2 r = 129.0

## 1. The boundary, stated before measuring

Per `AGENTS.md` §1 this thread has both kinds of boundary.

**The reference** is Pollard rho with the `<-1> x <pi>` speed-up on the same
subgroup: `2^60.9` iterations.

**The floor** is a counting argument. A factor base of `B` abscissae admits
about `B^m / m!` unordered `m`-multisets, each summing to a fixed point with
probability about `1/r`, so `m`-point relations begin to exist only at

    B_m = (m! r)^(1/m).

Below that threshold the expected yield is under one and a *complete* sweep
finds nothing, however the base is chosen. This part is independent of
structure: it counts multisets, and a structured base of size `B` has exactly
as many as a random base of size `B`.

Finding a relation once they exist is a second cost, best done by
meet-in-the-middle on a balanced split. Recovering a logarithm is a third,
and it is the one that is easy to drop: homogeneous relations among base
points give linear dependencies between *unknown* logarithms, not `log_P(Q)`.
For a logarithm one needs about `B` relations whose targets are known
combinations `[a]P + [b]Q`, so the relation cost is paid `B` times over.

| `m` | `log2 B_m` | split | `log2` memory (bytes) | `log2` cost, one relation | vs rho | `log2` total | vs rho |
|---:|---:|:--:|---:|---:|---:|---:|---:|
| 3 | 43.86 | 1+2 | 48.86 | 86.72 | +25.82 | 130.58 | +69.68 |
| 4 | 33.40 | 2+2 | 70.79 | 65.79 | **+4.89** | 99.19 | +38.29 |
| 5 | 27.18 | 2+3 | 58.36 | 78.96 | +18.06 | 106.14 | +45.24 |
| 6 | 23.08 | 3+3 | 71.66 | 66.66 | +5.76 | 89.74 | +28.84 |
| 7 | 20.19 | 3+4 | 62.97 | 76.16 | +15.26 | 96.34 | +35.44 |
| 8 | 18.04 | 4+4 | 72.56 | 67.56 | +6.66 | **85.60** | **+24.70** |

No `m` reaches the reference. The friendliest accounting imaginable — price
a single relation and ignore the `B` relations a logarithm actually needs —
still puts the best case, `m = 4`, at `2^65.79` against rho's `2^60.9`, a
factor `2^4.89`, while asking for `2^70.79` bytes of storage. That is about
two zettabytes, and it is not a constant that engineering moves. Priced
honestly at `B` relations the smallest gap is `2^24.7`, at `m = 8`.

Quoting the single-relation row as the method's cost is exactly the §3
**relabelling** error: the work has not gone away, it has moved to a column
the headline does not look at.

**Falsification target, stated in advance.** This thread would be a success
if some support of size `B` produced homogeneous relations at a rate above
`B^m / (m! r)` by a factor large enough to bring the total below `2^60.9`,
with every relation verified on the curve and with construction equations
excluded from the count. It should be abandoned if measured yields match the
counting prediction on every support tried, since then no structure is
evading the floor.

## 2. The one way structure could have won

The floor rests on an independence assumption: that each `m`-multiset sums to
the identity with probability `1/r`, as a random set would. That assumption
is not a theorem, and it is false for at least one family — a base built
deliberately from `P` and `Q` carries relations that hold by construction,
such as `R_{1,0} + R_{0,1} = R_{1,1}`.

So the experiments below are aimed at the assumption, not at the count. The
question each one asks is whether a structured support produces relations, or
pair-sum collisions, above the random-base prediction.

## 3. What a relation is, and what it is worth

Three points sum to `O` exactly when they are collinear. Substituting
`y = lx + m` into the curve equation gives

    x^3 + (a + l^2 + l) x^2 + m x + (b + m^2) = 0

so a collinear triple has `e2 = m` and `e3 = b + m^2`, i.e.

    (x1 + x2)^2 x3^2 + (x1 x2) x3 + (x1 x2)^2 + b = 0

which is Semaev's `S3`, already present in the repository as `s3_in_last`.
A pair therefore *fixes a quadratic*, and the third abscissa is two
candidates to look up rather than a base to scan. That is what makes a
complete sweep of a multi-million-pair base affordable. The curve's `a` does
not enter this condition at all; it enters only the solvability of
`l^2 + l = a + e1`, which is checked separately.

Writing each base point as `R_i = [u_i] P + [v_i] Q`, a relation gives

    sum(u_i) + log_P(Q) * sum(v_i) = 0   (mod r)

which determines `log_P(Q)` only when `sum(v_i)` is invertible mod `r`. A row
with `sum(u) = sum(v) = 0` is a **construction equation**: it holds whatever
the logarithm is. Such rows are real relations, they contribute real rank to
a relation matrix, and they carry exactly zero information.

**This is the accounting result of the thread.** Relation rank is not
evidence until the construction span has been quotiented out. A base
manufactured from `P` and `Q` can be made to produce as many equations as one
likes, of as high a rank as one likes, none of which say anything about
`log_P(Q)`.

## 4. The four-point experiment, and a superseded premise

The merged form of this study (#404) recorded a four-point sweep whose design
rested on a claim that is true but was applied too widely:

> The weight-two set is **not** Frobenius-stable in the polynomial basis.
> [...] This rules out normalising a relation by Frobenius so that one
> summand sits at an orbit representative *of the base*, because the base is
> not a union of orbits.

The first sentence is correct. The second does not follow, and the word doing
the damage is *the*. Weight-two is a property **of a basis**, not of the
field. In the challenge's polynomial basis `z^i + z^j` squares out of weight
two as soon as `2i >= 131`, and the merged note is right that only the
`i < j <= 65` part survives. In a **normal** basis it cannot happen at all,
because squaring is a cyclic shift of coordinates. Writing `alpha` for a
normal element and taking the support

    x_{i,j} = alpha^(2^i) + alpha^(2^j),     0 <= i < j < 131

Frobenius sends `x_{i,j}` to `x_{i+1, j+1}` with indices read mod 131, so the
support is closed under `sigma` by construction. The orbit-representative
normalisation the merged note ruled out is available after all; it was ruled
out for the polynomial basis and the conclusion was carried to a support that
does not have to be built in one.

Measured here over eight random normal elements of `F_2^131`: the weight-two
normal support is `sigma`-stable in all eight, with no exceptions and no
boundary cases.

### 4.1 The support is a union of orbits, and the orbits are all of size 131

`sigma` acts on index pairs by `(i, j) -> (i+1, j+1) mod 131`, so an orbit is
determined by the difference `d = j - i` read up to sign, and `131` is prime,
so every orbit has exactly `131` members. There are

    C(131, 2) / 131 = 8515 / 131 = 65

orbits, one per `d` in `1..65`. **Not 66** — an earlier count of this said 66
and was wrong by one; there is no `d = 0` class, since `i < j`.

An abscissa lifts to a curve point when the relevant trace condition holds,
and `sigma` preserves the curve, so lifting is constant on an orbit: an orbit
lifts entirely or not at all. Measured over the same eight bases, with zero
mixed orbits in any of them. A normal-basis weight-two support is therefore
always `131 * (number of lifting orbits)` points.

### 4.2 Support size is basis-dependent, over a wider range than previously recorded

An earlier note put the lifting-orbit count at `28` of `66`, giving the
`B = 3668` support this study has been quoting, and said the count "runs
28–36 across normal elements". Both the denominator and the range are wrong.
Eight random normal elements, measured, with the run of §4.6 alongside:

| `alpha` (leading bits) | lifting orbits (of 65) | `B` | pair-orbits swept | candidates | certified | `log2` expected usable |
|:--|---:|---:|---:|---:|---:|---:|
| `0x4f7aaf673fb46bcf…` | 23 | 3013 | 69,276 | 46 | 0 | -84.36 |
| `0x2de8f4e94249b03f…` | 25 | 3275 | 81,850 | 50 | 0 | -83.878 |
| `0x57eea059506531b2…` | 26 | 3406 | 88,530 | 52 | 0 | -83.652 |
| `0x4dea6516dc1ed48b…` | 29 | 3799 | 110,142 | 58 | 0 | -83.022 |
| `0x1a652f4bb357676c…` | 31 | 4061 | 125,860 | 62 | 0 | -82.637 |
| `0x45afb07e8cb20ea5…` | 32 | 4192 | 134,112 | 64 | 0 | -82.453 |
| `0x33471bc1d825ddd3…` | 37 | 4847 | 179,302 | 74 | 0 | -81.615 |
| `0x4d44a5aac26a00c4…` | 42 | 5502 | 231,042 | 84 | 0 | -80.884 |

The spread is `23–42` orbits, `B = 3013–5502`, not `28–36`. `B = 3668` is one
draw from this distribution and carries no special status. Nothing in the
conclusion turns on which draw is used — every one is astronomically below
the `B_4` of §1 — but a quoted support size should not be mistaken for a
property of the curve.

### 4.3 Frobenius acts as a scalar, so a support carries one unknown per orbit

`E` is a Koblitz curve, `#E(F_2) = 4`, so the base-field trace is `t = -1` and
`sigma` satisfies `sigma^2 + sigma + 2 = 0`. On the order-`r` subgroup this
has the root

    s = 196511074115861092422032515080945363956

verified here directly: `sigma(P) = [s]P` and `sigma(Q) = [s]Q` for the
challenge points, and `s^131 = 1 mod r`.

So for a support point `R` and its orbit, `log(sigma^k R) = s^k log(R)`. The
131 points of an orbit share **one** unknown. A 28-orbit support presents
3668 points to a sweep and 28 unknowns to the linear algebra.

**This is what gives the negative its force.** The usual reason an index
calculus attempt fails inconclusively is that the linear algebra was never
reached — not enough relations to fill a `B`-column matrix. Here `B` is
effectively the orbit count. About 29 relations would close a 28-orbit
system, and the solve is free at that size. The method still produces none,
so the failure cannot be attributed to an unaffordable second stage. It is a
supply failure, and supply is what §1 counts.

The same collapse is exercised where it can be checked: `planted.py` plants
`d`, recovers it from four-point relations against known targets
`[a]P + [b]Q`, and verifies `[d]P = Q`. Recovered at `m = 13` (`B = 52, 78`)
and `m = 19` (`B = 114, 152`), with one unknown per orbit plus `d`. No
discrete logarithm enters the solve: the Frobenius scalar comes from the
characteristic equation, with the root chosen by a point identity.

### 4.4 Where the sums actually live, and what `E[4]` is worth

The cofactor is 4, so a four-point sum landing anywhere in `E[4]` still gives
a usable equation once the cofactor is cleared. Demanding `sum = O` throws
usable relations away. The first correction of this priced the gain at `4x`,
by assuming the four cosets are hit uniformly. They are not.

Every point of a normal-basis weight-two support lies in the **index-2
subgroup** `H` of order `2r` — measured on the challenge curve and on the
small analogues at `m = 11, 13, 19`, with no exceptions. Four-point sums
therefore lie in `H` always, and

    E[4] ∩ H = {O, the order-2 point},   of size 2.

The two order-4 points are unreachable, and exhaustive counts put **exactly
zero** relations on them across every support tried. So the coset factor is
exactly **2** — not the 4 of the uniform guess, and not the `2.33` recorded
earlier, which was a measurement taken before the counter producing it had
been checked against brute force.

Note which way this cuts. Exact relations run at *twice* the rate a uniform
model over `E` predicts, which looks like structure finally beating the
independence assumption of §2. It is not: the usable rate, the one that
matters, lands exactly on `8 C(B,4) / r`. The bias moves relations between
cosets without making any more of them.

### 4.5 The supply model, and what it is a model of

The counting floor of §1 is about *multisets*. A search does not enumerate
multisets. It enumerates four distinct support abscissae with a sign on each,
modulo a global negation, so the search space is

    configurations = C(B,4) · 2^4 / 2 = 8 · C(B,4)

and, with §4.4,

    expected usable = 8 · C(B,4) / r,     expected exact = 8 · C(B,4) / 2r.

The abscissae must be **distinct**. Allowing a repeat admits
`R + (-R) + T + (-T) = O`, which is a relation, holds for every `R` and `T`,
and says nothing about any logarithm — the construction equation of §3 in its
purest form. Counting those is what made an earlier version of this model
miss measured yields by `2x`–`43x`, and it did so *in the validation code*,
where an inflated model makes a correct search look broken.

Checked against exhaustive ground truth at `m = 11, 13, 19`, over every
support with a predicted yield above 100: mean measured/predicted **1.03**
for usable relations and **1.08** for exact, with the coset factor measured at
**1.92**. At `B = 3668` the model gives `2^-83.2` expected usable relations.
An earlier figure of `2^-78` was wrong by `2^5`.

### 4.6 The run

Eight independent normal bases, the full canonical pair-orbit enumeration on
each, no budget exhausted — the largest basis used 71 s of a 3600 s budget,
and every basis completed.

    pair-orbits swept        1,020,114   (= 133,634,934 signed pairs, 2^27.0)
    candidates               490
    certified relations      0
    summed expectation       2^-79.4 usable relations

Zero is the predicted answer, to within a factor the experiment has no power
to resolve. The run is not evidence because it found nothing; it is evidence
because the same pipeline recovers planted logarithms at `m = 13` and
`m = 19` (§4.3) and reproduces exhaustive relation counts at three small
degrees (§4.5).

### 4.7 What the collisions are

The 490 candidates are not noise, and calling them degenerate undersells
them. Every one, on every basis, is an instance of

    2A + sigma A + sigma^2 A = O        or        2A - sigma A - sigma^3 A = O

and the count is exactly `2 x (lifting orbits)` on every basis — 50 at 25
orbits, 84 at 42, and so on down the table in §4.2. Both patterns are
consequences of `sigma^2 + sigma + 2 = 0`: the first is the characteristic
equation itself, the second is it again after `sigma^3 = 2 - sigma`. They
hold on **every point of the curve**, which is how the run classifies them —
by re-evaluating the pattern at random curve points having nothing to do with
the support, rather than by deriving it.

So the orbit-grouped search, run to completion on a million pair-orbits,
recovers the endomorphism ring's own relation and nothing else. These are
construction equations in the exact sense of §3 — real relations, real rank,
zero information — and they are the *only* thing a Frobenius-quotiented
four-point search on this support turns up.

There is a lesson in the shape of that. Quotienting by an endomorphism makes
the search cheaper and makes the endomorphism's own identities appear as
hits. A pipeline that counted collisions rather than certifying them would
report a yield here that rises linearly with support size and means nothing.

### 4.8 What the Frobenius quotient is worth against `m = 4`

Section 1 prices `m = 4` at `2^65.79` for one relation and `2^99.19` for a
logarithm. That accounting predates the quotient, which changes three of its
inputs: the floor moves, because the search space is `8 C(B,4)` and the
cofactor is worth 2; storage and streaming divide by 131, with an 8-byte
class key in place of a 32-byte point; and the relations needed divide by 131
as well, because the unknowns are orbits. Keeping §1's convention throughout
so the rows can be read against each other
(`results/four_point_orbit_boundary.json`):

| accounting | `log2 B` | `log2` one relation | vs rho | memory (EB) | `log2` total | vs rho |
|:--|---:|---:|---:|---:|---:|---:|
| section 1, as published | 33.40 | 65.79 | **+4.98** | 2,044.8 | 99.19 | **+38.38** |
| corrected floor | 32.65 | 64.29 | **+3.48** | 723.0 | 96.94 | **+36.13** |
| + orbit-grouped pair sums | 32.65 | 57.26 | **-3.55** | 1.4 | 89.91 | **+29.10** |
| + one unknown per orbit | 32.65 | 57.26 | **-3.55** | 1.4 | 82.87 | **+22.06** |

Two things to read off this.

The single-relation row goes **below** the reference — `2^57.26` against
rho's `2^60.81`. That row is precisely the one §1 warns about. The work has
not gone away; it has moved into the count of relations a logarithm needs,
and the honest row is the last column.

Priced there, the quotient is worth `2^14` — a factor of `131^2`, `131` from
the storage and `131` from the unknowns — and it makes orbit-grouped `m = 4`
the best accounting anywhere in this study, at `2^+22.06` over rho, ahead of
the `2^+24.70` that `m = 8` holds in §1. It is still `2^22` short. The
Frobenius quotient is a real and quantified saving, and it is about a quarter
of the way, in the exponent, to a result.

### 4.9 Six defects, each of which produced plausible output

Recorded because each returned numbers that looked like results, and three of
them were in validation code rather than in the search — the code whose job
is to catch the other kind.

1. **Double-counted pairs.** Every unordered pair was enumerated twice, once
   anchored at each element's orbit, and the first orbit-grouped run reported
   125,953 "relations". A pair whose elements lie in different orbits can be
   rotated to put either at position 0, and both results are anchored. The
   enumeration is now by `(t1, t2, delta)`, which names each orbit exactly
   once; a test checks the count three ways.
2. **Supply model never validated.** The expected-relation model was
   asserted, not checked. Its first test missed measured yields by
   `2x`–`43x`.
3. **Cancelling multisets counted as relations** — the cause of (2), and the
   §3 error reappearing inside the code written to validate against it.
4. **Detector too strict.** The search required the sum to be exactly `O`
   when the attack only needs `E[4]`.
5. **Coset factor guessed, then mismeasured.** The correction to (4) was
   first priced at `4x` by assuming uniform cosets. A later measurement gave
   `2.33x`. Both are wrong: it is exactly 2, for the structural reason in
   §4.4, and the `2.33` came from the same unvalidated counter as (2).
6. **Orbit points signed independently.** Support points were chosen one per
   abscissa by least `y`. An abscissa carries two points, so this flips the
   sign at about half the positions of each orbit. The collision search
   cannot see it — it reads only abscissae — but it silently breaks
   `log(sigma^k R) = s^k log(R)` and so corrupts the solve, in a system that
   still looks perfectly well-formed. Found by the planted-recovery control,
   which returned a wrong answer rather than an error.

Defects 2, 3 and 5 were in validation code. That is the pattern worth
carrying forward: on this thread the instrument has been more reliable than
the things built to check it.
