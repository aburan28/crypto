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
| 3 | 43.86 | 1+2 | 48.86 | 86.72 | +25.91 | 130.58 | +69.78 |
| 4 | 33.40 | 2+2 | 70.79 | 65.79 | **+4.98** | 99.19 | +38.38 |
| 5 | 27.18 | 2+3 | 58.36 | 78.96 | +18.15 | 106.14 | +45.33 |
| 6 | 23.08 | 3+3 | 71.66 | 66.66 | +5.85 | 89.74 | +28.93 |
| 7 | 20.19 | 3+4 | 62.97 | 76.16 | +15.35 | 96.34 | +35.53 |
| 8 | 18.04 | 4+4 | 72.56 | 67.56 | +6.76 | **85.60** | **+24.79** |

No `m` reaches the reference. The friendliest accounting imaginable — price
a single relation and ignore the `B` relations a logarithm actually needs —
still puts the best case, `m = 4`, at `2^65.79` against rho's `2^60.809`, a
factor `2^4.98`, while asking for `2^70.79` bytes of storage. That is about
two zettabytes, and it is not a constant that engineering moves. Priced
honestly at `B` relations the smallest gap is `2^24.79`, at `m = 8`.

(The ratios above are against the repository's frozen reference `2^60.809`,
as `results/boundary.json` computes them. An earlier printing of this table
quoted them against a rounded `2^60.9` and so ran `0.09` low in every "vs
rho" column.)

Quoting the single-relation row as the method's cost is exactly the §3
**relabelling** error: the work has not gone away, it has moved to a column
the headline does not look at.

**This table has since been superseded, in margin but not in verdict.** It
takes no Frobenius quotient, and it measures cost at a support size whose
expected yield is one relation while charging for `B` of them. §4.8 redoes
it with the quotient applied at every `m` and on a support that actually
yields what it is charged for: `m = 8` is still cheapest, and the best row
moves from `2^+24.79` to `2^+13.59` over the reference. Read the numbers
above as this thread's starting accounting, not its final one.

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

### 4.8 What the Frobenius quotient is worth, at every `m`

Section 1 prices `m = 3..8` with no quotient and finds `m = 8` cheapest at
`2^+24.79` over the reference. The orbit-grouped search changes three of
that accounting's inputs, and they apply at every `m`, not only at `m = 4`:

* the **floor** moves. A search enumerates `m` distinct abscissae with a
  sign each, modulo a global negation — `2^(m-1) C(B,m)` — and by §4.4 the
  sums reach only the two elements of `E[4] ∩ H`. Expected usable relations
  are `2^(m-1) C(B,m) / r`, so `B_m = (m! r / 2^(m-1))^(1/m)` rather than
  §1's `(m! r)^(1/m)`;
* **storage and streaming** divide by 131, both sides of the
  meet-in-the-middle coming in `sigma`-orbits, with an 8-byte class key in
  place of a 32-byte point;
* the **relations needed** divide by 131, the support carrying one unknown
  per orbit.

Keeping §1's convention — cost measured at the floor — so the rows read
against each other (`results/four_point_orbit_boundary.json`):

| `m` | `log2 B` | split | `log2` one relation | vs rho | memory (EB) | `log2` total | vs rho |
|---:|---:|:--:|---:|---:|---:|---:|---:|
| 3 | 43.19 | 1+2 | 78.36 | +17.55 | 0.0 | 114.52 | **+53.71** |
| 4 | 32.65 | 2+2 | 57.26 | -3.55 | 1.4 | 82.87 | **+22.06** |
| 5 | 26.38 | 2+3 | 69.53 | +8.72 | 0.0 | 88.87 | **+28.06** |
| 6 | 22.25 | 3+3 | 57.13 | -3.68 | 1.3 | 72.34 | **+11.53** |
| 7 | 19.33 | 3+4 | 65.70 | +4.89 | 0.0 | 77.99 | **+17.18** |
| 8 | 17.16 | 4+4 | 57.03 | -3.78 | 1.2 | 67.16 | **+6.35** |

**A correction.** An earlier reading of the `m = 4` row called
orbit-grouped `m = 4` the best accounting in this study, "ahead of the
`2^+24.79` that `m = 8` holds in §1". That compares a quotiented row
against an unquotiented one, and it is wrong. Applied uniformly the
quotient helps *larger* `m` more — `2^16.07` at `m = 3` rising to
`2^18.44` at `m = 8` — because the storage saving of `131` is a larger
fraction of a smaller stored side. So `m = 8` remains cheapest, §1's
ranking is unchanged, and what the quotient moves is the margin, not the
order.

The `m = 4`, `m = 6` and `m = 8` single-relation rows all fall **below** the
reference. Those are precisely the rows §1 declines to quote. The work has
not gone away; it has moved into the count of relations a logarithm needs,
and the honest column is the last one.

### 4.8.1 Section 1's convention is generous, and it now matters

§1 measures cost at `B_m`, where the expected yield is **one** relation,
while charging for `B` relations. That is internally inconsistent, in the
method's favour. At the margins §1 was quoting — `2^+24.79` at best — the
slack did not change any conclusion. With the quotient applied it does, so
the honest floor is the support size at which the yield actually reaches
the number of relations needed:

    2^(m-1) C(B,m) / r = B / 131   =>   B = (m! r / (131 · 2^(m-1)))^(1/(m-1))

| `m` | `log2 B` | split | `log2` one relation | vs rho | memory (EB) | `log2` total | vs rho |
|---:|---:|:--:|---:|---:|---:|---:|---:|
| 3 | 61.28 | 1+2 | 114.52 | +53.71 | 0.2 | 168.76 | **+107.95** |
| 4 | 41.18 | 2+2 | 74.33 | +13.53 | 190,500.0 | 108.48 | **+47.68** |
| 5 | 31.22 | 2+3 | 84.04 | +23.23 | 0.2 | 108.22 | **+47.41** |
| 6 | 25.29 | 3+3 | 66.26 | +5.45 | 705.2 | 84.51 | **+23.71** |
| 7 | 21.38 | 3+4 | 73.89 | +13.08 | 0.2 | 88.24 | **+27.43** |
| 8 | 18.61 | 4+4 | 62.82 | +2.01 | 65.1 | 74.40 | **+13.59** |

**The result.** `m = 8` at `2^+13.59` over rho, needing
`65` exabytes of
class keys. No `m` reaches the reference, so §1's answer stands — but the
gap at the best row is `2^11.20` smaller than §1 records, and that is a
material narrowing rather than a rounding.

Two things keep it a negative. The margin is still a factor of about 12,000
in operations, and it is bought with `65` exabytes of storage that no
engineering removes — the same wall §1 identified, moved but not breached.
And the quotient is not a free lever: it requires the support to be a union
of `sigma`-orbits, which is what §4.7 shows makes the endomorphism ring's
own identities appear as hits. A search that took the `2^18.4` saving and
counted collisions without certifying them would be measuring its own
symmetry.

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

## 5. Pushing it: what a logarithm costs, and where the table was flattering

§4.8 stopped at `m = 8` because §1 did. That was an arbitrary stopping point,
and the obvious thing to try is more of what was working. The result is worth
recording in full, including the part that looked like a win.

### 5.1 Extended past `m = 8`, the homogeneous accounting goes below rho

Continuing §4.8.1's self-consistent, quotiented table to larger `m`, with the
support constrained to be a union of `sigma`-orbits (so `B >= 131`):

| `m` | `log2 B` | `log2` total | vs rho |
|---:|---:|---:|---:|
| 4 | 41.18 | 108.48 | **+47.68** |
| 6 | 25.29 | 84.52 | **+23.71** |
| 8 | 18.61 | 74.39 | **+13.59** |
| 10 | 14.97 | 68.86 | **+8.06** |
| 12 | 12.71 | 65.41 | **+4.60** |
| 14 | 11.18 | 63.05 | **+2.25** |
| 16 | 10.08 | 61.36 | **+0.56** |
| 18 | 9.26 | 60.10 | **-0.71** |
| 20 | 8.63 | 59.12 | **-1.69** |
| 22 | 8.14 | 58.34 | **-2.47** |
| 24 | 7.74 | 57.71 | **-3.10** |
| 26 | 7.41 | 57.19 | **-3.62** |
| 28 | 7.14 | 56.76 | **-4.05** |

It crosses at **`m = 18`** and reaches `2^-4.05` at `m = 28`. Taken at face
value that is index calculus beating optimised rho on ECC2K-130.

It is not. It is §3, in the strongest form this thread has produced.

### 5.2 Those are homogeneous relations, and they determine nothing

Every row above prices relations *among base points only*. §3 says what such
a relation is worth:

> Writing each base point as `R_i = [u_i] P + [v_i] Q`, a relation gives
> `sum(u_i) + log_P(Q) * sum(v_i) = 0 (mod r)`, which determines `log_P(Q)`
> only when `sum(v_i)` is invertible mod `r`.

The base points here are not built from `P` and `Q` — they are whatever the
normal basis produces — so their logarithms are unknown and a homogeneous
relation is one linear equation among unknowns. Collecting `T + 1` of them
gives a consistent homogeneous system whose solution space contains the truth
and says nothing about which point of it is true. The cost falls below rho
because the thing being bought has been quietly swapped for something
cheaper. §1 warned about exactly this and the warning still bites at `m = 18`.

A logarithm needs relations against **known targets** `[a]P + [b]Q`.

### 5.3 The cost of a logarithm

Unknowns are the `T` orbit logarithms plus `d`, so `T + 1` relations are
needed. Each requires decomposing a random known target into `n` signed base
points, by meet-in-the-middle.

**The table is built once, not once per target.** An earlier form of this
section charged the whole meet-in-the-middle on every target attempt. That is
wrong, and ordinarily so: the stored side holds `s`-subset sums *of the
support*, which do not depend on the target at all. It is built once and
streamed against for every attempt afterwards.

Amortising it moves the optimum from `2^74.945` to
`2^67.287`. Those are **two separately optimised
configurations**, not one attack repriced — the correction changes which
attack is cheapest, and quoting the difference as the price of a fixed attack
would be wrong. Held at this section's own configuration, rebuilding per
attempt costs far more than the gap between the optima.

**A probe is not free, and there are two of them per streamed point.** The
first form of this model counted probes and not their cost. Each probe is a
group addition, and against a quotiented table it also needs the probe point
canonicalised over its own `sigma`-orbit before it can be looked up — measured
on this container at **13.6 group-operation
equivalents** (29.20 us per batched addition against 398.52 us per
canonicalisation). And each streamed point is tested against both elements of
`E[4] ∩ H`, not one. Together those omissions were worth `2^2.65`.

That makes the table layout a real choice, so both are searched:

* **quotiented** — one canonical class per orbit, 131x fewer entries to
  build, every probe canonicalised;
* **full** — all 131 rotations stored, probes are bare lookups, 131x the
  entries to build.

With memory free it is a pure build-versus-probe trade. Quotiented wins here,
but only by `2^1.4`.

Optimising over support size, relation length, split and table layout
(`results/target_boundary.json`):

    support           2 Frobenius orbits, B = 262
    relations         12 points each, split 11+1
    stored table      quotiented
    build (once)      2^66.777
    probes / attempt  2^10.033 at 2^3.868 each
    target attempts   2^51.64
    relations needed  3
    total             2^67.287      vs rho  +6.48

`2^6.48` is a factor of about 89.
The error bars that remain:

* **memory is charged at zero** — 1,011 exabytes,
  free and instantaneous. Charging it makes the negative larger;
* the build assumes one representative per `sigma`-class is enumerable in
  constant amortised time (necklace enumeration over `Z/131`); naively it
  costs `2^3.3` more;
* the total is dominated by *table construction*, not search, and a table
  write is counted as one rho step, which is generous to the table.

### 5.3.1 The structure, run for real

The optimum's shape — one orbit, two unknowns, two relations, a table built
once — is unusual enough to be worth exercising rather than trusting.
`validate_amortised_attack.py` builds exactly that table, streams targets
against it, solves, and checks `[d]P = Q`:

| `m` | `B` | `T` | `n` | table (built once) | relations | target attempts | `[d]P = Q` |
|---:|---:|---:|---:|---:|---:|---:|:--:|
| 13 | 13 | 1 | 3 | 20 | 2 | 5 | yes |
| 13 | 26 | 2 | 3 | 82 | 3 | 5 | yes |
| 13 | 13 | 1 | 4 | 98 | 2 | 3 | yes |
| 19 | 19 | 1 | 3 | 32 | 2 | 17 | yes |
| 19 | 38 | 2 | 3 | 140 | 3 | 13 | yes |
| 19 | 19 | 1 | 4 | 282 | 2 | 18 | yes |

Six of six recover the planted logarithm. The `T = 1` rows close a
two-unknown system from two relations, which is the shape the ECC2K-130
optimum uses, so that shape is not an artefact of the cost model.

**An eighth defect, in the script that produces that table.** Its first
form filtered the cells down to those carrying a recovery and then reported
"N of N" over the survivors. A support that could not be found, or a cell
that ran out of targets, returns without a recovery and so vanished from
the *denominator* — one real recovery beside one skipped support and one
outright failure reported as `all_recovered: true`, verdict "1 of 1".
Reproduced, then fixed: every cell is accounted for, the denominator is the
number of cells asked for, and the script exits non-zero if any cell did not
recover. Found by an external review agent on the pull request, not by this
thread.

That is four of nine defects now living in validation rather than in the
thing being validated. A ninth surfaced while fixing the eighth: the
counterfactual in `target_boundary.py` charged `attempts x max(build,
stream)` for rebuilding the table every time, which drops the build entirely
whenever streaming dominates — making the counterfactual *cheaper* than the
thing it exists to be worse than. Caught by this thread's own control
asserting amortisation is never dearer, which is the first time on this study
a test found a defect before a reviewer did. The pattern is stable enough to state as a finding of
its own: on this study, code written to check a result has been less reliable
than the code producing it, and the failure mode is always the same — the
check reports success over a subset it quietly chose.

### 5.4 The measured input, and a seventh defect

The model has one free parameter, the decomposition rate, so it is measured
rather than asserted — the same discipline that caught the `2x`–`43x` miss in
§4.5.

An earlier form of this model gave `sigma` a factor of 131 in that rate, on
the reasoning that a decomposition of `sigma^k(target)` is as useful as one
of `target`, so there are 131 acceptable right-hand sides. Measured on small
analogues, that model is optimistic by about **7x**:

| `m` | `T` | `B` | `n` | measured | `2^n C(B,n)/r` | ratio | with `sigma` factor | ratio |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 19 | 2 | 38 | 2 | 0.0117 | 0.0215 | 0.54 | 0.4082 | 0.03 |
| 19 | 3 | 57 | 2 | 0.0433 | 0.0488 | 0.89 | 0.9268 | 0.05 |
| 17 | 2 | 34 | 2 | 0.0650 | 0.0685 | 0.95 | 1.0000 | 0.07 |
| 17 | 3 | 51 | 2 | 0.1200 | 0.1558 | 0.77 | 1.0000 | 0.12 |
| 19 | 4 | 76 | 2 | 0.1033 | 0.0871 | 1.19 | 1.0000 | 0.10 |
| 17 | 4 | 68 | 2 | 0.2417 | 0.2783 | 0.87 | 1.0000 | 0.24 |
| 19 | 2 | 38 | 3 | 0.3633 | 0.5157 | 0.70 | 1.0000 | 0.36 |
| 17 | 2 | 34 | 3 | 0.7817 | 1.0000 | 0.78 | 1.0000 | 0.78 |

The reason is structural, and checked directly: **the set of `n`-subset sums
is itself `sigma`-closed** — true in every cell above. The support is
`sigma`-stable, so if a sum is reachable then so is every rotation of it.
Accepting `sigma^k(target)` is therefore one chance taken 131 times, not 131
independent chances.

So the Frobenius quotient buys memory (one canonical class per orbit) and
unknowns (`B/131` rather than `B`) — the latter is the large saving — and it
does **not** buy hit rate. That is the seventh defect in this thread's
running list, and like three of the first six it was in a model rather than
in code, and it was caught by measuring an assumption rather than by a test
failing.

The plain rate `2^n C(B,n)/r` holds to a mean of
`0.844` over 7 unsaturated cells —
slightly below 1, because subset sums collide — and that factor is carried
into the cost so the total is not quoted optimistically.

### 5.5 Why the meet-in-the-middle cannot be replaced

The remaining lever would be a better `k`-sum algorithm. Wagner's
generalised birthday would find an `n`-sum hitting a target in about
`k · 2^(129/(1+log2 k))` — `2^25.80` at `k = 16`, far below the reference.

It does not apply. Wagner needs partial matching: the low `t` bits of a sum
must depend only on the low `t` bits of the summands, so partial collisions
can be merged level by level. That holds for XOR and for addition mod `2^m`.
It fails for elliptic-curve addition, where `x(P+Q)` is not determined by any
truncation of `x(P)` and `x(Q)` — there is no prefix structure on curve
points to recurse on. The scalar domain would support it and the scalars are
exactly the unknowns; the abscissa domain supports it only through summation
polynomials, which is a Gröbner solve this repository has already measured
and closed (`RESEARCH_ECC2K130_DECOMPOSITION.md`).

So two lists are what is available, and a two-list meet-in-the-middle costs
the square root of the space it must cover. That square root, against a group
whose rho already takes its own `sqrt` with the same 262 automorphisms, is
the whole of the remaining gap.

### 5.6 Where this leaves the thread

`2^6.48` over the reference, with memory free, every input either
measured or derived, and the three directions that looked open — larger `m`,
a `k`-tree in place of the meet-in-the-middle, and rebuilding the table per
target — closed, corrected, and corrected again.

**One table, one unit**, as §2 of `AGENTS.md` asks. `S = operations / sqrt(r)`,
`sqrt(r) = 2^64.5`:

| variant | `log2` ops | `S` | vs rho | class |
|:--|---:|---:|---:|:--|
| Pollard rho, `<-1> x <pi>` *(reference)* | 60.81 | 0.077 | — | baseline |
| §1 as published, no quotient | 85.60 | 2,199,000 | `2^+24.79` | accounting |
| + Frobenius quotient at every `m` | 74.40 | 954 | `2^+13.59` | accounting |
| + priced as a logarithm, not a relation | 71.94 | 173 | `2^+11.13` | accounting |
| + stored table built once | 74.94 | 1.39e+03 | `2^+14.14` | accounting |
| **+ probes and `E[4]` translates priced** | **67.29** | **6.9** | **`2^+6.48`** | accounting |
| homogeneous relations, `m = 28` | 56.76 | 0.018 | `2^-4.05` | **relabelling** |

The last row is the one to read twice. It is *below* the reference, and it is
the only row in the table that is not an honest accounting of a logarithm —
it prices relations that determine nothing. §1's warning, drawn.

**The trajectory is the finding, and it is not about the curve.** Every step
above was a correction to this thread's own bookkeeping. Four of them made
the attack look better; the last made it look worse by `2^2.65`, and it was
found by an external review agent rather than by this thread. The remaining
`2^6.48` is small enough that one more error of the size of any of these
would move it materially in either direction. What it is not is evidence that
the curve is weak: the only lever anyone has applied here is arithmetic about
costs.

What would actually move it is a decomposition oracle cheaper than a
square-root search, and §5.5 closes the one generic candidate. Until such an
oracle exists, the floor is `sqrt` of a space that rho already square-roots
with the same 262 automorphisms.
