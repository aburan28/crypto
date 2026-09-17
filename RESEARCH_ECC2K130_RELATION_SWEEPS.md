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

Measured here over six random normal elements of `F_2^131`: the weight-two
normal support is `sigma`-stable in all six, with no exceptions and no
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
lifts entirely or not at all. Measured over the same six bases, with zero
mixed orbits in any of them. A normal-basis weight-two support is therefore
always `131 * (number of lifting orbits)` points.

### 4.2 Support size is basis-dependent, over a wider range than previously recorded

An earlier note put the lifting-orbit count at `28` of `66`, giving the
`B = 3668` support this study has been quoting, and said the count "runs
28–36 across normal elements". Both the denominator and the range are wrong.
Six random normal elements, measured:

| `alpha` (leading bits) | lifting orbits (of 65) | `B` |
|:--|---:|---:|
| `0x2de8f4e94249b03fd0…` | 25 | 3275 |
| `0x4dea6516dc1ed48bd5…` | 29 | 3799 |
| `0x1a652f4bb357676c68…` | 31 | 4061 |
| `0x4d44a5aac26a00c420…` | 42 | 5502 |
| `0x57eea059506531b2fd…` | 26 | 3406 |
| `0x4f7aaf673fb46bcf47…` | 23 | 3013 |

The spread is `23–42` orbits, `B = 3013–5502`, not `28–36`. `B = 3668` is one
draw from this distribution and carries no special status. Nothing in the
conclusion turns on which draw is used — every one of them is astronomically
below the `B_4 = 2^33.4` of §1 — but a quoted support size should not be
mistaken for a property of the curve.

### 4.3 Frobenius acts as a scalar, so the support carries 28 unknowns, not 3668

`E` is a Koblitz curve, `#E(F_2) = 4`, so the base-field trace is `t = -1` and
`sigma` satisfies `sigma^2 + sigma + 2 = 0`. On the order-`r` subgroup this
has the root

    s = 196511074115861092422032515080945363956

verified here directly: `sigma(P) = [s]P` and `sigma(Q) = [s]Q` for the
challenge points, and `s^131 = 1 mod r`.

So for a support point `R` and its orbit, `log(sigma^k R) = s^k log(R)`. The
131 points of an orbit share **one** unknown. A support of 28 lifting orbits
presents 3668 points to a sweep and 28 unknowns to the linear algebra.

**This is what gives the negative its force.** The usual reason an index
calculus attempt fails inconclusively is that the linear algebra was never
reached — not enough relations to fill a `B`-column matrix. Here `B` is
effectively 28. About 29 relations would close the system, and the solve is
free at that size. The method still produces none, so the failure cannot be
attributed to an unaffordable second stage. It is a supply failure, and
supply is what §1 counts.

### 4.4 Five defects, each of which produced plausible output

Recorded because each one returned numbers that looked like results, and two
of them were in validation code rather than in the search — the code whose job
is to catch the other kind.

1. **Double-counted pairs.** Every unordered pair was enumerated twice, once
   anchored at each element's orbit, and the first orbit-grouped run reported
   125,953 "relations". Caught by the count itself being implausible, not by a
   test.
2. **Supply model never validated.** The expected-relation model was asserted,
   not checked. Its first test missed measured yields by `2x`–`43x`.
3. **Cancelling multisets counted as relations** — the cause of (2). The
   validation counter admitted formally cancelling multisets such as
   `R + (-R) + T + (-T)`. This is the construction-rank error of §3 reappearing
   *inside the code written to validate against it*.
4. **Detector too strict.** The search required the sum to be exactly `O` in
   the full group, when the attack only needs the sum to land in `E[4]`: the
   cofactor is 4, and a sum in `E[4]` still yields a usable equation on the
   order-`r` subgroup. At `m = 13`, 455 usable relations against 195 exact —
   57% of the supply was being discarded.
5. **Coset factor guessed.** The correction to (4) was first priced at `4x` by
   assuming the four cosets are hit uniformly. Measured, it is `2.33x`.

### 4.5 Status of the quantitative results

The corrected implementation and its run are being rebuilt in this branch.
The claims of §4.1–§4.3 are measured above and stand on their own. The
run-dependent figures — the measured expected-usable-relation rate, the
collision counts, the `m = 4` cost exclusion and the `2.33x` coset factor —
are restated from a prior session whose container was reclaimed before its
outputs were committed, and are **not** reproduced by anything in `results/`.
They are recorded here as prior observations pending re-measurement, and this
section is updated with measured values and frozen evidence files when the
rebuilt run lands.
