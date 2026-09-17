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
