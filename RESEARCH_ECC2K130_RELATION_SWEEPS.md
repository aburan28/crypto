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

## 4. The three-point sweeps

Every row is a *complete* search over its support: all unordered pairs with
repetition, the third abscissa solved from `S3` rather than scanned. A zero
means the support has no homogeneous three-point relation at all.

| support | `B` | pairs | relations | `log2` predicted |
|:--|--:|--:|--:|--:|
| weight-two, complete | 4262 | 9,084,453 | 0 | −95.41 |
| weight ≤ two | 4336 | 9,402,616 | 1 (cofactor) | −95.34 |
| weight-two, σ-stable | 1090 | 594,595 | 0 | −101.31 |
| random matched #1 | 4262 | 9,084,453 | 0 | −95.41 |
| random matched #2 | 4262 | 9,084,453 | 0 | −95.41 |
| constructed, random `u,v` | 4000 | 8,002,000 | 0 | −95.69 |

Six supports, one relation between them, and it is in the cofactor. Measured
yield matches the counting prediction everywhere.

**The weight-two support is 4262 abscissae**, not the 3,668 an earlier
summary of this thread recorded; 8515 candidates `z^i + z^j` with `i < j`, of
which 4262 pass the trace condition and carry a point.

**The single hit** is `(0,1) + (1,0) + (1,0) = O`: the 2-torsion point plus
twice `(1,0)`. All three points are `F_2`-rational, none satisfies `[r]R = O`,
and the relation is the `Z/4` structure of `E(F_2)`. It lives entirely in the
cofactor and says nothing about `log_P(Q)`. It appears only because the
weight-at-most-two support admits `x = 0` and `x = 1`, which is the whole
difference between that row and the one above it.

**The constructed base with random coefficients is the control that matters.**
Points `[u]P + [v]Q` with uniform random `u, v` produce nothing, exactly as
counting predicts. Being built from `P` and `Q` is not by itself enough to
manufacture relations; the coefficients have to be small enough to make them.

## 5. Where structure actually wins, and what it wins

The four-point sweep is a collision search over pair sums: `P1+P2+P3+P4 = O`
means `P1+P2 = −(P3+P4)`, and a point and its negative share an abscissa, so a
relation forces two pairs whose sums collide. That is `B^2/2` work instead of
`B^4/24`, with one 64-bit digest stored per pair rather than the pair itself.

On the σ-stable weight-two support this found **88 relations where the
counting floor predicts 2^−89**. Structure does break the independence
assumption, by roughly `2^100`.

Every one of them is an identity of the endomorphism ring.

`#E(F_2) = 4` gives Koblitz trace `t = −1`, so Frobenius satisfies its
characteristic equation `σ² + σ + 2 = 0`, and therefore

    σ²(R) + σ(R) + [2]R = O

for **every** point `R` of the curve — verified here on random points and on
the challenge `P` and `Q` themselves. On a support closed under squaring that
is a free four-point relation with a repeated summand, one per base element
whose square and fourth power are also present.

It is not one shape. From `σ² = −σ − 2` comes `σ³ = 2 − σ`, hence
`2R − σ(R) − σ³(R) = O`, a different abscissa pattern and equally empty. A
first classifier matched only `{x, x, x², x⁴}` and caught 71 of the 88; the
remaining 17 were `{x, x, x², x⁸}`. The free relations are the whole ideal
generated by the characteristic equation, so matching shapes undercounts them
and invites reporting the remainder as a discovery.

The classifier therefore tests the definition rather than a pattern: read the
relation as a signed sum of Frobenius powers and apply that pattern to
independent random points. A pattern that annihilates every point of the curve
is an identity; one that annihilates only these points is a relation about
these points. Under that test:

| support | repeated-summand relations | Frobenius identities | genuine |
|:--|--:|--:|--:|
| weight-two, σ-stable | 88 | 88 | **0** |

**This closes the loophole the boundary left open.** In `(u,v)` coordinates a
Frobenius identity gives `(2 + λ + λ²)·(u,v) = (0,0)`, where `λ` is σ's
eigenvalue mod `r` and `λ² + λ + 2 ≡ 0`. Both `sum(u)` and `sum(v)` vanish:
these are construction rows. Structure can beat the counting floor, and what
it produces carries no information — the same accounting as a base
deliberately built from `P` and `Q`, arriving this time from the curve itself
rather than from how anyone chose the base.

## 6. The accounting result, measured

The small-coefficient constructed base is the positive control for §3, and it
is emphatic. Points `[u]P + [v]Q` with `u, v` drawn from a small range, base
size 4000:

| quantity | value |
|:--|--:|
| relations found (complete sweep, 8,002,000 pairs) | **2,045,254** |
| predicted by random counting | 2^−95.69 |
| verified on the curve (uniform random sample) | 20,000 of 20,000 |
| construction rows, `sum(u) = sum(v) = 0` | **20,000** |
| rows determining `log_P(Q)` | **0** |
| useful rank | **0** |

Two million relations. Every sampled one holds whatever `log_P(Q)` is. A
relation matrix built from this base would have enormous rank and would
determine nothing.

Set beside the two controls, the picture is complete:

| base | relations | useful rank |
|:--|--:|--:|
| constructed, small `u,v` | 2,045,254 | 0 |
| constructed, random `u,v` | 0 | 0 |
| weight-two, σ-stable (four-point) | 88, all Frobenius identities | 0 |

**Relation count is not evidence, and relation rank is not evidence.** Both
can be driven as high as one likes by choosing the support, and neither moves
`log_P(Q)` at all. The only quantity that means anything is the number of rows
whose `sum(v)` is a unit mod `r`, and across every experiment in this note
that number is zero.

This is why the `useful_rank` split exists rather than a rank count, and why
the thread's earlier framing — relation rank from `P,Q`-constructed bases —
needed the construction span quotiented out before any of it could be read as
progress.

## 7. Verdict

**Index calculus built on homogeneous point relations does not beat optimised
parallel Pollard rho on ECC2K-130, and the margin is not close.**

The reason is the counting floor of §1, not a failure to optimise. Against the
frozen reference of `2^60.809`:

- priced at a single relation, ignoring the `B` relations a logarithm needs,
  the best `m` is 4 at `2^65.79` — a factor `2^4.98` above rho, and asking
  `2^70.79` bytes of storage, about two zettabytes;
- priced honestly at `B` relations, the best `m` is 8 at `2^85.60`, a factor
  `2^24.79` above rho.

Class: **accounting**, by the §3 test. No `S` fell; a cost that earlier
framings of this thread had left in a column the headline did not look at has
been priced.

The falsification target of §1 was not met and its abandonment condition was.
Measured yields matched the counting prediction on every support where the
independence assumption holds, and the one family of supports that breaks the
assumption — σ-stable ones, via the Frobenius characteristic equation — breaks
it only into free relations.

**What this does and does not establish.** It closes homogeneous-relation
index calculus at `m ≤ 8` on this curve, including every support tried here:
low-weight, Frobenius-stable, random, and deliberately constructed from the
public `P` and `Q`. It says nothing about decomposition-based index calculus,
where the target is a known combination `[a]P + [b]Q` and the work is a
Semaev/Gröbner solve rather than a collision search; that family is priced
separately in `RESEARCH_ECC2K130_DECOMPOSITION.md` and on the scoreboard, and
it is not addressed by anything here. Nor does it threaten any deployed curve.

**What would reopen it.** A support whose relations are *not* consequences of
the endomorphism ring and *not* consequences of how the base was built, at a
rate above `B^m/(m! r)`. The `is_frobenius_identity` test in
`run_four_point.py` and the `useful_rank` split in `relations.py` are the two
screens any such claim has to pass, and both are cheap to apply.
