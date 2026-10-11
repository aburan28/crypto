# Contribution relative to Joux–Vitse

The contribution developed here is an arithmetic restriction and a reproducible
class analysis for the **norm-one cubic subfamily**. The cover-and-decomposition
algorithm and the broader weak-cover families belong to earlier work.

Joux and Vitse's [EUROCRYPT 2012 paper, §4.1](https://www.iacr.org/archive/eurocrypt2012/72370010/72370010.pdf)
combines Weil descent with decomposition index calculus. It treats
`y²=h(x)(x−alpha)(x−alpha^q)` with **h of degree one or two over F_q**.
Its isogeny-reach discussion explicitly assumes independence of weak form and
isogeny class, and conjectures reach for curves whose cardinalities are divisible
by four. The underlying weak-cover classification is credited there to prior
work, including Thériault and Iijima–Momose–Chao.

## A scope correction verified by exact counts

Our original census enumerates `h=x`, up to the recorded translations, square
scalings and twists. It does not enumerate the nonsplit quadratic h branch.
The conductor obstruction and the cubic zero labels remain correct in that
domain. They cannot be read as an exclusion of every family in Joux–Vitse.

We tested 100 constructed quartics
`y²=(x²−d)(x−alpha)(x−alpha^49)` over `F_(7^6)`, with d a nonsquare in
F_49 and alpha outside F_49. The seed, exact modulus, d, and all 100 parameters
are retained in [the scope check](quadratic_scope_check.txt).
Of these published-family quartics, **46 have ordinary traces excluded by the
cubic norm-one predicate**. Three separate all-x enumerations verify the
Jacobian point counts and an explicit quartic-to-Weierstrass change of coordinates:

| Trace | Exact cardinality | Frobenius conductor | Cubic-family status | Quadratic-family status |
| ---: | ---: | ---: | --- | --- |
| −38 | 117,688 | 18 | Excluded by conductor depth one | Explicit model, exact count |
| −10 | 117,660 | 26 | Excluded by conductor depth one | Explicit model, exact count |
| 610 | 117,040 | 72 | Zero in the complete cubic census | Explicit model, exact count |

In particular, **trace ±610 is a cubic-family zero, while the broader
Joux–Vitse family has a quadratic representative at this trace pair**.
The opposite sign is supplied by a nonsquare base-field quadratic twist,
which retains h in F_49[x]. This example therefore does not refute their
conjecture for the full family.

For the explicit model, let f(x) be the displayed monic quartic, let
`c3=f'(alpha)`, `c2=f''(alpha)/2`, and `c1=f'''(alpha)/6`. The map
`X=c3/(x−alpha), Y=c3*y/(x−alpha)²` gives
`Y²=X³+c2*X²+c3*c1*X+c3²`; its inverse is
`x=alpha+c3/X, y=c3*Y/X²`. The source has the rational point `(alpha,0)`
and two points at infinity. Each independent enumeration includes all affine
points and those two infinity points. Three additional mapped-point checks
per model verify the equation and inverse.

## What belongs in the contribution statement

| Work component | Appropriate assessment |
| --- | --- |
| Genus-3 cover and decomposition method | Reproduction of the published method; credit Joux–Vitse and the underlying cover literature |
| Native solver and enumeration implementation | Engineering and measured reproduction; any speed claim requires its own matched complete-pipeline evidence |
| Norm-one conductor restriction | Proved restriction for the cubic subfamily, stronger than the order-divisible-by-four observation in their §4.1; priority against the wider literature is still to be checked |
| Exact cubic class census and certified zeros | New project evidence about the stated subfamily, with exact/probabilistic point-count status retained |
| Absolute-Frobenius orbit formula | A proved enumeration formula and reproducible reduction in required point-count calls |
| 192–252-bit population | New arithmetic and bounded-search measurements; 63.1% is cubic admission, with passing cubic class labels unresolved |
| One-edge full-4-torsion result and density bound | Explanatory arithmetic using standard torsion/isogeny facts; the direct-density order 1/q is already discussed in the prior paper |
| Exact CM-intersection track | A candidate stronger class-decision contribution; the present large-field support preflight fails at the recorded library integer limit |

The earlier [Auer–Top theorem](https://pure.rug.nl/ws/portalfiles/portal/14410242/2002JNumberThAuer.pdf)
addresses isogeny to arbitrary Legendre curves, with the stated extremal
exception; it does not select our norm-one subfamily. Its Lemma 2.1 and
Proposition 2.1 already describe rational halving and full 4-torsion in terms
of root differences. These elementary ingredients should be credited rather
than presented as newly invented techniques.

A defensible contribution statement is: **We prove and measure arithmetic
class restrictions for the norm-one cubic branch of the Joux–Vitse setting,
derive exact Frobenius-orbit census formulas, and distinguish that branch from
the nonsplit quadratic branch by independently counted representatives.**

The next substantial step is a complete two-branch class criterion and an
explicit construction with a verified cost on independently sampled classes.
That will address the full-family reach question rather than assigning the
cubic subset's admission rate to all published weak covers.
