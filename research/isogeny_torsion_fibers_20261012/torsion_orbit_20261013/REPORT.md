# Rational three-torsion orbit relation-support screen

2026-10-13 local / 2026-10-14 UTC, Conductor T-2611. This is the exact
combinatorial screen frozen in [PROTOCOL.md](PROTOCOL.md), following the
[four-size cold solve](../full_cost_20261013/REPORT.md). The [raw TSV](raw-v1.tsv),
[native source](screen.rs), [run log](run-v1.log), and [run receipt](RUN.md)
are retained. The [artifact manifest](ARTIFACTS.sha256) hashes the study
records. [The graph](graph-support.svg) displays the exact support
counts; [the PDF](report.pdf) embeds it.

An independently frozen [follow-on protocol](EQUIVARIANT_PROTOCOL.md)
tested a distinct degree-three coordinate. Its [raw counts](equivariant-v1.tsv),
[native verifier](equivariant.rs), [run log](equivariant-v1.log), and
[orbit diagram](equivariant.svg) are also part of this result. The
[follow-on run report](reports/run-20261014T0551Z/README.md) records the
same-host replay and its retained first infrastructure failure.

## Question and result

For the four previously measured non-isomorphic degree-three pairs, we
compared two-sum support of 12 subgroup points with all 36 points obtained
by adding rational order-three torsion. The orbit arm allowed the known
torsion correction, so a relation was counted whenever it legitimately
returned to the prime-order subgroup. We also tested 36 *distinct* subgroup
points selected by the same factor-base scan.

| `p` | subgroup `r` | arm | rational representatives | log unknowns | raw pairs | distinct reachable subgroup elements | exact hit probability |
| ---: | ---: | :--- | ---: | ---: | ---: | ---: | ---: |
| 617 | 103 | 12 subgroup | 12 | 12 | 78 | 64 | 64/103 = 0.621359 |
| 617 | 103 | 36 torsion lifts | 36 | 12 | 666 | 64 | 64/103 = 0.621359 |
| 617 | 103 | 36 distinct subgroup | 36 | 36 | 666 | 103 | 103/103 = 1.000000 |
| 2,633 | 439 | 12 subgroup | 12 | 12 | 78 | 74 | 74/439 = 0.168565 |
| 2,633 | 439 | 36 torsion lifts | 36 | 12 | 666 | 74 | 74/439 = 0.168565 |
| 2,633 | 439 | 36 distinct subgroup | 36 | 36 | 666 | 351 | 351/439 = 0.799544 |
| 9,677 | 1,613 | 12 subgroup | 12 | 12 | 78 | 75 | 75/1613 = 0.046497 |
| 9,677 | 1,613 | 36 torsion lifts | 36 | 12 | 666 | 75 | 75/1613 = 0.046497 |
| 9,677 | 1,613 | 36 distinct subgroup | 36 | 36 | 666 | 559 | 559/1613 = 0.346559 |
| 38,561 | 6,427 | 12 subgroup | 12 | 12 | 78 | 78 | 78/6427 = 0.012136 |
| 38,561 | 6,427 | 36 torsion lifts | 36 | 12 | 666 | 78 | 78/6427 = 0.012136 |
| 38,561 | 6,427 | 36 distinct subgroup | 36 | 36 | 666 | 642 | 642/6427 = 0.099891 |

The ratio of orbit support to 12-point support is exactly `1.000` at all
four sizes. The 36-distinct/36-lift support ratios are `1.609`, `4.743`,
`7.453`, and `8.231`. These are **one-probe relation-hit ratios**, not
complete-solve speed ratios. Orbit construction, scalar multiplication,
row work, and memory were not timed in this screen.

The source code makes 78 group-addition calls to form the 12-point pair
table, 1,998 calls to form and torsion-correct the 666 orbit pairs, and
666 calls for the 36-distinct table, per curve. Constructing the 36 orbit
representatives makes another 36 addition calls as implemented. These
semantic counts exclude base selection and validation arithmetic; they
are not a complete cost unit. The 36 orbit representatives collapse to
twelve quotient classes on every curve.

## Mathematical explanation and proof status

Let `K=ker(phi)={O,T,2T}` and `H=<G>` have prime order `r` coprime to three.
For any `A,B∈H` and `t,u∈K`, homomorphism gives

`phi((A+t)+(B+u)) = phi(A+B)`.

Likewise `[6](A+t)=[6]A`, because `[6]K=O`. After correcting the known
torsion residue, `(A+t)+(B+u)−(t+u)=A+B`. Thus the 666 raw pairs made from
36 torsion lifts have the same quotient two-sum support as the 78 raw pairs
from their twelve orbit classes. Extra representations do not add distinct
subgroup targets. Since `gcd(3,r)=1`, a fiber has exactly one point in `H`:
if two `H` points differ by a point of `K`, that difference has order
dividing both `r` and three and is identity.

The native program checked the exact torsion order, all 36 lifts per curve,
their distinctness, both source and destination curve equations, the
isogeny invariance, cofactor invariance, and all subgroup orders. It
enumerated the full `r`-element scalar index and matched point-addition
support to modular-scalar support. The conditional quotient-support
implication is compiled as
[`torsion_orbit_pair_support`](lean/OrbitSupport.lean) with no axioms;
see [proof map](PROOF_MAP.md). The Lean module assumes the curve map is a
homomorphism and that `K` maps to identity. A formal curve-level bridge and
the unique-subgroup-lift theorem remain open; native checks cover only these
four finite curves.

## Interpretation and next discriminating test

The free three-lift yield multiplier is false for this pair-table collector.
This does **not** refute torsion symmetry as a way to reduce a decomposition
polynomial or its linear algebra. Faugère, Huot, Joux, Renault and Vitse
developed that different mechanism for extension-field index calculus in
[EUROCRYPT 2014](https://www.iacr.org/archive/eurocrypt2014/84410158/84410158.pdf).
It changes the oracle representation and factor-base quotient; the cost of
the oracle must be measured against a same-target source system. The
[next gate](NEXT_GATE.md) specifies the experiment and its controls. The
current prime-field scan has no nontrivial subfield factor-base condition,
so it cannot test the polynomial mechanism.

This result corrects the unconditional `×|K|` collection-gain reading in
[`RESEARCH_ISOGENY_FIBERS_ECDLP.md` §4.4–4.5](../../../RESEARCH_ISOGENY_FIBERS_ECDLP.md):
quotienting a torsion-closed base reduces duplicate log columns, but a
same-size *independent* base can have much higher relation-hit probability.
A full comparison must hold the input target and total method cost fixed.
The exact support screen has no `S` or rho speed ratio. The standing generic
two-class floor is `sqrt(pi/4)=0.886227` in whole-method `S`; the measured
source negation rho rows remain in the prior [full-cost report](../full_cost_20261013/REPORT.md).

## A distinct equivariant coordinate

The quotient-invariant coordinate `x∘phi` has the same value on all three
points of a fiber. The follow-on found an explicit **nonfactorable**
degree-three coordinate

`chi(P)=(y(P)+1)/(y(P)−1)` and `f(z)=1/(1−z)` on `P¹`,

with `chi(P+T)=f(chi(P))` for `T=(0,1)`. At the kernel points,
`chi(O)=1`, `chi(T)=infinity`, and `chi(2T)=0`, while `phi` sends all three
to identity. Thus `chi` does not factor through `phi`, and `f` has order
three. The native verifier exhausted all `p+1` rational points for each
of the four fields: 618, 2,634, 9,678 and 38,562 points. It checked the
`+T` and `+2T` identities on each, for 1,236, 5,268, 19,356 and 77,124
coordinate identities, and checked `f³=id` on every `P¹(F_p)` element.
All assertions passed. These are finite checks of a proposed identity,
not an extension-field decomposition or cost measurement.

For the algebraic proof away from the three exceptional points, set
`u=(y+1)/x` and `v=(y−1)/x`. The curve equation gives `uv=x`; the slope
for adding `T` is `v`. Substitution into the group law gives
`u(P+T)=−v(P)` and `v(P+T)=u(P)−v(P)`. Taking the ratio yields
`chi(P+T)=−v/(u−v)=1/(1−chi(P))`. The kernel values above check the
exceptional points directly. The transformation `(u,v)↦(−v,u−v)` has
order three. Because the coordinate is a Möbius transform of `y`, which
has degree three on this Weierstrass curve, `chi` also has degree three.
The argument uses characteristic different from two and three; at
`O,T,2T`, the order-three group law supplies the projective cycle
`1→infinity→0→1` directly.
The [proof map](PROOF_MAP.md) records the formalized order-three
consequence and the remaining curve-level Lean bridge.

An equivalent generic calculation makes the denominator conditions
explicit. For `x≠0`, the addition slope is `(y−1)/x`, so the group law
gives `y(P+T)=y−2−(y−1)³/x³`. On `y²=x³+1`, with `y≠−1`, this reduces to
`g(y)=(y−3)/(y+1)`. The denominator-cleared reduction is the polynomial
identity in [IDC1h41adc1f4e1c5522c](group-law.identity.certificate.json);
the conjugacy `chi∘g=f∘chi` is certified, after clearing denominators,
by [IDC1ha0c3d83775e32754](mobius-conjugacy.identity.certificate.json).
Each `identity.certificate/v1` record passed its 64 sealed and 64 extra
checker-seed evaluations using the native
[issuer and checker](../../../examples/isogeny_torsion_identity.rs).
The first certificate's first three terms simplify over the integers to
`(y−1)²(x³−y²+1)`, and its final term is their negative. Both sides of
the conjugacy certificate simplify to `−4(y−1)`. These expansions give
an exact written polynomial proof in addition to the finite-field checks.
These records certify polynomial steps over the certificate field;
the written curve-law argument supplies their curve assumptions, and the
three projective kernel cases were checked separately. Neither certificate
measures a decomposition oracle.

This is a candidate for a torsion-equivariant *polynomial oracle* over an
extension field. The [2014 primary paper](https://www.iacr.org/archive/eurocrypt2014/84410158/84410158.pdf)
distinguishes nontrivial equivariance from coordinates factoring through
an isogeny. It also notes that quotienting a torsion orbit can decrease
decomposition probability; any polynomial-cost improvement must outweigh
that loss and all recovery costs. The [next gate](NEXT_GATE.md) uses `chi`
as the candidate and `x∘phi` as a control.

## Requirement status

| Requirement | Status | Evidence or remaining gap |
| :--- | :--- | :--- |
| Non-isomorphic source/destination paired relation solves and matched rho | verified for the earlier four-size public workload | [Full-cost report](../full_cost_20261013/REPORT.md): 640/640 measured arm scalars verified; local cold CPU only. |
| Rational three-torsion lift relation-support effect | verified for the exact pair-table screen | Twelve rows in [raw-v1.tsv](raw-v1.tsv), complete subgroup scalar-index cross-check, same-host replay, and conditional Lean theorem. |
| Nonfactorable cubic equivariant coordinate | algebraically derived and finite-instance verified; full Lean bridge partial | [Raw identity sweep](equivariant-v1.tsv), two [polynomial certificates](CERTIFICATES.md), [proof map](PROOF_MAP.md); projective curve theorem remains to formalize. |
| Extension-field decomposition, relation rank, target recovery, and fully charged matched-rho comparison for `chi` | not attempted | Construct the native four-arm oracle and execute the [pending gate](NEXT_GATE.md). |
| Portable `S`, isolated wall admission, and independent-host replay | pending | The previous full-cost round is shared-host L0; no new timing was measured in this screen. |

The study-local support graph and equivariance diagram are new. The existing cold-CPU graph
and the canonical index-calculus scoreboard were checked and left unchanged:
this screen adds no complete-solve timing, scalar-recovery arm, or portable
counted `S` value.
