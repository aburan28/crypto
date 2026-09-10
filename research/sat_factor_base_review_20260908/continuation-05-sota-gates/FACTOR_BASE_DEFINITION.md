# Algebraic factor-base definition and proof boundaries

Reviewed on 10 September 2026 against crypto commit
`6289769cbb781b242275f45c76803ef28e0ec2c8`. This note supplies an explicit
mathematical definition, proof obligations, and interpretation corrections.
It records no new solver run, benchmark, external review, or novelty result.

## Definition and exact kernel dimension

Let K = F_(2^n), let sigma(x) = x^2, and take a nonzero monic divisor
f(T) = sum_i f_i T^i of T^n - 1 in F_2[T]. Define

\[
L_f(X)=\sum_i f_iX^{2^i},\qquad V_f=\{x\in K:L_f(x)=0\}.
\]

The following is a direct algebraic proof. Each squaring power is F_2-linear,
so V_f is an additive F_2-subspace. Since the coefficients are in F_2,
L_f(x^2) = L_f(x)^2. Squaring therefore maps V_f into itself, and its
bijectivity on a finite field makes this an equality.

Write T^n - 1 = f(T)g(T). Composition of linearized polynomials gives

\[
L_g\circ L_f=X^{2^n}-X.
\]

Every root of L_f in an algebraic closure consequently lies in K. The
constant coefficient of f is 1, so L_f has derivative 1. Its ordinary degree
is 2^(deg f), and all those roots are distinct. Thus

\[
|V_f|=2^{\deg f},\qquad \dim_{\mathbb F_2}V_f=\deg f.
\]

This proof works for even n as well. For odd n, T^n - 1 is squarefree, allowing
the familiar subset-of-distinct-factors description. For even n the divisor
multiplicities must be retained; that mathematical extension does not certify
the implementation's distinct-factor enumeration for even n.

The degree of f is an additive-space dimension. It is different from the
ordinary degree of L_f, the degree of a field-representation polynomial, and
the field degree of an individual root. V_f need not be a subfield. The edge
cases f=1 and f=T^n-1 give {0} and K respectively; the reviewed constructor
intentionally rejects these endpoint dimensions.

The divisor/linearized-coordinate construction is established prior art.
GGMP section 4.1 describes it, and Lemma 4.1 explicitly assumes that n is an
odd prime for the equal-degree assertion about nontrivial irreducible
factors. [Primary conference paper](https://sacworkshop.org/SAC20/files/preproceedings/18-IndexCalculus.pdf).
The proof and implementation distinctions here are explanatory analysis,
not a proposed new construction.

## Rational points and public group maps

For the nonsingular Koblitz model

\[
E_a:y^2+xy=x^3+ax^2+1,\qquad a\in\mathbb F_2,
\]

define the affine factor base

\[
B_f=\{(x,y)\in E_a(K):L_f(x)=0\}.
\]

The identity O has no affine x-coordinate and is excluded. Coordinate
squaring preserves the curve equation and the coordinate predicate.
Negation, (x,y) -> (x,y+x), preserves the predicate as well. Hence B_f is
stable under Frobenius and negation. Neither property makes B_f an
elliptic-curve subgroup.

The coordinate count does not equal the rational-point count. The zero
coordinate has the unique lift (0,1); each nonzero coordinate has either zero
or two rational lifts. With d=deg f,

\[
|B_f|=1+2\,\#\{x\in V_f\setminus\{0\}:x\text{ has a rational lift}\}
\le 2^{d+1}-1.
\]

Using |B_f| approximately 2^d requires a density assumption or measurement.
It is not the exact kernel-dimension theorem.

Suppose a subgroup G has order r and index h with gcd(h,r)=1. Multiplication
by h kills the quotient E_a(K)/G and is invertible on G, so it maps the full
curve group onto G. It commutes with Frobenius and negation. On the restricted
set B_f it can create collisions and the identity; the resulting coordinates
need not satisfy the original L_f predicate. Calling this operation
"projection" is shorthand for cofactor clearing: it is not generally the
identity on G or an idempotent map.

These public group operations do not provide scalar labels for factor-base
points relative to a generator. Likewise,

\[
L_f(x(P))=0
\]

is an equation in the additive field. It does not imply f(pi)P=O under
elliptic-curve addition. Coordinate-linear algebra and the point-group law
must not be interchanged.

## Correspondence to the reviewed code

- `subspace_basis_for_divisor` multiplies the selected ordinary polynomial
  factors, constructs the associated linearized exponents, and rejects a
  kernel whose returned dimension differs from the ordinary divisor degree.
- `build_frobenius_factor_base_from_divisor` passes that basis into rational
  point materialization. Enumerating the selected factor-base coordinates
  is different from enumerating all elements of the target subgroup, and
  its cost still has to be charged.
- `koblitz_public_factor_base_discovery` receives n, a, m, and a requested
  dimension. Its candidate scoring uses public point counts, cofactor-class
  reachability, and cofactor-cleared signed-orbit counts. The inspected
  selection path has no target input, relation-yield input, or solver-timing
  input. Its public curve setup also performs point search and group-order
  work; the abstract definition does not make those implementation costs
  disappear.
- Output fields saying that forbidden inputs were not used are not by
  themselves a proof. The supporting evidence is the scoped source review
  and the pinned source/artifact provenance. This note is not a formal
  whole-program information-flow proof.
- For a target R in G, a decomposition R=P_1+...+P_m with P_i in B_f,
  allowing repeated points, requires sum_i [r]P_i=O. The cofactor gate tests
  whether O is reachable with exactly m such classes. Failure rules out
  that decomposition; success does not certify any prescribed target, useful
  relation yield, independent relations, matrix rank, or a runtime bound.

## Scope of the result

This proof explains the existing algebraic construction. It does not establish
a good end-to-end factor base, efficient equation solving, a complete benchmark,
licensed Magma resources, or an unaffiliated novelty assessment. Those are
separate evidence gates. No executable behavior changes accompany this note.
