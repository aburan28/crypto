# Exact contract and work model

## C1. Cached square in characteristic two

Let K be any field of characteristic two, H a monic quintic, and
A=K[X]/(H). Write x for the class of X. The five powers 1,x,...,x^4
form a K-basis even when H has repeated roots. Define
C6=x^6 and C8=x^8 in A.

For r=r0+r1*x+r2*x^2+r3*x^3+r4*x^4, the exact square is
r0^2+r1^2*x^2+r2^2*x^4+r3^2*C6+r4^2*C8.
All cross terms vanish in a commutative ring of characteristic two.
The coefficients must be squared in K. No reducedness, separability,
squarefree hypothesis, inversions or root extraction is used.

Thus the cache preserves every intermediate vector and final support
predicate of the unmodified recurrence, including nilpotent elements.
It does not enumerate or recover original multivariate roots.

## C2. Exact setup

Write H=X^5+h4*X^4+h3*X^3+h2*X^2+h1*X+h0.
In A, x^5=h0+h1*x+h2*x^2+h3*x^3+h4*x^4.
Multiplication by x sends v=(v0,...,v4) to

(v4*h0, v0+v4*h1, v1+v4*h2, v2+v4*h3, v3+v4*h4).

Three applications starting at (h0,...,h4) give x^6,x^7,x^8.
Retain the first and third. The setup uses at most 15 coefficient
multiplications plus additions. Cache identity includes the complete
canonical H and field modulus/basis. Reusing it for a different H is invalid.
Nonmonic or different-degree input must use a proved fallback or report
unsupported; normalization, if used, is charged.

## C3. Domain support

For a binary independent basis v_i, beta_i=L_i(v_i) is nonzero and
L_(i+1)=L_i^2+beta_i L_i. Starting at X+a for an affine domain
a+span(v_i), the quotient recurrence is r_(i+1)=square(r_i)+beta_i*r_i.
C1 proves equality with literal polynomial squaring/reduction by induction.

For monic quintic H and squarefree L whose roots are precisely the requested
finite domain, all algebraic roots of H lie in that domain, retaining
multiplicities, iff H divides L^5. In characteristic two H divides L^8 is
equivalent. Three quotient squares test this nilpotence. If r=0, acceptance
is immediate. If r is nonzero and H squarefree, rejection is immediate.
Otherwise the full nilpotence test is required.

For prime fields or odd extensions the characteristic-two square formula
does not apply. A p-Frobenius cache needs a separate derivation and cost
model. In particular a prime-field interval is not a proper Fp-subspace.

## C4. Conditional dense operation model

S is one coefficient-field square and M one coefficient-field multiply.
Additions, allocation, parsing, serialization and hardware effects are
excluded from this symbolic model and must be included in wall measurements.

Use as reference a sparse characteristic-two square, not generic dense
polynomial multiplication: 5S and up to four monic reduction steps of 5M.
The stated dense model charges 5S+20M per quotient square, with no zero skips.
The cached model charges 5S+10M and at most 15M setup.

For 18 full recurrence steps and three tail squares, both charge 90M for
beta_i*r_i:
reference = 21*(5S+20M)+90M = 105S+510M;
cached = 15M+21*(5S+10M)+90M = 105S+315M.

For S,M nonnegative and M positive the model ratio is between 1 and
510/315=34/21. The upper endpoint is the limiting model with S=0.
This statement concerns this fixed dense count. Zero coefficients, sparse
moduli, different algorithms, cache traffic and actual field implementations
can change real costs. It is not a runtime bound or a universal lower bound.

The falsifiable prediction is a lower complete-support cost after paying
setup. A fast predicate alone does not establish that consequence.
