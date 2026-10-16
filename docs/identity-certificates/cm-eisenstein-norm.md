# Exact degree of a compact CM expression

Certificate `IDC1hf56b3716b0a7eb8b` is the evaluated supplement in
[cm-eisenstein-norm.json](cm-eisenstein-norm.json). It is checked over the specified
field of order 2^61-1, with no formal-proof pointer. The native producer and
independent verifier additionally expand and compare every integer coefficient,
using different integer implementations and sparse versus dense polynomial storage.
The exact check rejects a coefficient error of 2^61-1 even when the evaluated
certificate accepts it.

The identity, before geometric specialization, is

```
(a+b*t)(a-b-b*t) = a*a-a*b+b*b - b*b*(t*t+t+1).
```

On the exact SEC2 secp256k1 curve, a nontrivial cube root beta defines the
automorphism phi(x,y)=(beta*x,y). Endomorphism-ring cancellation gives
phi^2+phi+1=0; its dual is phi^2. Duality is additive and fixes integer
multiplication, so alpha=[a]+[b]phi has dual [a-b]-[b]phi. The exact polynomial
identity and the dual-degree theorem then give deg(alpha)=a^2-a*b+b^2 for the
original signed integers. The standard positive-characteristic facts and the
injectivity of integer multiplication in End(E) appear in
[MIT 18.783 Lecture 12](https://math.mit.edu/classes/18.783/2023/LectureNotes12.pdf).
The independent pre-implementation review is archived in
[Crypto Autoresearcher PR #2221](https://github.com/aburan28/crypto-autoresearcher/pull/2221),
at docs/reports/isogeny-cm-norm-design-20261016/independent-review.md.

The constructor uses ell=2^e-1 for an odd prime exponent in 3..=65536. It tests
exponent primality exactly, computes all e-2 Lucas-Lehmer transitions from s0=4,
and requires terminal zero. Its independent verifier repeats every transition
with repository limb integers and general modular reduction; the producer uses
num-bigint and Mersenne folding. The selected theorem is
[Conrad, Theorem 2.2](https://kconrad.math.uconn.edu/blurbs/ugradnumthy/lucas-lehmer-test.pdf).
Prime exponent 11 with composite M11 and composite exponent 9 are separate
negative controls. Exponent 2 is outside this implementation's documented domain.

Cornacchia supplies x,y with x^2+3*y^2=ell and unreduced coefficients a=x+y,
b=2*y. The verifier replays the normalized root, Euclidean divisions, stopping
condition, exact square test and coefficient transformation. It recomputes the
norm from the serialized integers. It also checks nonzero differential a+b*beta
mod p, full-point beta/lambda pairing, the nonzero rational-subgroup scalar and
its inverse, eight fixed public point images and both dual compositions.

The exact standardized p,n,G and cofactor are bound to the certificate and
checked against the registered source. Primality of p,n is explicitly trusted
from [SEC2 v2, section 2.4.1](https://www.secg.org/sec2-v2.pdf); this artifact does
not contain self-contained p,n primality certificates. Target-degree primality
has its separate independently replayed Lucas-Lehmer certificate.

The decisive alias control adds n to a. Rational-point outputs agree, while
the test recomputes both geometric norms from their unreduced integers and
checks a strictly positive difference. The old degree certificate is rejected.
Signed and zero-expression fixtures are checked through the native geometric
evaluator against independent affine arithmetic. The command's generated
certificate schema is canonical Cornacchia form; its subgroup evaluator is a
separately named optimization and never replaces a,b in the degree certificate.

Pure Frobenius and 1-Frobenius are theorem-level controls outside this canonical
Mersenne/Cornacchia executable schema. No executable acceptance/rejection of
those two inputs is asserted here. Their distinct issues are derivative zero
versus zero rational-point action; the native certificate independently checks
both differential and subgroup invertibility.

The compact certificate is a known-endomorphism self-map representation. It
contains neither expanded degree-sized rational functions nor a kernel
polynomial, and it earns no score in the expanded-map challenge. Implementation
controls do not instantiate the separately proposed M8191 scientific experiment.
