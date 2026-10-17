# Bounded finite-field reduction diagnostic

Frozen before execution. This diagnostic targets the original C model, not a
replacement curve: p=1637958981148213249, a=571607911220311176,
b=29329524900995129, N=1637958981148213248. Native complete factorization
gives N=2^26*3^20*7 and largest prime r=7, with ord_r(p)=1. The 61-bit field
does not make this a 61-bit prime-order DLP. This is one model, not C's entire
1,934,901,248-vertex class, its other prime-power components or all strata.

Implement standard Miller evaluation and the reduced Tate pairing over F_p,
using the shifted divisor (Q+S)-(S) and final exponent (p-1)/r. Require a
nontrivial r-th root image of the generator; reject a degenerate image or
exhausted shift search. No pseudo-pairing or generic DLP fallback is allowed.
Test the p=113,r=7 control and all six nonzero control/C scalar inputs against
independently constructed curve multiples. These finite tests are not a new
bilinearity theorem or a proof for untested curves.

The finite-field destination method uses two distinct inverse-paired images
of small integer primes raised to (p-1)/r. Its unknowns are both factor-base
logs and the target log relative to the pairing image. Collect checked
two-summand multiplicative relations beta*alpha^s=product(z_i^e_i) and solve
the resulting system modulo r with native incremental Gauss-Jordan. Do not
provide any factor-base log or planted target scalar to the solver. Verify
all recovered base logs, the destination answer and [x]G=Q. Rank failure is
exhaustion, not permission to invoke rho/BSGS/enumeration.

This is a projected prime-subgroup factor-base relation method, not NFS or
full multiplicative-group integer-smoothness index calculus. r=7 is an
implementation diagnostic, with no scaling, sub-rho or security claim.

Run two cold public targets, made by independently searching and cofactor
projecting points, on original C. Include setup, projection, pairing, relation
collection, linear algebra and verification in end-to-end elapsed time.
Record phase operation ledgers, actual query counts, candidate pair tests,
relations, rank, missing/unpriced work and unqualified elapsed observations.
No admitted runtime comparison is planned; builds/tests may precede this
diagnostic and finish before the resumed survey. Preserve every failure.
