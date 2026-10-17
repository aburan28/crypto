# Exact actual-family rank profile

Input contract: the source-identified GF(2^83) fixture with modulus
X^83+X^45+X^2+X+1, zero source offset, independent dimension18 source basis,
curve Y^2+XY=X^3+a X^2+b and nonzero target abscissa z. This component
returns exact compact affine profiles, not original roots or source tuples.
It performs no outer parameter enumeration and no search of2^73 candidates.

The literal norm of P(X)+(h X+k)Y is P^2+X P(h X+k)+(h X+k)^2*C(X),
C=X^3+aX^2+b. Choosing P(z)+(h z+k)yT=0 forces division by X+z.
Direct native polynomial arithmetic and exact division independently check
all five elementary-coefficient formulas on2304 small inputs and on a
particular plus every kernel direction of each actual affine family.
The native port initially omitted z^2*k from e4; its failing preparation
log is retained, and the corrected port passed literal norms before freezing.

At fixed h,k each coefficient is F2-affine in r,s because Frobenius squaring
is F2-linear. Product spans W_j are constructed by actual field products
and Gaussian basis reduction, not assumed from bit bounds. The source
polynomial's j-th elementary coefficient must lie in W_j. Enforce e1 in
W1 as an outer condition; concatenate linear quotient residuals for e2..e5.
Column elimination gives every solution (r,s) as a particular mask plus an
independent kernel basis. Tiny independent controls enumerate every 4-bit
(r,s) mask and compare direct literal norms against this affine solution set.
For the four declared actual h,k profiles, independent rank93/nullity73 and
span dimensions18,35,52,69,83 must be reproduced or the worker fails.

U=e3+e1*e2 and V=e5+e1*e4 are themselves F2-affine on a fixed-h,k family.
Evaluate their differences on every kernel generator; ordinary binary ranks
are the exact U-image, V-image and joint-image dimensions. A fixed-U slice
has dimension73-rank(U). This does not count roots or satisfy Delta=0.
The proposed bounded partition policy admits at most2^8 U-values; if U
rank exceeds8 the worker refuses this global extrapolation and retains the
rank evidence for a different variable-U incidence formulation.

Correct transfer identities:
U=h*s+k*r+k^2+(h^2+h)*e2,
V=(h^2+h)*e4+k*(z*s+z^2*r+k*yT+z^3+h*z*yT).
The second expression includes the missing term identified in ERRATUM.md.
The native small oracle explicitly refutes the old k=0 simplification.
No prior frozen source, input, report or measurement is edited.

Work is polynomial-size field arithmetic, at most166 original affine columns,
73 mapped kernel directions per admitted profile, bitset elimination and small
polynomial products; storage polynomial in n. All actual family basis points
are replayed against the literal norm. There is no all-parameter/root scan.
Review is one-model with independently formulated arithmetic, not Lean or
independent-human review. Classical Gaussian elimination is not novel.
90seconds/4GiB/oneCPU; final serialization is checked and late completions
are demoted. The supervisor additionally charges the complete child lifetime.
