# Added exact exceptional and zero-root branches

All v1 proof obligations remain; frozen v1 is unchanged. This version replaces
the U=V=0 fallback by H=(X+e1)*B(X^2), B(T)=T^2+e2*T+e4.
Every algebraic root of H lies in D iff e1 lies in D and every algebraic
root of B lies in D^2. Frobenius is a bijection on the algebraic closure,
and on the finite field; this equivalence retains all multiplicities and
nonrational quadratic factors. Test B with the squared annihilator coefficients
beta_i^2 and squared affine offset a^2, using the exact quadratic nilpotence
criterion. No square roots are computed.

For Delta=0, U!=0, V=0, t=0. Then e4=0 and e5=0, so H=X^2*Q.
All-support iff 0 belongs to D and Q has all roots in D. The cubic is
obtained by dropping the first two zero coefficients; no inversion or
rational-membership ladder is needed. Repeated-zero-outside is rejected
only as all-support, never as an empty original root set.

Source-map transfer, derived from research/advance/five_divisor.py:
V=e5+e1*e4=k*(z*s+z^2*r+k*yT+z^3+h*z*yT).
For fixed h,k, both U and V are F2-affine in r,s. In particular k=0
forces V=0; thus the zero branch has an actual structural use. This alone
does not quantify its frequency across all outer parameters or solve the
norm-family search. Next measure U-map ranks on the actual73-dimensional
relaxed affine family before proposing fixed-U partitioning.

The v1 comparison showed a component signal but GF409 A/A exceeded10%;
no timing qualification is inherited. This is a distinct frozen candidate
and seed. Five internal rounds and A/A are diagnostic; complete original
all-roots matched cold F4/F5/F6 and external-fastest comparisons remain open.
