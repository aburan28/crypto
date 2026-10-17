# Correction to the prior transfer formula

Preserve all frozen versions and raw evidence. The transfer-only formula in
quintic_rational_support_20261017_v2/PROOF.md and
autolab_cycle_20261017_v2/TRANSFER_REVIEW.md omitted a term.
The correct expression is
V=(h^2+h)*e4+k*(z*s+z^2*r+k*yT+z^3+h*z*yT).
Consequently k=0 does NOT generally imply V=0. The simplification holds for
h=0 or1, but the actual recorded outer h values are neither.

Counterexample over GF4 with omega^2=omega+1: choose the valid curve
Y^2+XY=X^3+1, target(z,yT)=(1,0), h=omega,k=r=s=0.
The literal norm quotient is H=X^5+X+1, e1=0,e4=e5=1, so V=1.
The erroneous expression gives0. This refutes the claimed general k=0 transfer.

The exact generic factorization code, small oracles, and recorded support
predicate timings are unchanged: they compute U,V directly from H and never
use the erroneous transfer identity. The conditional fixed-U worker also
makes no assumption that actual k=0 cells have V=0. Do not promote the
claimed subfamily relevance; test the corrected source map independently.
This correction is not a new timing measurement or a novelty claim.
