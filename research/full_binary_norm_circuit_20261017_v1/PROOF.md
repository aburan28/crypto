# Complete binary norm coefficient predicate

For E:Y^2+XY=X^3+aX^2+b, b!=0, T=(z,yT), z!=0, let
F=P+(hX+k)Y, P=X^3+rX^2+sX+t, t=z^3+rz^2+sz+(hz+k)yT.
Its literal norm N=P^2+XP(hX+k)+(hX+k)^2*C is divisible by X+z.
For H=N/(X+z), the following native representation is an exact predicate
for matching ALL six monic quintic coefficients, with no generic branch
assumption on h,k, repeated sources, zero sources or nonrational factors.
It is not a completed original S3-chain equivalence or an all-root solve.

Five source coordinates are d-bit polynomial-basis elements. Their monic
source polynomial is computed directly at true degrees and each coefficient
is reduced into K. Write its coefficients(e5,e4,e3,e2,e1,1).
The norm identity forces h^2+h=e1+z. The binary Artin-Schreier map has
kernel{0,1} and rank n-1. A fixed linear section and its one-dimensional
residual provide h=section(e1+z)+eta, eta in F2, together with the exact
residual-zero equation. No h roots are discarded, even in even extensions.
The section/residual are precomputed and checked on every polynomial basis
vector by literal field squaring; no half-trace parity shortcut is used.

Set beta=sqrt(b), c=(yT+beta)/z. The curve identity gives
c^2=yT/z+a+z. Then coefficient e2 uniquely determines
k=e2+r^2+hr+a*h^2+z*e1.
Coefficient e5 uniquely determines, without any variable inverse,
w=sqrt(e5/z), s=zr+c*k+yT*h+z^2+w.
Perfectness of finite fields makes square roots unique. The complete two
remaining coefficient conditions are
R3=e3+h*s+k*r+k^2+z*e2=0,
R4=beta*k^2+(z*w)*k+e5+z*e4=0.
The last identity follows by substituting s and using c^2 above. Thus the
residual equation plus R3/R4 are necessary AND sufficient for all norm
coefficients; k and s are not omitted independent witnesses.

Every primary assignment is five ordered d-bit sources, n-bit r, one eta.
For n83,d18 this is174 bits. The predicate includes every outer h,k value
consistent with the source-sum relation, not just four sampled fixed cells.
It still searches only the norm chart; original polynomial roots may have
additional nonrational/exceptional branches. Complete original S3 coverage,
all source orderings/auxiliary fibers, and invocation-bound ghost/exception
closure remain required before a global all-roots claim. No source or
original exact count is produced by this compile-only worker.

All Frobenius squarings/square roots, constant multiplication and AS sections
are linear Boolean circuits. Source products and the four variable field
products use the same selected arithmetic cutoff. Both arms include every
norm-generation operation. Reachability removes dead intermediate gates
symmetrically before counting. This differs from earlier subsystem counts;
do not multiply gains or compare incompatible totals. Each arm uses the
same standard3-AND/4-XOR/unit Tseitin model and a separate gate arena.

Independent controls exhaust73216 primary assignments on8 small curve/target
models over GF4/GF8/GF16: both target signs, a0/a1 where rational, c=0,
zero/repeated sources and even extension degree. They evaluate the generated
predicate and compare to a separately evaluated literal norm and exact
division by X+z. They do not reuse eliminated residual formulas. Both cutoff
arms pass. This is mathematical/implementation evidence, not Lean or external
independent review; Karatsuba and norm identities are classical prior art.

Prediction: on the actual frozen n83 fixture cutoff8 reduces reachable
AND gates>=25percent and standard clauses>=15percent relative to the full
schoolbook predicate. Include parsing/constants/field verification/compile
in diagnostic cold cost, no solver call or output roots.90seconds/4GiB/oneCPU;
final serialization checked and late profile flags demoted; supervising full
process lifetime also charged. No qualified8x or matched F4/F5/F6 comparison.
