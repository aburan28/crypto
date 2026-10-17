# Exact repeated Q-root branch on the actual family

Retain the corrected literal-norm identities and affine-family proof from
actual_norm_u_rank_20261017_v1. Work over the actual GF(2^83) fixture and
four declared fixed h,k values only. Source offset is0 and z is outside D.
This is not original multivariate-root or full norm-family completion.

Let Q(X)=hX+k, h!=0, x=k/h. If x is outside D the branch with a repeated
source root at Q(x)=0 has no all-support candidate. If x belongs to D,
x!=z by the admitted target condition. For every r,s, the norm polynomial
N=P^2+XPQ+Q^2 C has N(x)=P(x)^2. Hence H=N/(X+z) has x as a root iff
P(x)=0. When this holds, both P and Q are divisible by X+x, so N and H
are divisible by (X+x)^2. Conversely that divisibility forces P(x)=0.
This algebraic argument covers x=0, ramified fibers, and nonrational
curve ordinates without choosing or discarding y coordinates.

For fixed h,k,x, P(x) is affine over K in r,s:
P(x)=r*(x^2+z^2)+s*(x+z)+x^3+z^3+(h*x+k)*y(x)+(h*z+k)*yT.
In the Q-root branch h*x+k=0, so there is no unspecified y(x).
Use the literal P evaluation in code, not this eliminated expression.
Intersect P(x)=0 with the complete73-dimensional relaxed family by exact
binary linear elimination. Composition of two independent affine bases
is injective. Therefore all generated (r,s) masks are distinct and exhaust
this branch's relaxed candidates. A family dimension above12 is refused
without enumeration; report its exact rank but leave branch completion false.
For admitted dimensions<=12, test every candidate with literal norm replay,
exact (X+x)^2 division and the multiplicity-safe nilpotence support oracle.
Every all-support H is returned, but source tuple permutations and original
Semaev auxiliary fibers are not reconstructed and no original-root count is set.

The residual cubic H/(X+x)^2 exposes the remaining lower-arity support task.
Keep all other repeated tangencies with Q(x)!=0, squarefree H, pole exceptions
and outer h,k cells open. h=0 is outside this worker's declared four-profile
contract and may not be silently skipped in a general solver.

Tiny independent controls compare literal norms and all GF4 r,s masks,
including144 repeated-pair divisibility checks with h!=0,x!=z. No rational
curve-lift filter is used. Native support is exact R^8=0 modulo monic H
with the source-domain linearized recurrence. Classical norm factors and
linear algebra are prior art; no novelty or qualified timing claim.
90seconds/4GiB/oneCPU; complete serialized output and supervisor lifetime
are charged, and late completion flags/counts are demoted.
