# Hoisted exact tangency search

The v1 branch coverage proof is unchanged: fixed four h,k profiles, source
x!=0,z with Q(x)!=0, exact Artin-Schreier fibers, unique tangent r,s, full
product-span and literal-norm/support checks. No original tuple, squarefree
or outer-parameter completion is inferred from this component.

x^2,x^3,(x+z)^-2 are identical for both ordinates and all h,k arms. Compute
them once for each x with a rational curve fiber. z^3 and (h*z+k)*yT are
invariant across the entire scan. hx+k is shared by both ordinates.
The cached helper uses exactly
s=(h+1)*x^2+k*(x+w),
r=[s*(x+z)+x^3+z^3+(hx+k)y+(hz+k)yT]/(x+z)^2.
The same new helper is checked against the independent exhaustive norm/
derivative oracle on15664 small parameter states, not just against v1.
All repeated Q-root and x=0 cases stay in their separate proven branch.

Search denominator counts are exact scoped primitive counts: v1 computes
one inverse of x for each nonzero visited source plus one inverse per
admitted tangent arm; v2 computes one inverse of x and at most one inverse
of (x+z)^2 per source with rational fibers. Count field/input verification
and certificate costs separately; cold process time includes them all.
No time ratio is formed against v1's incomplete20-second execution.

New scientific prediction: all2^18 source coordinates for four profiles
complete within30 diagnostic seconds after hoisting. This worker has a
60-second native search cap and75-second supervisor; increasing the budget
is not itself claimed as an improvement. The concrete changed assumption
is shared denominator/invariant reuse. If partial, null exact counts and
false branch flags remain; redesign arithmetic instead of an identical retry.
4GiB/oneCPU. Preserve source/input hashes and independent small oracles.
This is native implementation engineering/classical algebra, no novelty or
qualified8x/F4/F5/F6 claim. Final serialization and supervisor charge remain.
