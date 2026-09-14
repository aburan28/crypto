# Next experiment: reject coefficient blocks, not just individual functions

The branch lower bound identifies a control-flow problem. Both frozen
hybrids pay for every conditioned-root branch on complete enumeration; an
early check in recovery cannot remove that work. The next candidate must
reject a set of (h2,z) choices before evaluating each member.

## A useful exact equation, including subfield coefficients

The existing `pullback_incremental` experiment normalizes a restricted linear
equation on the a2=0,a6=1 curve. For the general curve

    y^2 + x*y = x^3 + A*x^2 + B,

take a nonzero target abscissa r, target ordinate s, and set

    h2 = b^2+b+r,    u=h2+z,    t=r+z,
    gamma=t/(b*r),   eta=z*u/b, D=(r+s)/r.

Here b,z,t,u are nonzero, h2,z lie in V, and z differs from r. Let
I_u={w^2+u*w:w in V}. The exact condition on the residual pair product v is

    v in I_u,
    gamma^2*v^2 + (z/r)*v = eta^2 + b^2*B/r^2.             (1)

Recover a=gamma*v+r+b*D+eta. The cancellation uses the target curve equation:

    D^2+D+r+A = B/r^2.

Thus the right-hand side need not retain D or A explicitly. For the
reconstructed a, if E is left side plus right side of (1), direct substitution
in the general norm gives

    (H_1 + z*u) + v = E,       H(z) = t*E.

This proves equivalence to the support-product and conditioned-root
equations. The original nonzero-root, distinctness, target exclusion and
signed-curve verification still apply. The r=0 chart remains separate.
The test checks these residual identities directly on all F4 coefficient
choices at GF64 and fresh larger-field samples; it does not time a solver.

## Concrete proposal and rejection criteria

1. Express v=w^2+u*w with w in V and retain equation (1) as an F2-linear
   equation in w for fixed h2,z,b. This fixed-branch linearization already
   exists in the pullback approach; it is not the new scaling claim.
2. Partition h2 and z coordinate bits into blocks. Build a sound linear
   relaxation of (1) for the unfixed bits, treating nonlinear monomials as
   independent variables where necessary. An inconsistent augmented matrix
   certifies that the entire block has no solution. Consistency is only a
   necessary condition; it never certifies a relation.
3. Expand only surviving blocks; verify every recovered relation. Record
   the number of branches eliminated, preprocessing/rank work, failed blocks,
   residual solving, extraction and exact checks. The constraint b^2+b=h2+r
   must be retained; b values cannot be independently discarded.
4. Compare a reference with no block pruning, the block candidate, and the
   direct S3 table on identical cold and named-batch targets. Keep both
   first-hit and complete enumeration modes, old inputs and fresh holdouts.

No block relaxation or pruning gain has been implemented in this round.
This proposal is falsified as an engineering candidate if its bound/rank
work costs as much as the eliminated work, or if most blocks survive until
singletons. Any false rejection invalidates it. A new pre-measurement
contract and matched benchmarks are required before testing a speedup.

The separate trace-fiber count is useful as an exact boundary calculation.
If every fiber has positive mass, it supplies no whole-target rejection;
it must not be promoted into a claimed shortcut. A richer trace constraint
inside the coefficient equation remains optional and unproved.
