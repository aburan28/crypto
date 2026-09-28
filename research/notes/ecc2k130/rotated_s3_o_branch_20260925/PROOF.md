# Branch-complete rational S3-chain and one-hot CNF semantics

This extends the finite-fibre theorem in [PR #777](https://github.com/aburan28/crypto/pull/777).
Let `E/K: y²+xy=x³+1`, `L(a)` be the complete affine rational fibre
at x=a, and let every factor slot contain entire fibres. Let a state be
O or a finite x. A finite state denotes *both* rational signs in `L(x)`;
x=0 denotes the singleton point `(0,1)`. A path records the state after
each addition, including the terminal state.

For an O predecessor, `O+L(a)=L(a)`, so the only next state is a.
For a finite predecessor u and factor a, [#777's local theorem](../rotated_s3_fibre_theorem_20260925/PROOF.md)
identifies the finite next states with exactly the roots of
`S3(u,a,X)`. An O next state occurs exactly when the chosen points are
inverses, which is possible iff u=a because equal x is the inverse pair
for this binary curve. This includes u=a=0: `(0,1)` is its own inverse,
`S3(0,0,X)=1` has no finite root, and the sole outcome is O. For equal
nonzero x, equal signs give the finite doubling root and opposite signs
give O. For distinct x, O is impossible; one zero input gives the unique
finite double root and two nonzero inputs give the two distinct finite
roots. These cases exhaust every transition without division by a
possibly zero coordinate in the exporter.

Inductively, a valid state path is realized by a signed factor tuple.
Suppose a prefix realizes point U at its recorded state. For a finite
next transition, the local theorem supplies some representative U' of
the same fibre and a factor P with the requested next state; for an O
transition, the inverse condition supplies such U',P. If U'=-U, negate
**all earlier factors**. Their previous finite x states are unchanged and
O remains O, and complete factor fibres retain every negated factor.
Then append P. When the old state is O, U'=U=O and either allowed factor
sign can be appended. Thus every local transition can be made consistent
with the entire earlier path. Conversely, every signed point tuple
projects to one allowed state path by the local cases.

For a finite terminal state x, globally negating all factors retains the
factor-x tuple and every prefix state while swapping the final point
between the two members of `L(x)`; when x=0 the final point is its own
negative. Thus a terminal-state literal accepts a specified exact
rational target R iff there is a signed tuple summing to R. Terminal O
accepts the identity and no finite point. The claim requires rational
complete fibres; formal nonliftable x or a single permitted sign breaks
this correspondence.

In the exported DIMACS, every factor and state group is exactly one hot.
For each transition and each triple of its group values that is *not* in
the local table, one clause excludes precisely that triple. Therefore a
CNF assignment chooses exactly one factor-x tuple and one state path,
and it satisfies all clauses iff every local transition is valid. Adding
the final-state unit literal for R enforces exactly its terminal state.
Combining this elementary CNF equivalence with the induction proves the
branch-complete semantics. The independent verifier checks the actual
clause list with a separately implemented full-point group law, so this
proof is not used as a substitute for exporter validation.

No claim is made about a direct S6/S7 resultant, an ANF/CNF obtained by
clearing denominators, a sign-incomplete factor base, SAT solver behavior,
relation rank, n131 cost or the public discrete logarithm.
