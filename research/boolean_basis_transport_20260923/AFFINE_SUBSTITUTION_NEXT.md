# Prospective contract: eliminate general affine relations

Status: design only. The current performance arms transport row spaces and use
unit affine consequences, or solve an entirely affine residual system. They do
not substitute a general affine relation while nonlinear equations remain.
`affine_probe` counts such opportunities on discovery traces; its counts are not
a speedup or a count of independent relations across the search tree.

Linear-variable elimination is established prior art, including the ElimLin
steps described in [the primary polynomial-system treatment](https://arxiv.org/abs/2304.07820)
and [this explicit ElimLin description](https://eprint.iacr.org/2014/893.pdf).
The proposed experiment concerns its complete cost on this frozen generic
Boolean workload, not a novelty claim or a cipher-specific application.

## Coefficient and recovery identities

Let a consistent affine tail have rank r in k remaining variables. RREF gives an
injective parametrization x=A y+b of its solution space, with k-r free variables.
For

\[
f(x)=c+\sum_i\ell_i x_i+\sum_{i<j}q_{ij}x_ix_j,
\]

the transformed Boolean polynomial has coefficients

\[
c'=c+\sum_i\ell_i b_i+\sum_{i<j}q_{ij}b_i b_j,
\]

\[
\ell'_s=\sum_i\ell_i A_{is}
+\sum_{i<j}q_{ij}(A_{is}A_{js}+b_iA_{js}+b_jA_{is}),
\]

\[
q'_{st}=\sum_{i<j}q_{ij}(A_{is}A_{jt}+A_{it}A_{js}),\qquad s<t.
\]

All arithmetic is in GF(2). The diagonal contribution A_is A_js belongs in the
linear coefficient because y_s squared equals y_s in the Boolean quotient.
Omitting it silently changes the equations. The degree remains at most two,
although support may expand and new quadratic coefficients may cancel.

Recovery maps compose as

\[
(A_0,b_0)\circ(A_1,b_1)=(A_0A_1,\ A_0b_1+b_0).
\]

The solver must reconstruct original variable values from its final free-variable
assignment and verify the original equations. This bijection does not require a
unique nonlinear solution. An inconsistent affine tail is UNSAT; a resource cap
is UNKNOWN and cannot stand in for that contradiction.

## Required experiment

Implement a direct polynomial substitution reference and an exact packed
coefficient implementation. Charge construction of A and b, any monomial-image
schedules, cancellations, density growth, row reduction, branching, recovery and
validation. If original-equation images are retained for branch selection, charge
their transformation too. The previous immutable-quadratic-coefficient shortcut
does not apply after general affine coordinate changes.

Correctness should include all quadratic polynomials in three variables under
every affine map, checked against an independent truth-table/Mobius-transform
oracle; repeated map composition; zero-rank and inconsistent tails; Boolean
diagonal terms; and exact full-solve results against the current controls. Include
recovery maps in the new policy's trace so equality checks cover coordinate changes.

Freeze the new cost gate before performance measurements. Use the fastest
retained complete method, including packed_state, as the reference; retain current
failed cases and add unused holdouts. A faster substitution kernel or fewer free
variables alone does not establish a complete-solve gain. Full index-calculus and
rho comparisons remain outside this proposed generic experiment.
