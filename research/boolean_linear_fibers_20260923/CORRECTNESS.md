# Producer correctness argument and validation boundary

This is a bounded mathematical argument plus producer tests, not an external
audit or a machine-checked proof. The executable accepts the declared generated
n12/16/20/24 fixtures. No claim about production index calculus follows.

## Fiber decomposition

Let I be independent in the union graph of nonzero quadratic coefficients.
Every original equation then has the form

    f_e(y,z) = b_e(z) + sum(i in I) a_ei(z)*y_i,

where b_e is quadratic and every a_ei is affine in the outside variables z.
The absence of y_i*y_j terms is checked exactly. Packing equation e into bit e
commutes with GF(2) addition; a zero word means every represented equation is zero.
The worker never drops equations to fit a word.

The independent-set algorithm is exact. For each right-half subset, dynamic
programming compares excluding its first vertex with including it and excluding
its neighbors. Enumerating every independent left subset then queries the optimal
right subset compatible with all its vertices. Each independent set occurs among
these candidates. Cardinality ties choose the smallest binary mask. The Python
reference recomputes selection on every measured fixture.

## Coefficient transport and coverage

The outside constant and linear terms follow the retained exact Gray recurrence.
Flipping a high outside coordinate changes each low-variable column by its fixed
quadratic cross coefficient. Contributions from the four low outside coordinates
are precomputed affine images. The resulting b and every column A_i match direct
restriction of the original polynomial system at every outside point.

The outside order is bijective: reflected Gray order on its high coordinates,
and binary order on up to four low coordinates. Every outside assignment is
therefore covered exactly once on a completed UNSAT scan. Selection and coefficient
setup are charged even when a cap permits no scan. A cap is checked before each
whole block and returns UNKNOWN before any unsupported completion claim.

## Contradiction screens

A bit of b outside the union of all column bits is a row with coefficient zero
and right-hand side one. Such a linear system has no solution.

For d in the fixed sequence 1,2,3,4,5,7,8, define T_d(v)=v XOR (v>>d). This is a
linear map on equation-coordinate words, including the zero-filled upper rows.
If A*y=b, then T_d(A)*y=T_d(b). Applying the same zero-row test after T_d is therefore
a sound rejection rule. The maps need not be injective, and passing these tests
does not prove consistency. Every survivor receives a full solve of the original
linear system.

Scalar and native implementations use identical four-lane groups and identical
early stopping: another redundant-row round is skipped only when all four lanes
are already rejected. NEON compares unsigned words for equality to zero; SSE2 uses
bit-exact equality and lane masks. All loads and stores address complete arrays.
Wider selected sets and partial blocks use the scalar fallback. Tests compare
syndrome values, masks and round counts, including every forced hit lane and bit 31.

## Linear solving and canonical witnesses

The row reference eliminates all coefficient rows, retaining contradictory
constant rows, then recovers a witness with free variables zero. The column method
inserts columns in increasing inside-variable order, recording the column
combination used for each independent basis vector. Reducing b against that basis
decides span membership and recovers a valid coefficient vector y.

This greedy basis produces the smallest binary witness. Suppose another solution
differs at highest position j. If column j was independent when inserted, it cannot
be canceled by earlier columns, contradicting that both vectors solve the system.
Hence the highest differing position is an omitted column, whose greedy value is
zero. The alternative is larger. This also covers zero and dependent columns.
The row and column implementations are checked against full solution enumeration
for all **1,157,359** linear systems through four rows and four columns, counting
all right-hand sides. Additional tests cover full 32-bit equation words and up to
24 columns.

## Recovery, work and evidence

Inside and outside labels form a disjoint partition of original variables.
Scattering a consistent fiber witness through those labels therefore recovers a
solution of the complete original system. Every returned model is independently
checked. Benchmark UNSAT also requires the completed retained search reference;
the streamed trace is diagnostic, not a cryptographic certificate.

An entire screened block is charged even when its first lane succeeds. Rejected
fibers each rule out exactly 2^|I| extensions, but those extensions are not labelled
evaluated full assignments. On completed UNSAT, the total rejected extensions
equals 2^n. Filter rounds count evaluated lane-rounds, including already-rejected
lanes whose group continues. Conditional matrix rank sums are not counts of
independent relations. All setup, failed screens, solves, recovery and validation
remain inside cold totals; missing calibrated-operation and IC costs stay null.
