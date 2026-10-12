# Conditioned quadratic support on fresh targets

This frozen campaign implemented the linear-image proof from solver_03.
For a fixed nonzero u in V, Q=X²+uX+v splits into two distinct roots in V
exactly when v belongs to the image of w -> w²+uw on V. Gaussian elimination
constructs that image, of dimension d-1. Membership becomes binary linear
constraints on v. The circuit derives z=h2+u and retains H(z)=0, z in V,
H(0)!=0, H(r)!=0 and Q(z)!=0. These conditions are equivalent to the
restricted cubic support condition; see solver_03 for the proof.

The code builds a separate specialized circuit for every u branch. Branch
construction, solver loading, searching, extraction and verification all
share a three-second soft budget; remaining time is divided among remaining
branches. First-relation search and complete projected enumeration are
separate runs. Unknown branches cannot establish global UNSAT.

All 144 trials completed without validation errors. The 6,944 coefficient
and branch checks matched scalar output vectors and existential cubic
support. The chosen validation target had zero admissible functions, so
that exhaustive check alone supplies no positive acceptance coverage.
Positive models in the target campaign are independently checked against
CNF, recovered curve witnesses and the exhaustive group oracle.

Four uniform and four supported targets per size exclude every target from
solver_02/03. Strata may overlap and are reported separately. The candidate
resolved all nine-bit first-relation cases, but only one of eight eleven-bit
first-relation slots and none of its full enumerations. This implementation
fails the reliable-eleven-bit objective. It is retained as negative evidence.

The combined result tables, strongest measured controls and remaining goal
are in ../goal_round_20260913.md. All raw rows, contract, source hashes and
publication mapping are preserved here. Timing is preliminary; comparable
SAT operation accounting is absent. This is an engineering experiment and
changes no counting boundary.
