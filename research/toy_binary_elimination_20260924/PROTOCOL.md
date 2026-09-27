# EXP-TOY-BINARY-ELIMINATION-001

Frozen before execution. Dependency: PR #708's constant-coefficient audit.
Class: accounting / bounded algebraic stage diagnostic.

Question: within equal-point-count classes of tiny binary curves, do
polynomial-zero counts, rational lifts, and a fixed Boolean Macaulay filtration
differ? Which comparisons are merely Frobenius conjugacy?

Inputs: GF(8), GF(16), GF(32), moduli 0b1011, 0b10011, 0b100101;
curves y^2+xy=x^3+c with every c != 0; subspace basis [1,2]; its squared
basis as a transported control. Every ordered triple in V^3 is enumerated.
No seeds, target scalar, DLP algorithm, larger-field mode, or cloud workers.

Curve grouping: enumerate every affine point and add the identity. Over a
finite field equal point counts certify isogeny of elliptic curves (Tate).
Record Frobenius orbits of c separately. This supplies isogeny-class
membership, NOT explicit neighbor maps or volcano levels. Do not label these
pairs descending neighbors.

Rational lifting: for every polynomial-zero triple, enumerate rational y
coordinates and verify the group sum. Count x-triples admitting a lift, not
sign multiplicities. Repeated points and zero x are retained and labeled.

Elimination: six Boolean input variables. Work in the 64-dimensional Boolean
algebra F2[b]/(b_i^2+b_i). For each nonzero coordinate generator g, insert
all squarefree monomial multiples m*g when deg(m)+deg(g) <= D, for D=0..9.
Reduce Boolean products before bitset Gaussian elimination. Use fixed numeric
monomial order and record matrix rank, generated rows, and row XORs at each
D. Stop only after the fixed D=9 cap. Report the first D where rank equals
64 minus the exhaustive zero count as the saturation degree of THIS
filtration. It is not F4/F5 solving degree or degree of regularity.

Correctness gates: validate the group law on all point pairs; polynomial
zeros containing a rational point in each x-fiber must admit some rational
zero-sum lift; Frobenius transport preserves point and lifted-triple counts;
each generated row vanishes on all enumerated polynomial zeros; full rank at
D=9 equals 64 minus zero count. Fail immediately on a mismatch. Preserve raw
results and all stages; do not discard unsatisfiable cases.

Reference is exhaustive enumeration, and the full Boolean ideal dimension
identity is the exact boundary. A successful audit has zero mismatches.
Different XOR counts alone are not evidence of a cryptanalytic improvement.
Separate within-class comparisons by solution count and Frobenius orbit.
End-to-end operations, S, rho ratio, speedup, and regularity remain null.
The frozen WDSat suite is inapplicable to this correctness-only toy diagnostic.

Cost scope: field setup, point enumeration, lift enumeration, and polynomial
construction are outside the reported matrix row/XOR counters. These counters
describe only elimination, not total cost or a runtime comparison.
