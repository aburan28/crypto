# Chained S3 at the larger sizes

This control completes the solver_07/08 comparison on the exact same targets,
first/enumeration modes, x-subspace bases, three-summand restrictions, and
five-second all-phase budgets. The implementation reuses the previously
validated chained S3 frontend and charges model rejection and exact lifting.

Before the campaign, complete projected-set equality is checked on all 43
five-bit affine targets using the independent pair oracle. Larger-field
SAT models are checked against CNF, signed curve sums and the oracle; full
enumeration must match its complete projected set. Timeouts remain unknown.

See ../scaling_23_29.md for the combined outcome. This is a sequential control
run after the image-space campaign, not a claim of repeated-run statistical
significance. Comparable SAT operation accounting remains absent. Natural
yield and full ECDLP cost cannot be inferred from the supported-target stratum.
