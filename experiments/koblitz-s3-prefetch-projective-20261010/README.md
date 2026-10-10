# Combined S3 target prefetch and projective scalar checks

This candidate retains the full-rank N53 C1200 factor base, relation
collector, root table, and one-worker target search. It prepares up to 32
canonical partner keys per eight-state batch, prefetches their exact root
table home slots, and then performs the original witness checks in order.
Target subgroup validation and final scalar replay use the portable
López–Dahab projective multiplier. Both changes are charged within the
single-target online interval.

The source is derived from the verified late-prefetch kernel and the
verified projective scalar-check kernel. Its 16 release tests include exact
projective/affine agreement on N41 and N53 subgroup-boundary scalars. The
T025 control compares the same public point with the frozen C1200 affine,
non-prefetch reference. It requires the same first four-point relation,
factor-base digest, full rank, recovered scalar, and all unchanged report
fields. Prefetch look-ahead may increase the final batch's S3 call count;
the complete charged count remains in the report. A fresh N53 holdout screen
will follow only if the correctness and 32 GiB memory gates pass.
