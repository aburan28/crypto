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

The T025 control passed: all unchanged report fields match the reference,
including its first relation, and checked Sage independently replayed scalar
`429728658711`. The target made 3,280 charged S3 calls versus 3,272 before
look-ahead prefetch. Its exploratory online interval was 0.949 ms versus
1.701 ms for the single earlier reference run. These T025 runs select the
candidate; they do not enter the fresh holdout ratio.

`SCREEN_FREEZE.json` binds previously unused N53 points T189–T196, both
candidate manifests, source and release binary hashes, 14 preparation
workers, one target worker, 32 GiB peak RSS, and alternating AB/BA order.
Every arm receives one point per fresh process with full reusable preparation
outside the online clock. The prespecified advancement gate is all 16 scalar
recoveries verified and a held-out geometric reference/combined online ratio
of at least **2×**. The host-isolation and 10% setup-noise gates still apply
before a controlled wall-time claim.
