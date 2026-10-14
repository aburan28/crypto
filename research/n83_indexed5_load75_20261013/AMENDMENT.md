# Dense indexed-root capacity amendment

The first version of this branch requested `ceil(2*S/3)` entries for `S`
retained states. That assumed about one distinct canonical root key per
state. The completed K2543 symmetric-S3 cloud receipt instead counted
536,737,856 keys from 268,375,505 states, almost two keys per state.
`ceil(2*S/3)` can therefore leave the indexed table without a free slot.

The initial local focused test, preserved in `first_witness.log`, exhausted
its indexed table after this request was added to the six-representative
fixture. The corrected test in `first_witness_backtrace.log` and the full
16-test example suite in `example_tests.log` passed with `ceil(4*S/3)`.
The constructor allocates the next power of two above twice the request;
at K2830 the checked source-loop bound enters a 2^30-slot, five-GiB table.

The large K2543 density comparison checks equal base hashes, root-key
counts, recovered scalars, and independently valid four-point relations.
Different rank-worker completion order can change the subsequent target
query sequence, so target-row identity is retained as a diagnostic rather
than a correctness condition. First-witness equality for the full large
index remains a source-level proof obligation; the N83 source test checks
it for the six-representative fixture at both table densities.
