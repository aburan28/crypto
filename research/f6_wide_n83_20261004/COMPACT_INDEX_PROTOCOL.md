# K_0 n=83 compact pair-sum representation

Registered before fresh baseline runs on 2026-10-04. The current index keeps
each pair sum as a full `BinaryPoint` with two inlined multiword field
elements. The candidate stores each sum as two packed `u128` coordinates
plus compact point indexes, and reconstructs field elements only when a
target query scans a chunk. The lookup remains exact and retains one
representative pair per group sum. No query, factor-base, or witness rule
changes.

Run the existing release index probe on the pinned K_0 dimension-8 and
dimension-10 bases and public T001, with its one warmup and three build/query
observations per size. Run the full dimension-12 index probe once to capture
all 8,219,485 pairs, one exact query, and peak process RSS via in-process
`getrusage`. Record these baseline rows before changing code. Then build and
run the same probes after the compact representation change on the same
host and compiler. Preserve raw rows, exit status, source hashes, and
failures. Exact witness/no-witness outcomes and pair counts must match.

Keep the compact implementation only if correctness matches, dimension-8
and dimension-10 exploratory median build times each stay within 1.10 times
their baseline, full-base query time stays within 1.10 times its baseline,
and full-base peak RSS is at most 0.70 times its baseline. A failed gate
causes a code revert while preserving both receipts. All timing ratios are
unisolated component diagnostics, not controlled CPU speedups or complete
F6/F4/F5 or index-calculus comparisons.
