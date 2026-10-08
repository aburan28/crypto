# K_0 n=83 pair-index capacity experiment

Registered before baseline on 2026-10-04. Compare the one-representative
pair-index implementation with and without preallocating its hash lookup to
the exact maximum unordered-pair count. The query algorithm and point lists
remain identical. Retain the preallocation only if exact correctness agrees,
and exploratory median build time is no worse than 1.10 times the baseline
at either size. Query times and memory tradeoffs must be reported.

Run `f6_wide_n83_index_probe` in release mode on the same local host and
compiler for the pinned K_0 curve, actual standard-subspace cofactor-projected
bases of dimensions 8 and 10, and public T001. Warm up once per size; record
three fresh index builds and corresponding queries per arm. Preserve raw
durations, no-witness/witness status, source hashes, pair counts, and any
failures. The baseline is the retained one-representative index before this
capacity change.

This host has no isolation receipt. The result is an exploratory stage
diagnostic only, not a controlled CPU speedup, natural-yield measurement,
complete F6/F4/F5 comparison, or one-target index-calculus result.
