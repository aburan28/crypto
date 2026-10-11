# K_0 n=83 pair-index representation experiment

Registered before the first baseline run on 2026-10-04. Compare the
current pair index at commit `7ea5f6470` with a change that stores one
representative pair per distinct group sum instead of a heap vector of
pairs. No four-summand search semantics change: repetitions are allowed,
and any pair with the same sum can complete the same target residual.

Use the pinned K_0 curve, standard polynomial-subspace bases of dimensions
8 and 10, and their public cofactor-projected point lists, in deterministic
order. Build the full unordered pair index with its cap equal to the exact
pair count. Use public T001 as the query and verify its on-curve and
subgroup properties. For each dimension and each code version, warm up one
build and one query, then measure three fresh index builds and one query on
each index, preserving exact witness/no-witness status and every duration.
Build both versions in release mode with the same compiler on the same host.
Run the baseline first, then the candidate. No input substitution after an
outcome. Retain source hashes, all raw rows and failures.

The comparison is an unisolated arithmetic/index diagnostic, not a
controlled CPU speedup or a complete IC candidate comparison. It cannot
establish natural relation yield, a target logarithm, or an F6/F4/F5 or
IC/rho speedup. Keep the implementation change only if all correctness
results match and index-build time does not regress on either size by more
than 10% in exploratory medians.
