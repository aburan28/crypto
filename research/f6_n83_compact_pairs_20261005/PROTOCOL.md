# n83 F6 signed-pair compact storage gate

Registered before changing runtime code or running the candidate. This is a
four-summand meet-in-the-middle component of F6-IC, not a complete n83
point-decomposition or discrete-logarithm method.

## Hypothesis and frozen comparison

The retained `F6SignedPairIndex` stores each of 4,108,723 signed pair sums as
a heap-backed `BinaryPoint`, then clones 256 such objects for each query
batch. Store each sum as two packed field words and four `u32` source indices
instead. On the pinned degree-83 ARM64 PMULL path, compute residual x keys
directly from those packed words. Keep the portable path, exceptional group
cases, sign check, stable search order, and independent full-group witness
verification exact. If the base is too large for `u32` indices, construction
must reject it rather than truncate an index.

Baseline is PR #1399 head `a6577e05039851bedd101dd3f38dba11b3ee4829`.
Freeze registered K0 curve `icv1-f2m83-tm6151469093347-debefd74`, its
standard cofactor-projected bases at dimensions 8, 10, and 12 (258, 1,048,
4,054 usable points), and public T001
`(355fb5df7a905f16921eb,5900a390f42d290f1bbe)`. The full base has
8,219,485 unordered pairs and 4,108,723 signed-sum representatives. Use
the same release compiler, target CPU, single Rayon thread, resource limits,
and frozen probe sources for both arms. The reference is the same exact
target query and pair index; no candidate or workload identity is inferred
from this component comparison.

## Correctness, accounting, and decision

Before performance runs, pass the focused geometry tests, small-base
exhaustive four-sum tests, a full-base planted witness replay, and the
ordinary T001 exact-miss replay. Check zero, same-x, opposite-point, and
identity cases. Compare representative counts and output status exactly.
Record `size_of` for the old and new entry layout, code and binary SHA-256,
toolchain, raw exit status/stdout/stderr, build and query wall intervals,
and peak RSS. Run small bases three times per arm. Run the full base in
baseline, candidate, candidate, baseline order after both binaries are
built. Preserve failures and timeouts as rows.

Retain the candidate only if all exactness checks pass, full-base peak RSS
falls at least 20%, full-base query median falls at least 15%, and full-base
build median does not rise more than 10%. If it misses, restore baseline
runtime code, archive the rejected patch and raw results, and still publish
the decision. A query improvement is a target-dependent stage diagnostic;
index construction is target-independent and is reported separately.

This unisolated Apple silicon host cannot promote a CPU wall-time speedup.
The K0 four-summand uniform-target coverage ceiling is `4.662e-12`, so even
a successful component gate does not provide an ordinary relation, complete
F6/F4/F5 comparison, one-target IC interval, or paired rho speedup.
