# Grouped recursion / shared algebra optimization evaluation

Upstream: https://github.com/danadran01/exact-dft-power-saving
Companion: https://github.com/aburan28/crypto-autoresearcher/pull/2063

This is an experiment contract, not a claimed implementation or speedup. The upstream exact complex DFT theorem does not directly apply to finite-field convolution.

## Engineering arms
1. Identify Rust finite-field polynomial multiplication and batched products; implement opt-in grouped request scheduler with deterministic order and baseline fallback.
2. Evaluate expression-DAG hash-consing across polynomial products, sparse polynomial reduction and Macaulay row construction. Record cache memory and hits.
3. Evaluate exact convolution for supported moduli and lengths with root-of-unity validity and coefficient-bound checks; use a correctness-preserving fallback elsewhere.
4. Integrate opt-in kernels into Gröbner F4/F5 and summation polynomial workloads only after exact differential tests.

## Validation
- Baselines: existing Rust algorithms; same inputs, seeds, thread pinning, compiler flags, hardware.
- Sizes: m=31,51,53,83 binary-field toy/experimental curves where applicable; dense/sparse and adversarial polynomials.
- Tests: zero, singleton, non-power-of-two, repeated expressions, high degree, characteristic two, randomized comparisons.
- Metrics: wall time including preprocessing, memory, operation counts, allocations, misses, timeouts, relation verification, independent rank.
- Primary: total seconds per new independently verified relation including failed searches and duplicates.
- Do not promote without >=2x improvement against fastest valid compiled end-to-end baseline and independent validation.
- Preserve immutable receipts, failed trials and exact revision IDs. Stop on correctness mismatch.
