# ECBench n−1 / large-prime correctness protocol

Status: frozen implementation and correctness protocol. This document does not make a performance claim.

## Hypothesis

For an explicit binary curve and target imported from the ECBench corpus, exact relations

`[a]P + [b]Q = F_1 + ... + F_m`

with `m = n − 1` when requested can be collected over a nested factor base. After projection into the prime-order subgroup, zero-, one-, and two-large-prime partial relations can be eliminated exactly and the remaining modular system can recover the target logarithm. The answer must agree with independent scalar multiplication; the planted logarithm is never solver input.

## Frozen inputs and references

- Source baseline: `db8703203234949ce650b62d27a7a4b698fa8ee5`.
- External corpus: `/workspace/ecbench/corpus/manifest.json`, SHA-256 `ff1174299d4d89c88501482b2cc61dd5ca783591f703f90aaf95609b96a9d362`.
- Corpus families: all eight records in that manifest. The smallest family is the required end-to-end target; the other seven must at least import and validate exactly.
- Generic reference: `rho.signed_frobenius` on the same explicit Koblitz curve, subgroup, generator, and target through native `ecbench`.
- Independent verifier: the general binary-curve implementation, separate from the fast relation-search implementation.

## Success conditions

1. Exact manifest import records the raw file's SHA-256 and byte length, rejects duplicate family ids, and checks the field polynomial, curve membership, `#E = h r`, `[r]P = O`, and the planted public equality without giving the planted scalar to the solver.
2. Bounded `m = n − 1` decomposition runs on a tractable small fixture and refuses, with a state-count explanation, when its configured state cap is exceeded.
3. Zero-, one-, and two-large-prime modes recover the same scalar on the smallest imported family; the double-large-prime mode records and eliminates genuine two-large-prime rows.
4. The recovered scalar passes both fast and independent point verification.
5. Native `ecbench` can run `ic.large_prime` and a generic reference on one exact `BinaryExplicit` workload.
6. The installed `ic` reports build-time source provenance even when invoked from another repository.

## Stop conditions

- Relation collection stops at `max_trials`.
- Meet-in-the-middle table construction stops before `max_states` is exceeded.
- No family is silently downscaled and no `n − 1` request is substituted with a smaller summand count.
- Larger-family runtime is reported as unestablished when the frozen resource caps are insufficient.

## Accounting boundary

Every elliptic-curve addition and multiplication used by the search is charged through the repository's operation counters. Hash-table probes, combination enumeration, canonicalization, and modular row operations are reported as named unpriced native counters. Wall time is diagnostic only. No speedup, asymptotic-transfer, or ECC2K-130 claim follows from this protocol.
