# Stage 178 protocol: reused dense reducer index in current F4

## Hypothesis

Stage 177 proved that adaptive submask enumeration reduces symbolic reducer
lookups by 97.5%, but a freshly allocated hash map made paired time worse. The
frozen Phase B fixed-X1 systems have only 18 Boolean variables. A direct array
indexed by the exact monomial mask can avoid hashing, use one cache-friendly
load per submask, and be reused across every F4 step in one solve.

This follow-up keeps the current repository F4 algorithm and every Stage 177
tie and correctness rule. It replaces the rejected hash map with a reusable
dense array only for domains of at most 20 variables; larger domains retain the
current linear reference. The rejected hash implementation remains only as the
Stage 177 patch artifact.

## Candidate and control

- Candidate default: allocate one `u32` entry per possible monomial mask once
  per F4 call, reset only touched entries between steps, and adaptively
  enumerate submasks when `2^deg(m)-1 <= active.len()`.
- Same-binary control: `F4_F2_INDEXED_REDUCERS=0`, the current deterministic
  linear scan.
- Export dense-index peak bytes, submask probes, linear tests, stable
  row-equivalent XORs, actual table-assisted XORs, full process CPU/wall, and
  RSS. Dense index memory is included in the solver peak-memory upper bound.

## Correctness gates

1. Exhaustive differential lookup tests compare the dense index and forced
   linear reference over every monomial in randomized small domains, including
   duplicate leading monomials and both adaptive paths.
2. All current Boolean-F4 certified-basis, Buchberger-agreement, budget, and
   `BlockTables` equivalence tests pass.
3. The fixed-X1 specialization proof and Phase B backend tests pass in both
   candidate and control modes.
4. On the frozen target both modes return exhaustive UNSAT, authenticate the
   same source, visit 512 masks, complete 242 systems, and reproduce equation
   fingerprint
   `02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb`.
5. Pair, matrix, basis, extraction, row-equivalent XOR, and actually-performed
   XOR counters agree exactly. Only divisor-index and timing/memory counters
   may differ.

## Frozen benchmark and decision

Use the same opened `n=59, ell=9, m=3` target, twelve Rayon workers, X1 batch
512, 300-second internal budget, and 360-second watchdog. Run three interleaved
pairs in fixed order:

```text
linear, dense, dense, linear, linear, dense
```

Select the dense index only if every correctness gate passes and the median
paired dense/linear wall and total-core ratios are both below `0.97`. Report
RSS and absolute ratios to Stage 174 and direct MITM separately. Charge the
fresh exact-commit build and every process. Single-core time remains null.

This remains one-target solver engineering. It cannot establish relation yield,
unknown-scalar recovery, a rho crossover, independent reproduction, novelty,
or Koblitz index-calculus SOTA.
