# Contiguous solver matrix for exact linear-tail F4

The selected solver previously allocated one `Vec<u64>` per Macaulay row. This
iteration packs the solver-only matrix into one row-major buffer and performs
forward high-column elimination directly on contiguous row slices. Public
full-RREF F4 matrices remain unchanged. The flat path is selected at 24 or more
variables and `KIC_F4_FLAT_MATRIX=0` restores the nested control.

## Final same-binary pair

| Measurement | Nested control | Flat selected | Ratio |
|---|---:|---:|---:|
| Complete process wall | 113.528 s | 105.112 s | 1.080x |
| Child user CPU | 112.958 s | 104.558 s | 1.080x |
| x-chained found median | 11.930 s | 10.109 s | 1.180x |
| Symmetrised refuted median | 7.454 s | 6.105 s | 1.221x |
| Peak RSS | 482,017,280 B | 480,067,584 B | 0.996x flat/nested |

Both policies use the same executable and public targets. Verdicts, effort,
built degree, FFD, gates, and inconclusive counts are identical.

The sym solver profile improves from 7.256 s to 6.217 s. Its reduction phase
improves from 3.763 s to 2.795 s while rows, columns, F4 calls, and 2,415,432,155
word XORs remain exactly equal. The x profile improves from 0.875 s to 0.771 s,
with reduction improving from 0.440 s to 0.330 s and identical work counts.

The scoped final subprofile prices flat construction at 1.437 s for shifted-row
formation and 1.784 s for verified dense packing. Reduction costs 2.765 s.

All frozen and holdout verdict digests and word-operation counts match. The
above-gate rows improve about 3–6%; below-gate rows are unchanged apart from
timing noise.

Relative to the retained pre-linear-tail legacy, cumulative whole-process wall
is more than twice as fast, and the x and sym medians improve by more than five
times.

Rejected controls:

- Exact per-polynomial row-product caching records 253,921 misses and zero hits,
  regressing wall from 6.93 s to 8.25 s.
- Paired child specialisation changes the three-process median by only 0.17%;
  it shares a support scan but not shifted-row construction.

Classification: `N31_BOOLEAN_F4_FLAT_MATRIX_PAIRED_T2_PASS`.

This is bounded public-synthetic decomposition-linear-algebra evidence. It does
not measure relation yield, a complete index-calculus run, asymptotic
complexity, or a crossover against Pollard rho, and it is not a key-recovery
result.
