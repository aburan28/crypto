# Fused verified packing and trusted-u64 column hashing

On exact layout hits, the solver now forms, cancels, verifies, and packs each
shifted polynomial product directly into its final flat row. It accepts the
cached layout only when every product monomial maps and every cached column is
observed; otherwise it discards the partial matrix and performs full support
discovery. `KIC_F4_DISABLE_FUSED_PACK=1` restores materialized support rows.

The cached layout's monomial index now uses a deterministic SplitMix-style
hasher for trusted internal `u64` masks. Hash collisions remain exact-map
collisions and are resolved by key equality. `KIC_F4_STD_COLUMN_HASH=1`
restores standard hashing.

## Same-binary paired result

| Measurement | Merged-policy control | Selected | Ratio |
|---|---:|---:|---:|
| Complete process wall | 110.673 s | 100.891 s | 1.097x |
| Child user CPU | 110.256 s | 100.123 s | 1.101x |
| x-chained found median | 11.171 s | 8.022 s | 1.393x |
| Symmetrised refuted median | 6.895 s | 4.874 s | 1.415x |
| Peak RSS | 464,715,776 B | 464,338,944 B | 0.999x selected/control |

All verdicts, effort, built degree, FFD, gates, and inconclusive counts are
identical.

The sym profile improves from 6.268 s to 4.810 s. Construction falls from
3.344 s to 1.903 s while reduction, calls, rows, columns, and 2,415,432,155 word
XORs remain equal. Disabling fusion raises the profile to 5.133 s; retaining
fusion with standard hashing raises it to 5.875 s. The x profile improves from
0.774 s to 0.597 s with equal work counts.

All frozen and holdout ladder digests and word-operation counts match. Every
row improves, from 1.07x to 1.37x.

Relative to the retained pre-linear-tail legacy, cumulative process wall is
more than 2.16x faster; x and sym arm medians improve by about 6.78x and 6.36x.

Rejected directions remain recorded: exact per-polynomial row caching produced
253,921 misses and zero hits, and paired child specialisation changed the
three-process median by only 0.17%.

Classification: `N31_BOOLEAN_F4_FUSED_HASH_PAIRED_T2_PASS`.

This is bounded public-synthetic decomposition-construction evidence. It does
not measure relation yield, a complete index-calculus run, asymptotic
complexity, or a crossover against Pollard rho, and it is not a key-recovery
result.
