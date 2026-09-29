# Measured results: complete toy index calculus

All **360/360** full-DLP runs verified. The post-run audit checked **3,306** decomposition calls across **172** distinct intermediate targets, including **498** complete empty answers. All **215** exhaustive five-bit adapter/target checks passed. **18** frozen historical scalar and operation records reproduced exactly.

Classification: **accounting/integration**, not a cryptanalytic advance. Cached Nagao was the fastest tested IC implementation; plain rho was faster overall. No common-unit operation speedup, ECC2K-130 throughput, or runtime extrapolation is established.

Host: Intel Xeon Platinum 8370C at 2.80 GHz, Python 3.12, CryptoMiniSat/pycryptosat 5.14.7, one solver thread, sequential runs. Ten seeds × three repetitions per fixture and variant. Five-bit runs use seven distinct final targets; nine-bit runs use ten. Fresh seeds can repeat targets in a tiny group. Process peak RSS was 34,752 KiB; this is a campaign high-water mark, not per-variant memory.

## Cold elapsed-time diagnostics

Per-row medians below include all full-DLP phases, including exhaustive toy group setup. Each cell has 30 verified runs and zero incomplete runs. Field bits are not subgroup bit lengths: the subgroup orders are only 11 and 127. The ratio is a geometric mean of paired per-seed median times; it need not equal the ratio of the two displayed global medians. Lower is better.

| Variant | GF(2^5), N=11: median ms | GF(2^9), N=127: median ms | Paired nine-bit time / rho (95% interval) | Full-DLP S | Calibrated cost/rho | Calibrated cost/floor | Class |
|---|---:|---:|---|---|---|---|---|
| Nagao: quadratic image | 39.46 | 409.70 | 4.22 (3.49–5.12) | null | null | null | Integration/accounting |
| Nagao: coefficient pullback | 28.46 | 299.75 | 2.89 (2.51–3.34) | null | null | null | Integration/accounting |
| Nagao: cached pullback | 22.89 | 222.53 | 2.33 (2.09–2.61) | null | null | null | Integration/accounting |
| Symmetric S4 + CryptoMiniSat | 49.75 | 1434.49 | 14.52 (10.65–19.48) | null | null | null | Integration/accounting |
| Chained S3 + CryptoMiniSat | 49.93 | 1435.44 | 14.67 (10.80–19.62) | null | null | null | Integration/accounting |
| Plain Pollard rho | 4.72 | 91.50 | 1.00 (1.00–1.00) | null | null | null | Reference |

These wall-time intervals describe these implementations on these ten seeds. They do not establish a calibrated algorithmic speedup or generalize to real key sizes. The toy rho reference omits Frobenius/negation acceleration. No F4 or WDSat result was measured in this panel.

## Where the nine-bit time went

Mean milliseconds per complete run, derived from the summed exclusive phase counters. Relation collection includes generation, per-call setup, solving and independent relation verification. Individual descent includes its decomposition calls. Adapter subtimers in raw records are nested and are not added again.

| Variant | Setup | Factor base | Relation collection / rho walk | Matrix | Individual descent | Final check + orchestration |
|---|---:|---:|---:|---:|---:|---:|
| Nagao: quadratic image | 84.366 | 4.623 | 293.051 | 0.280 | 21.977 | 0.830 |
| Nagao: coefficient pullback | 83.934 | 4.791 | 169.620 | 0.253 | 13.509 | 0.853 |
| Nagao: cached pullback | 81.091 | 4.889 | 119.918 | 0.263 | 9.642 | 0.892 |
| Symmetric S4 + CryptoMiniSat | 80.905 | 4.604 | 1338.587 | 0.375 | 68.210 | 0.910 |
| Chained S3 + CryptoMiniSat | 82.232 | 4.911 | 1331.827 | 0.394 | 70.158 | 0.856 |
| Plain Pollard rho | 81.714 | 0.000 | 10.725 | 0.000 | 0.000 | 0.803 |

For cached Nagao, relation collection accounts for 55.3% of nine-bit cold time. Matrix elimination is tiny here; that does not predict large-matrix cost. Exhaustive curve setup dominates the tiny rho runs and is unsuitable for extrapolation. Process module import took 0.018396 seconds outside per-fixture cold times.

## Operation vectors, kept in their own units

Each row sums 30 complete cold runs. Binary-field inversions expand into counted field primitives. Scalar arithmetic is separately modulo N. SAT rows report only the field work actually counted; symbolic encoding and SAT search are **unpriced**, so their smaller field vector cannot be interpreted as a smaller total. No uncalibrated component sum is used as a performance score.

| Field bits | Variant | Field additions | Field multiplies | Field squares | Scalar add / multiply / inversion | Decomposition calls |
|---:|---|---:|---:|---:|---|---:|
| 5 | Nagao: quadratic image | 280,335 | 282,630 | 98,706 | 1,533 / 1,620 / 150 | 222 |
| 5 | Nagao: coefficient pullback | 165,165 | 229,947 | 47,097 | 1,533 / 1,620 / 150 | 222 |
| 5 | Nagao: cached pullback | 122,169 | 169,437 | 44,322 | 1,533 / 1,620 / 150 | 222 |
| 5 | Symmetric S4 + CryptoMiniSat (partial) | 63,567 | 115,500 | 18,702 | 1,824 / 1,875 / 150 | 234 |
| 5 | Chained S3 + CryptoMiniSat (partial) | 63,567 | 115,500 | 18,702 | 1,824 / 1,875 / 150 | 234 |
| 5 | Plain Pollard rho | 20,208 | 39,885 | 7,170 | 582 / 30 / 30 | 0 |
| 9 | Nagao: quadratic image | 2,029,215 | 1,851,174 | 949,377 | 10,479 / 12,636 / 330 | 426 |
| 9 | Nagao: coefficient pullback | 1,069,488 | 1,456,620 | 384,993 | 10,479 / 12,636 / 330 | 426 |
| 9 | Nagao: cached pullback | 835,764 | 1,103,148 | 377,877 | 10,479 / 12,636 / 330 | 426 |
| 9 | Symmetric S4 + CryptoMiniSat (partial) | 378,480 | 759,645 | 245,466 | 11,769 / 13,890 / 330 | 447 |
| 9 | Chained S3 + CryptoMiniSat (partial) | 378,480 | 759,645 | 245,466 | 11,769 / 13,890 / 330 | 447 |
| 9 | Plain Pollard rho | 264,393 | 472,146 | 198,690 | 2,037 / 360 / 30 | 0 |

All 3,306 decomposition calls completed without a timeout in this fixed toy campaign. SAT and Nagao can select different first valid relations; their later attempt paths therefore need not coincide. The generated full-DLP targets, supports and seeds remain matched.

## Support boundary and correctness

| Field bits | Rational factor abscissae M | Projected columns B | Exact covered affine targets | Exact coverage | Counting bound p_max | Expected-attempt floor B/p_max |
|---:|---:|---:|---|---:|---:|---:|
| 5 | 5 | 4 | 33/43 | 0.767442 | 1.0 | 4.0 |
| 9 | 11 | 10 | 440/507 | 0.867850 | 1.0 | 10.0 |

The attempt-count floor is deliberately distinguished from a total-operation floor. The counting bound saturates on both tiny fixtures. Rank filtering can only increase the required attempts. Every IC run achieved full column rank and a verified scalar. Independent post-run auditing took 5.548 seconds outside the solver timers.

## Evidence

- [Frozen contract](contract.json): accepted fixtures, variants, accounting and correctness gates.
- [Source and host provenance](results_v1/provenance.json): base commit and SHA-256 hashes of imported repository sources.
- [Raw main runs](results_v1/raw.jsonl): all 360 attempts with intermediate relation calls and scalar certificates.
- [Machine-readable summary](results_v1/summary.json): phase totals, vectors, paired wall-time diagnostics and audit.
- [Exhaustive validation](results_v1/validation.json).
- [Smoke evidence](smoke_v1/summary.json): the separate 12-run integration check.
- [Candidate review and reproduction instructions](README.md).

The implemented comparison is complete for its two fixed toy fixtures. A matched F4/WDSat comparison, automorphism-aware support handling, realistic sparse linear algebra, and calibrated total cost remain outside this result.
