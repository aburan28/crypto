# What coefficient reuse proves and cannot fix

The circuit vector is quadratic in binary coordinates of the affine coefficient set {b : b²+b+r is in V}. Polarization therefore reconstructs it exactly from 1+k+k(k−1)/2 direct evaluations. Binary doubling materializes all 2^k vectors with (2^(k+1)−k−2)(d+1)² field additions. This is an algebraic identity, independently checked against direct circuits and exhaustive relation truth; it is not a change to the relation-yield bound.

On the frozen complete eight-target n30 d8 batch, all three pruned traversals execute the same 9,794,908 field additions in block rejection. Coefficient-only changes cannot remove any of them. This rejection phase alone is 16.8677 times the S3 total of 580,690 field API calls. The bound applies to the identical traversal and this uncalibrated component unit only.

Reuse with batch inversions reduces total multiplications from 655,536 to 87,552, but total API calls only fall from 10,980,028 to 10,495,212. S3 uses 207,401 multiplications and 580,690 API calls. A multiplication-only claim would hide the dominant work.

The next experiment caches T_k = span(Z_i, M_ij : i<k), independent of each high-bit prefix. Its exact equivalence and frozen test are in [rank_01](../rank_01/README.md). There is a second, smaller obstacle: unchanged coefficient-table setup already uses 588,888 API calls, slightly exceeding the whole S3 batch. Even perfect free elimination cannot make that unchanged setup meet the 20% API target. Any later attempt must reduce both preparation and elimination, and still account for binary work, verification and calibration.

Possible follow-ups, not results: exploit Z=W and mixed-matrix symmetry to store and interpolate fewer circuit entries; restrict z prefixes to abscissas that lift to the curve using a charged support table; or eliminate the fixed free-column span once and maintain the remaining equations directly in its quotient. Each requires a frozen matched rerun. The current search also has an auxiliary root witness w, so it is a function/support hybrid rather than a root-free solver. No asymptotic, calibrated-total, or full-DLP crossover has been established.
