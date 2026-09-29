# Shape-selected M4RI for Boolean Koblitz F4

This package measures the row-reduction change in the Boolean matrix-F4 engine
used by the Koblitz Semaev decomposition path. It is based on the current
`origin/main` implementation, including active-suffix and sparse elimination.

## Selected policy

Macaulay matrices with at least 128 rows and 256 columns use block-4 Method of
Four Russians elimination when `columns <= 4*rows`. Other shapes retain the
active-suffix reducer. This routes the dense symmetrised degree-3 matrix through
M4RI while preserving the existing reducer for the wider x-chained matrix.
Fully assigned search leaves are checked directly against the untouched input
equations instead of constructing another matrix.

The frozen reference test checks byte-exact RREF output and identical rank on
dense, sparse, and rank-deficient matrices across 64-bit word boundaries. The
small-system solver test compares the complete root set with exhaustive
evaluation.

## Repeated kernel result

The public-synthetic fixture is `K_0/F_2^31`, dimension 16, two summands, and
published target scalar 66,142. Each row below contains 31 uncached degree-3
matrix constructions and reductions.

| System | Matrix | Active suffix | Selected | Counted word XORs | Result |
|---|---:|---:|---:|---:|---:|
| x-chained | 1,023 x 4,369, rank 1,022 | 8.189 ms | 7.960 ms (active suffix) | 8,930,014 | unchanged by policy |
| symmetrised | 992 x 3,872, rank 992 | 6.866 ms | 5.870 ms | 6,973,621 -> 3,700,406 | 1.170x median; 46.9% fewer XORs |

Block widths 4, 6, 8, and 10 are retained as 11-repeat ablations. Width 4 was
selected because it gives the lowest symmetrised median in that panel and keeps
the combination table small.

## Paired solver control

The suffix and selected-policy runs produced the same finite outputs:
x-chained `FOUND`, symmetrised `REFUTED`, efforts 7,056 and 4,095,
built degree 3, first-fall degrees 3 and 4, zero inconclusive outcomes, and
passing independent gates.

The complete selected process took 164.563 seconds versus 173.841 seconds for
the suffix control, a 1.056x result in this pairing. The symmetrised arm was
29,201 ms versus 34,196 ms, while the unchanged x arm was 55,675 ms versus
57,676 ms. One target and one process pairing do not establish an end-to-end
solver speedup. The retained result is the repeated symmetrised kernel
improvement and exact operation-count reduction; whole-solver timing remains
indeterminate until a multi-target paired distribution is collected.

Classification:
`N31_BOOLEAN_F4_M4RI_KERNEL_PASS_PAIRED_WALL_INDETERMINATE`.

This is bounded public-synthetic decomposition-kernel evidence. It does not
measure relation yield, rank accumulation, a complete index-calculus run,
asymptotic complexity, or a crossover against Pollard rho, and it is not a
key-recovery result.
