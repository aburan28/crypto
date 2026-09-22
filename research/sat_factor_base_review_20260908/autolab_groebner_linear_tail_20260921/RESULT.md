# Exact linear-tail F4 solver optimization

This iteration changes the work performed by the Boolean matrix-F4 splitting
solver without changing its ideal, roots, split tree, node budget, or external
gates.

For systems with at least 24 variables, the solver now builds only its maximum
configured Macaulay degree. That matrix includes all rows from lower degrees.
It then echelonises only the nonlinear columns, extracts the exact intersection
of the row space with the linear/constant tail, fully reduces that tail, and
materializes only `1`, `v`, or `v+1` consequences. The public full-RREF F4 API
is unchanged, and environment controls retain the previous solver behavior.

## Final same-binary paired result

Both arms used the same final executable, two public targets, one thread, and a
50,000-node budget. Every target passed the independent enumeration and group
gates.

| Measurement | Legacy | Selected | Ratio |
|---|---:|---:|---:|
| Complete process wall | 218.436 s | 150.964 s | 1.447x |
| Child user CPU | 215.335 s | 149.993 s | 1.436x |
| x-chained found median | 54.359 s | 27.696 s | 1.963x |
| Symmetrised refuted median | 31.013 s | 15.534 s | 1.996x |
| Peak RSS | 445,349,888 B | 462,749,696 B | 1.039x selected/legacy |

Both policies found both x relations and refuted both symmetrised systems. The
efforts remained 7,056 and 4,095, built degree remained 3, FFD remained 3 and
4, and there were no inconclusive outcomes.

## Mechanism profile

On the retained symmetrised solver profile, F4 calls fall from 16,382 to 8,191,
word XORs from 6,595,269,352 to 4,509,530,169, and nonlinear-row readback from
8.888 seconds to 7.2 milliseconds. Solver wall falls from 30.723 seconds to
15.951 seconds with the exact same 4,095 splits and 4,096 infeasible leaves.

On the x first-root profile, calls fall from 2,055 to 1,028 and word XORs from
853,440,117 to 584,547,728. Wall falls from 4.033 seconds to 2.189 seconds with
the same root, 515 splits, 19 propagations, and 511 infeasible branches.

Random dense, sparse, and rank-deficient matrices across word boundaries give
the same reduced linear-tail row space as full RREF. Forced-path exhaustive
small-system tests and every split rule return the exact brute-force root set.

The current-binary frozen ladder preserves all six verdict digests and is
faster on every rung, from 1.016x at `K_0/2^13, m=2` to 1.504x at
`K_0/2^9, m=3`. The untouched holdout ladder also preserves both digests:
`K_0/2^19` improves from 1.178 s to 0.593 s (1.986x), and `K_1/2^17`
improves from 0.444 s to 0.334 s (1.331x).

## Rejected directions

The root degree-4 sparse closure has 15,407 rows and 32,817 columns, takes 28.4
seconds, and neither refutes nor determines a variable. FFD 4 is not the
solving degree on this fixture. Forward-only M4RI widths 2, 3, and 4 reduce the
XOR count but lose on wall time to the simpler suffix-forward pass, so they are
not retained.

Classification: `N31_BOOLEAN_F4_LINEAR_TAIL_PAIRED_T2_PASS`.

This is bounded public-synthetic decomposition-solver evidence. It does not
measure relation yield, a complete index-calculus run, asymptotic complexity,
or a crossover against Pollard rho, and it is not a key-recovery result.
