# Results: degree-73 descent at 2^37

The SageMath 10.6 run constructed and checked a separable degree-73 isogeny defined over GF(2^37). It returned 74 base-field-defined cyclic kernels, with a squarefree degree-36 kernel polynomial dividing the degree-2664 division polynomial. The selected map has equal source/target group orders `137439487532`; its 73-kernel subgroup is Galois-stable, and its nonzero points are not assumed individually rational over the base field. The map includes exact rational functions in the raw certificate.

All 128 Boolean systems completed, with all planted signed decompositions verified on both source and codomain. The selected isogenous neighbor reduced aggregate F5 matrix XORs by at least 10% in 4 of the 8 registered seed/summand/encoding cells; the largest regression was 2.00x. Completion degree rose in 7 of 64 paired systems. The preregistered criterion therefore **failed**.

| Split | Summands | Encoding | Source F5 XORs | Degree-73 F5 XORs | Ratio | Source degree sum | Neighbor degree sum | Paired higher |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| development | 2 | ordered | 8243 | 9713 | 1.1783 | 32 | 32 | 0 |
| development | 2 | canonical | 8596 | 10079 | 1.1725 | 32 | 32 | 0 |
| development | 3 | ordered | 85864 | 71463 | 0.8323 | 51 | 50 | 1 |
| development | 3 | canonical | 88506 | 60881 | 0.6879 | 51 | 49 | 0 |
| holdout | 2 | ordered | 9781 | 8771 | 0.8967 | 32 | 32 | 0 |
| holdout | 2 | canonical | 10164 | 9196 | 0.9048 | 32 | 32 | 0 |
| holdout | 3 | ordered | 65073 | 130446 | 2.0046 | 49 | 55 | 6 |
| holdout | 3 | canonical | 67423 | 48716 | 0.7225 | 49 | 49 | 0 |

Every row is a Boolean matrix-stage diagnostic. The comparison does not solve a discrete logarithm and has no calibrated conversion to field operations. Full-DLP total cost, `S`, rho/floor ratios, and speedup remain null. The two deterministic repetitions are replay checks, not independent samples. No security or end-to-end attack claim follows.

## Reproduction and custody

Run the frozen contract with the command in `research/koblitz_isogeny_descent_37_20260925/README.md`. The raw result is `raw.json.gz` in this directory; its SHA-256 is `eb4773d556886672e8b80735fc486c16b73765fd5069ff81133100bfc3eafa90` and its compressed size is `40311` bytes. The summary is derived from those committed raw cases.

Provenance: workflow run [36087300313](https://github.com/aburan28/crypto/actions/runs/36087300313), head commit `2679fe8acec49750b1775c7dd8e85facc90bd7f4`; GitHub artifact 10843942663 has archive digest `sha256:7b5305016d2490b0cfac7accff90aa977ccf5a9e28daf95e341ebcf378126665`. Source hashes and the Sage map/kernel data are preserved inside `raw.json.gz`.
