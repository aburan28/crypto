# Knowledge note: sparse-output thin matrix multiplication (arXiv:2610.06783)

**Status:** research lead; theorem-level result, practical crypto benefit unproven.

The 2026 Alman–Vassilevska Williams work accelerates computation of a prescribed sparse set of entries of a thin matrix product. The reusable research pattern is: identify a huge pairwise computation with a genuinely small middle dimension D and a sparse wanted-output set W, then share algebra across wanted outputs rather than evaluating pairs independently or materializing the full product.

For cryptanalysis, automatically inspect relation collection, residual matching, sparse polynomial solving, and linear-algebra kernels for this shape. Highest-priority application is ECC2K index-calculus residual matching. Any candidate must be compared end-to-end using total time per new independent verified relation; include preprocessing, misses, duplicates, verification, and rank updates.

Reference: https://arxiv.org/abs/2610.06783
