# Wanted thin products for ECDLP research

## Source
Alman and Vassilevska Williams, arXiv:2610.06783 (2026), introduces a sparse-output thin matrix-product primitive: compute only a prescribed set W of entries of X*Y rather than the full product. The theorem should be treated as an asymptotic research lead, not a claimed practical speedup.

## Crypto application rule
When an expensive crypto kernel performs a large family of pairwise/bilinear compatibility computations, test whether it can be represented as

    X (N x D) * Y (D x N), restricted to wanted outputs W.

Do not count a reformulation as a speedup. Benchmark end-to-end cost, including feature construction, filters, misses, duplicates, verification, and rank updates.

## ECC2K-83 / ECC2K-130 experiments

### A. Relation residual matching (priority 1)
Split candidate decompositions into left/right partial decompositions. Construct compact feature vectors x(P), y(Q) encoding necessary algebraic compatibility. Cheap filters (Frobenius orbit, trace, support, phase, degree) define W. Evaluate only W.

Acceptance metric: total wall-clock time per new independent verified relation against the current matcher on identical targets and factor base.

### B. Degree-19 residual pipeline
Instrument residual generation separately from matching. Test whether residual predicates admit bilinear feature maps. Report N, D, |W|/N^2, feature-build time, product time, candidate count, verification time, and independent-relation yield.

### C. F4/F5 and Macaulay sparse algebra
Search multiplication/reduction stages for thin bilinear products where only selected monomials/rows are consumed. Do not apply the matrix method if materializing feature matrices or W dominates.

### D. Reduction graph
Maintain a machine-readable reduction graph connecting ECDLP, point decomposition, summation polynomials, sparse polynomial solving, residual matching, triangle/set problems, and matrix-product primitives. New algorithmic results should be propagated backward through this graph to identify accelerated kernels.

## ECC2K-130 extrapolation
Only extrapolate after a measured win on valid smaller instances (prefer m=83; m=51/53 may be used for scaling). Record exponents in powers of two where useful. For m=130/131, report memory, arithmetic count, sparsity, and expected verified-relation throughput. Never silently convert a microbenchmark win into an ECDLP work-factor claim.

## Experimental obligations
1. Compile and run the direct control baseline.
2. Same factor base, targets, timeout policy, and verification.
3. Include failed searches and preprocessing.
4. Count duplicates and rank contribution.
5. Report total cost per verified independent relation.
6. Separate theorem applicability from heuristic analogy.
7. A negative result is a result; do not weaken requirements to claim success.

## Implementation sketch
Introduce an experimental `wanted_thin_product` interface with:
- feature builders for left/right candidates;
- a sparse wanted-pair set W;
- scalar/reference backend first;
- optional optimized backend only after equivalence tests;
- instrumentation for N, D, |W|, arithmetic work, memory, and timings.

The first deliverable is a correctness-preserving backend and benchmark harness, not a claim that the asymptotic paper already improves ECC2K.
