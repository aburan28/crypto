# Correctness argument and test mapping

This is the producer's argument and validation map, not an external audit or a
machine-checked proof. The supported problem is a generated system of squarefree
Boolean polynomials of degree at most two. The syndrome path packs at most 32
equations and whole enumeration is bounded to 24 variables.

## Equation packing

Assign one u32 bit to each equation. XOR of coefficient words performs the same
GF(2) coefficient addition in every equation simultaneously. The resulting word
is zero exactly when all represented equations are zero. No equation is discarded
to fit the word: unsupported equation counts use the retained fallback. Tests
include all 32 bits, including bit 31.

## Gray differences

Split off four low variables. Their quadratic contribution is fixed, and their
linear coefficients depend affinely on the high variables. The remaining constant
term is a quadratic function of the high assignment.

For reflected Gray order g(t)=t XOR (t>>1), the transition from t-1 to t flips
j=ctz(t). When j>0, its lower Gray coordinates are zero except coordinate j-1,
which is one. Therefore the change in the high-only constant is

\[
D_j=\ell_j+Q_{j,j-1}+\sum_{i>j}Q_{ji}h_i,
\]

with the neighbor term omitted for j=0. Each D_j starts with its fixed linear and
neighbor terms. Flipping high variable j updates D_i only for i<j, because those
are the differences for which j is a higher coordinate. Each of the four low
linear coefficients changes by its cross coefficient with j. This maintains the
block's complete quadratic evaluation without treating the low variables as fixed.

The high Gray order is bijective, and each block visits all sixteen binary low
assignments. A completed whole enumeration therefore covers every assignment.
Reaching the point cap before the next complete block returns UNKNOWN, not UNSAT.
Widths below four use direct single-point evaluation.

## SIMD block and recovery

The four lanes encode the four assignments of the first two low variables. Four
groups encode the other two low variables. Each lane XORs the same constant,
linear contributions and fixed low quadratic contribution as the scalar path.
NEON uses unsigned minima to detect any zero syndrome; SSE2 uses equality masks.
Both read only complete four-u32 arrays. On a hit, the scalar block determines
the first satisfying low assignment and checks the wrapping sum.

The checksum folds block indices and u32 wrapping sums. It is deliberately called
a diagnostic: equal sums or equal hashes do not by themselves prove equal point
evaluations. `every_gray_block_matches_direct_quadratic_evaluation` checks every
point through twelve variables against direct evaluation. Separate tests force
each possible hit lane and include full-width u32 values and wrapping sums.

For a tree leaf, only variables occurring in its residual equations are enumerated.
Their labels are recorded in increasing order, and the local assignment is scattered
back into those labels and combined with the prefix assignment. Nonoccurring free
variables may be zero. Different tree branches fix disjoint prefixes, so their
enumerated assignment sets cannot overlap. Every returned model is evaluated on
the original equations by the benchmark driver.

## Other correctness controls

- Disabling the leaf cutoff reproduces the retained packed solver's model,
  logical counters and full trace.
- Scalar and SIMD members of each enumeration policy match models, assignment
  counts, block counts, prefix traces where applicable, and checksums.
- The untraced packed accounting control matches the original model and logical
  counters. Its absent trace is null; it is never presented as trace evidence.
- Affine substitution is checked over 524,288 small polynomial/map combinations
  against a truth-table/Mobius oracle, with separate recovery-composition checks.
- Whole-word affine products are compared over all 5,592,404 argument pairs in
  dimensions zero through ten.
- Degree-two products are admitted only after exact canonical-degree checks in
  the reference. The optimized common-variable rule is checked independently.
- Dense tests exercise coefficient-word boundaries through 36 variables.
- SAT assignments are checked directly. Benchmark UNSAT results require agreement
  with the completed retained search reference. Capped results remain censored.

The performance claim, if its gates pass, concerns complete solves of the frozen
generated systems on the recorded host. These checks do not establish production
index-calculus performance, a generic exponent change or a rho crossover.
