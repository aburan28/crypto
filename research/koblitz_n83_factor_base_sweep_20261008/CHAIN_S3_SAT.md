# Higher-arity chained-S3 SAT source gate

`binary_semaev_chain_sat.rs` adds a full-width native-XOR circuit for the
binary-curve relation `S3(x,y,z)=xyz+(xy+xz+yz)^2+b=0`. Four such nodes encode
five summands and five nodes encode six summands. Every field multiplication
is a set of two-input AND gates with exact polynomial-basis reduction; the
linear output and S3 equations are native XOR rows. Rewriting
`xy+xz+yz = xy+(x+y)z` needs only three field products per affine S3 node;
the square is a linear XOR transform. The module accepts field
degrees through 127, including the study's degree-83 polynomial, but this
change constructs **no retained N83 instance**.

For `m` summands, there are `m-2` internal point sums. The identity has no
x-coordinate, so the producer searches all `2^(m-2)` intermediate-identity
patterns: eight for m=5 and sixteen for m=6. An affine node gets an S3
equation. An identity output requires equal input x-coordinates; adding an
affine summand to an identity input copies its x-coordinate to the output.
The impossible identity-plus-affine-to-identity case is refuted immediately.
Every SAT model is re-evaluated over the original binary field, checked
against the exact finite factor-base coordinate set, and lifted over all
available point signs in the original curve group. Internal S3 roots or
trace-inadmissible algebraic solutions cannot become accepted relations.
If a model fails group lifting, all of its intermediate x assignments share
the same summand tuple, so one blocking clause excludes that tuple. A
refutation requires every pattern to finish; model, conflict, variable or
domain-clause caps return an inconclusive result.

The pre-domain all-affine source count at n=83 is 84,743 SAT variables for
m=5 and 105,908 for m=6. The respective circuits have 82,668 and 103,335
two-input AND gates, each represented by three CNF clauses, before exact
finite-domain constraints. These are construction counts, not memory or
solving-time observations. The finite-coordinate trie can dominate memory:
the producer computes its exact missing-prefix clause count on a given base
before installing any domain clauses and applies a caller-supplied cap. An
external hard memory and wall guard remains necessary for any retained-base
construction or search.

`sat_decompose_chained_s3` is a separate experimental library entry point
for ambient-basis explicit-orbit bases and m=5 or m=6. It is not silently
substituted for the frozen sweep's SAT disposition, and no retained-base
relation, natural rank, full cold time or winner follows from its small-curve
checks. The focused correctness suite also fixes both operands above bit 64
in the pinned degree-83 polynomial and compares every product-output bit
against independent field multiplication; it does not build a retained-base
model or time a relation search. Its next empirical gate is a source-pinned,
cgroup-guarded retained
model-construction receipt with the exact stored base domain. A subsequent
search must record unsuccessful patterns, solver work, original field/group
replay and rank novelty before the backend can enter a versioned sweep.
