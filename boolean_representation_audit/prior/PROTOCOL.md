# Frozen protocol: small Boolean closure audit

Frozen before implementation and measurements, 2026-09-28.

Scope: standalone synthetic finite Boolean algebra, at most 10 variables. No curve inputs, discrete-logarithm pipeline, or cryptanalytic integration. This is a correctness experiment and a candidate adaptation of the recovered storage component, not a claim of mathematical novelty or a scalable solver.

Hypothesis: replacing a lossy degree cap with a degree-prioritized queue of exact products preserves the generated Boolean ideal. Generating variable multiples only from accepted basis rows can avoid some redundant submissions relative to exhaustive monomial multiplication. It may still increase runtime or memory.

Ring: F2[x1,...,xn]/(xi^2+xi). Monomial multiplication is bitwise union; repeated terms cancel. All degrees through n are retained. Zero coordinates are allowed.

Reference: enumerate every monomial multiple of every input generator, using the recovered sparse echelon component at its full degree n. This is deliberately simple, exact, and not a state-of-the-art performance reference.

Candidate: degree-prioritized completion. Submit the original generators, then variable multiples of every newly independent row. Retain origin records. Success requires an independent dense-bitset verifier checking provenance, generator containment, and closure under every variable. Budget exhaustion means incomplete.

Independent controls: exhaustive truth-table solutions including zero coordinates; the Boolean ideal's vector-space dimension must be 2^n minus the number of solutions; compare the complete candidate and reference spans. Mutation tests must reject forged provenance and missing closure. Seeded random correctness controls use n=1,...,6, ten systems per n, fixed seed 20260928.

Measurement suite: n=6 and 8; monomial, block-linear, planted sparse-quadratic, and dependent/contradictory families. Fixed seed 131. Each algorithm runs three times with no tracing for median wall time. One separate traced run records peak Python allocations. Record submitted/accepted rows, XOR steps, stored terms, peak row terms, pending terms, and maximum scheduled degree. All runs use identical equations, row representation, reducer, and resource budgets. Include both positive and negative results.

Budgets: 100,000 terms per row, 1,000,000 total pivot terms, 1,000,000 queued terms, and 100,000 submitted rows. At most 10 variables is a deliberate validation limit. Full Boolean ideal closure can have dimension 2^n; sparsity does not remove that exponential ceiling. A speedup against this elementary reference is not evidence of a faster general algebra algorithm.

Stop: report correctness failures before performance claims. Stop after the declared suite; no post-hoc tuning to obtain a favorable result. Prototype completion is not an ECC attack exponent result.
