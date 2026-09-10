# Stage 12: indexed public projection map

The original public cofactor-projection pass found each projected factor-base point by scanning every signed-Frobenius representative. The replacement first identifies each orbit, builds one point-key hash table for all signed Frobenius translates, and then classifies every projected point by hash-table lookup. It uses the same public point operations and computes no discrete logarithm.

The verifier compared every degree-23 candidate field against the frozen Stage 10 outputs: divisor indices and polynomials, linearised exponents, dimensions, abscissa counts, rational-point counts, cofactor admission, projected point counts, projected signed-Frobenius counts, and selected candidates were identical.

| Curve/candidate | Reference projection | Indexed projection | Internal speedup |
|:--|--:|--:|--:|
| `K_0 [0,1]` | 1.841205 s | 0.145139 s | 12.69x |
| `K_0 [0,2]` | 2.716478 s | 0.166195 s | 16.35x |
| `K_1 [0,1]` | 1.553248 s | 0.080663 s | 19.26x |
| `K_1 [0,2]` | 1.504511 s | 0.064637 s | 23.28x |

Across both fresh discovery processes, total core time fell from 10.857294 to 5.537515 seconds, a 1.96x process-level reduction. The `K_0` process fell from 6.228588 to 3.258001 core-seconds; `K_1` fell from 4.628706 to 2.279514.

The matched degree-15 scalar-blind smoke retained the same target, predicate, seven relations, seven trials, 5,231 conflicts, zero invalid models, and verified scalar 53. Its fresh-process core time was 0.118330 seconds versus the frozen 0.203253-second reference. This single timing pair is supporting engineering evidence; the exact recovery counters are the correctness evidence.

This optimization reduces charged public setup cost. It does not accelerate SAT relation solving and is not a SOTA result.
