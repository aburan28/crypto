# Stage 14: orbit-derived cofactor admission

The public cofactor gate previously evaluated `[r]P` for every factor-base point. The replacement performs one scalar multiplication per signed Frobenius orbit, then derives the remaining torsion classes by Frobenius and negation. A dedicated degree-15/23 unit test confirms that the resulting distinct class set equals the full pointwise set exactly. Existing reachability tests retain their exact admissible and inadmissible summand counts.

All degree-23 candidate polynomials, point counts, projected point counts, projected signed-Frobenius columns, cofactor verdicts, and selections match the frozen Stage 12 reference.

| Curve/candidate | Pointwise admission | Orbit-derived admission | Internal speedup |
|:--|--:|--:|--:|
| `K_0 [0,1]` | 1.396340 s | 0.031098 s | 44.90x |
| `K_0 [0,2]` | 1.653789 s | 0.025585 s | 64.64x |
| `K_1 [0,1]` | 1.171100 s | 0.017853 s | 65.60x |
| `K_1 [0,2]` | 1.008953 s | 0.017415 s | 57.94x |

Across both fresh degree-23 discovery processes, total core time fell from the Stage 12 value of 5.537515 seconds to 0.700330 seconds, a 7.91x incremental reduction. Relative to the original Stage 10 processes at 10.857294 core-seconds, projection indexing plus orbit-derived admission provide a 15.50x cumulative process-level reduction. The selected bases and all public algebraic evidence are unchanged.

This optimization reduces charged factor-base discovery. It does not reduce SAT relation-solving cost and does not change the non-SOTA conclusion.
