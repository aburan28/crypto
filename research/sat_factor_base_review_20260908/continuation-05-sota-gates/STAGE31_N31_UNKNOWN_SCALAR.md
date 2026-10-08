# Stage 31: n=31 public unknown-scalar workflow

Stage 31 is gated on a successful Stage 30 factor-base artifact. It accepts only a schema-v1 degree-31 `K_0` divisor-kernel recipe or its explicit two-torsion saturation with the divisor retained as its parent. It does not select or modify the factor base and cannot fall back to a default recipe.

The run precomputes factor-base logarithms from public probes and collected relations. Sparse linear algebra must determine every projected column, and every computed log is certified by its group identity before descent begins. Thus factor-base logs are learned by the algorithm and are not known by construction.

Five domain-separated public targets use seeds 31001 through 31005. Hash-to-curve and cofactor projection construct their points without constructing or supplying target scalars. A descent result is accepted only after `[d]G = Q`; the same public target is then passed to the signed-Frobenius rho baseline, which uses the same group-identity check.

The wrapper charges the full select/spec-load, relation collection, pair table, sparse linear algebra, five descents, rho controls, and persistence envelope. It records child and outer CPU, wall time, average parallelism, process RSS, and sampled process-tree RSS; the exact binary build is retained separately.

This is a finite public n=31 experiment. Even a successful run does not supply n=41 end-to-end index calculus, licensed Magma, external reproduction, novelty, or SOTA evidence.
