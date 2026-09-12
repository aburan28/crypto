# Stage 37 algebraic-walk corroboration

Stage 37 preserves the second hosted execution of the Stage-35 algebraic walked descent. The repeat ran from commit `467b6b652c39730c6069442613e2397a6f173c48`, while the primary ran from `11dc73b8f1aed58d0d6f29d1d95bbc6a156a1387`.

The two commits are not whole-tree identical. Six unrelated quasi-subfield research paths were added between them. The verifier therefore calls this a Koblitz-path-equivalent corroboration and compares the Git objects for `Cargo.toml`, the frozen lock, Stage-35 parameters and runners, `src/binary_ecc`, the complete `src/bin/ic` tree, and the Koblitz factor-base, arithmetic, sparse-linear-algebra, and index-calculus modules. All those objects are identical.

The two outputs reproduce the same parameter digest, algebraic 4,759-point/29-column factor base, 39 group-certified factor-base logs, five public target coordinates and recovered scalars, 2,417,575 walked probes, and 425,616 rho iterations. Rho restart and group-addition counters also match exactly.

Timing varies while the narrow online result does not:

| Run | Source | IC faster than rho online | IC slower than rho amortized | Fresh build plus science / rho |
|---|---|---:|---:|---:|
| Primary `34693802792` | `11dc73b8` | 2.048767x | 28.857106x | 208.529527x |
| Repeat `34694400320` | `467b6b65` | 1.455658x | 32.136319x | 219.067345x |

The repeat supports a finite post-precomputation online crossover on project-authored hosted runners. It is not a whole-tree exact replicate, unaffiliated reproduction, full-cost crossover, asymptotic improvement, novelty verdict, licensed Magma result, or Koblitz index-calculus SOTA claim.

Custody:

- GitHub artifact `10297933814`
- Artifact digest `sha256:b1eec6f64f5dd70ee01a512892519f237397449a69c6970b1ba6dc96b7fbe744`
- Committed archive SHA-256 `f97578b9dbd9f4322fe57ed50f5920bc1bc95277897ea8452a7c929caa7af879`
- Result inventory SHA-256 `6c73c6e21b74170637877aa50a9b627d45646e645261122756489a3348f7106e`
