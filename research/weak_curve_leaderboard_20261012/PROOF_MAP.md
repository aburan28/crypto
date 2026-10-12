# Lean proof map

The proof artifact is `lean/PrimeExamples.lean`. Its own `lean-toolchain`,
`lakefile.toml` and manifest pin Lean 4.34.1 and Mathlib commit
`d13f23b723b8a846827a245b89c10fc7d3f11612`. The declaration
`affineCount` is the finite sum of all
solutions to `y²=x³+ax+b (mod p)` for `0≤x,y<p`.

| Mathematical result | Lean declaration | Scope and status |
| --- | --- | --- |
| `E/F_101: y²=x³+x+32` has 100 affine solutions | `anomalous_F101_affine` | Exact finite sum, checked |
| `E/F_101: y²=x³+x` has 99 affine solutions | `smooth_F101_affine` | Exact finite sum, checked |
| `E/F_103: y²=x³+x` has 103 affine solutions | `pairing_F103_affine` | Exact finite sum, checked |
| First model has trace one | `anomalous_trace_identity` | Follows from finite sum and trace arithmetic, checked |
| For `r=13`, 103 has multiplicative order two | `pairing_embedding_degree_two` | Exact modular arithmetic, checked |
| All three displayed short-Weierstrass models are nonsingular | `anomalous_nonsingular`, `smooth_nonsingular`, `pairing_nonsingular` | Nonzero discriminant factors modulo the stated prime, checked |
| All three selected points satisfy the curve equation | `anomalous_point_on_curve`, `smooth_point_on_curve`, `pairing_point_on_curve` | Exact modular arithmetic, checked |
| The selected pairing order divides the group order | `pairing_selected_order_divides_group_order` | Divisibility arithmetic, checked |

The group-order formalization still needs the group construction with its
unique point at infinity and a
formal equivalence to this finite affine sum. Point order, Smart's trace-one
algorithm, Pohlig–Hellman, and the MOV/Frey–Rück pairing reductions require
their respective group or pairing formalizations. The present Lean proof is
therefore an exact arithmetic component, not a formalization of the full
cryptanalytic transfer.
