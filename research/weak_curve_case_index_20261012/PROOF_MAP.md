# Proof map for the published-case index

| Mathematical statement | Source or Lean declaration | Exact scope and status |
| --- | --- | --- |
| `X^3+1=(X^2+X+1)(X+1)` over `F_2[X]` | `PublishedCases.cubic_frobenius_factorization` | Lean-checked; supports the cubic `n=2` mask polynomial, not a curve covering |
| `X^7+1=(X+1)(X^3+X+1)(X^3+X^2+1)` over `F_2[X]` | `PublishedCases.septic_frobenius_factorization` | Lean-checked; supports the degree-seven factor choices, not a curve covering |
| Binary mask rule `Phi=a Fhat` modulo `X^d+1` and cyclic equivalence | [IMC 2015, Sec. 3.3](https://eprint.iacr.org/2015/805.pdf), native enumerator | Published structural rule; Rust replay checks the finite table, Lean formalization of quotient by rotation remains open |
| Existence and completeness of the 36 elliptic covering cases | [IMC 2015, Sec. 5](https://eprint.iacr.org/2015/805.pdf) | Published result under odd characteristic and prime extension degrees 2,3,5,7; full Lean proof remains open: formal definitions of `k_d(C)`, base-field covering descent, branch divisor, genus, and classification equivalence are missing |
| An exact model realizes a listed case and transports a selected subgroup | Per-model field, curve, covering, and map certificates | Open for the 36 case metadata rows in this companion index; do not infer from a binary mask |
| Odd cubic trace excludes a genus-three hyperelliptic Jacobian in the cubic restriction isogeny class | [Costa et al., Theorem 2.8](https://arxiv.org/pdf/2002.02067) plus the Weil polynomial `T^6-tT^3+q^3` | Source-derived conditional observation only; the hyperelliptic-Jacobian theorem and Weil-restriction bridge are not formalized here. The parity check excludes zero rows in the current even-trace corpus |

The Lean toolchain is `leanprover/lean4:v4.34.1`, with mathlib commit
`d13f23b723b8a846827a245b89c10fc7d3f11612` pinned in
[`lean/lakefile.toml`](lean/lakefile.toml) and
[`lean/lake-manifest.json`](lean/lake-manifest.json). The axiom audit reports
only `propext`, `Classical.choice`, and `Quot.sound` for both declarations.
No research-specific axiom, `sorry`, `admit`, or native-evaluation proof is
used. These checked lemmas are explicitly partial relative to the covering
classification.
