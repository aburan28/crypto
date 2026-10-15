# Pairing batch cardinality proof map

Lean 4 toolchain: `leanprover/lean4:v4.34.1`, pinned in `lean/lean-toolchain`;
no external Lean dependency. `lake build` compiles
`lean/TargetCardinality.lean`; its `#print axioms` output must show only the
standard Lean foundations, with no custom research axiom or placeholder.

| Written statement | Lean declaration | Status and exact assumption |
| --- | --- | --- |
| An order-19 group has at most 18 distinct nonidentity targets | `order19_batch_bound` | Conditional arithmetic proof: assumes the distinct-target count `n ≤ 19 - 1`; the finite-group counting bridge remains written mathematics |
| An order-61 group has at most 60 distinct nonidentity targets | `order61_batch_bound` | Conditional arithmetic proof: assumes `n ≤ 61 - 1`; same missing group bridge |
| Sixty-four distinct nonidentity targets require order at least 65 | `order_for_64_nonidentity` | Conditional arithmetic proof: assumes `64 ≤ r - 1`; prime-order eligibility strengthens this to the next available prime |

The unformalized bridge is: a finite group of cardinality `r` contains one
identity and `r-1` other elements. It is immediate by deleting the identity
from the finite set, but the present repository Lean project imports only
`Init`, not the finite-group cardinality library needed for a faithful group
declaration. The handwritten proof and the original requirement stay open for
formalization; the conditional Lean arithmetic theorems do not claim it.
