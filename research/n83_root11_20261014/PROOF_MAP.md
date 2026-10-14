# Proof map: one retained N83 S3 root

Lean toolchain: `leanprover/lean4:v4.34.1`. The checked module is
`RootSum.lean`; run `lean RootSum.lean` from this directory and inspect
`proof_check.log`. Its `#print axioms` reports Lean's standard `propext`
and `Quot.sound` for the XOR lemma and no axioms for the closed arithmetic
inequalities.

| Result | Lean declaration | Exact assumptions | Status |
| --- | --- | --- | --- |
| Reconstruct second bit word from first and their XOR | `xor_reconstruct` | Arbitrary natural-number bit words | Checked |
| K3400 complete two-key bound fits `2^30` slots | `k3400_complete_key_bound_fits_fixed_slots` | `stateBound` as defined | Checked |
| K3400 state ID fits 29 bits | `k3400_state_id_fits_29_bits` | `stateBound` as defined | Checked |
| K3400 one-root and table arrays fit below 10 GiB | `k3400_one_root_and_slot_payload_below_10_gibibytes` | 11 bytes per complete state, five bytes per fixed slot | Checked |
| S3 root difference is `xy/(x+y)^2` on the regular branch; exceptional roots coincide | Open | GF(2^83) field laws, source solver `prepare`/`finish` semantics, nonzero denominator | Written proof; formalization partial |

For the regular branch, let `a=(x+y)^2`, `p=xy`, and `ps=p^2`.
The solver forms `q=p*ps/(a*ps)` where `a` and `ps` are nonzero.
Field cancellation gives `q=p/a`. Its returned pair is `[first,first+q]`,
and characteristic two makes `+` equal XOR in the normal basis.
For `a=0` or `ps=0`, the source delegates to the direct solver, whose
corresponding branches return `[root,root]`. Thus the retained first root
and state coordinates reconstruct the same ordered pair in every branch.

The field statement has not been kernel checked because this standalone
Lean module has no finite-field model or bridge from the Rust bit operations
to GF(2^83). The exact missing obligation is a formal definition of that
field representation and the `prepare`/`finish`/direct-solver equations.
Rust tests compare all roots in the six-representative index and exact
lookup parity, and the cloud replay checks the complete base and target.
