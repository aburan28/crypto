# Four-byte ID30 proof map

| Claim | Lean declaration | Exact hypotheses | Verification |
| --- | --- | --- | --- |
| The packed word is nonzero and fits `u32` | `N83FourByteId30.encode_nonzero_and_fits` | `state + 1 < 2^30`, `selector < 2`, `fingerprint < 2`, all natural numbers | Lean 4 v4.35.0 kernel checked; `#print axioms` records standard foundations |
| Low 31 bits and high bit recover the three fields | `N83FourByteId30.decode_roundtrip` | `state + 1 < 2^30` and `selector < 2`; the division identity holds for every natural fingerprint | Lean 4 v4.35.0 kernel checked; `#print axioms` records standard foundations |
| A power-of-two table with at least twice the ceiling request has at most 75% root-slot occupancy when each state contributes at most two roots | `N83FourByteId30.load75_bound` | Natural `states`, `slots`, and `2 * ((4 * states + 2) / 3) <= slots` | Lean 4 v4.35.0 kernel checked; `#print axioms` records standard foundations |
| N83 symmetric-state count fits the 30-bit ID mode through K5086 | `N83RootIndexCapacity.id30_exact_ceiling` in the parent proof module | `83*K*K + K < 2^31` for natural K, equivalent to K≤5086 | Parent Lean 4 v4.35.0 kernel checked; source-side bound audited separately |

The proof concerns the packed word and integer capacity. The full-key lookup,
first-witness behavior, field arithmetic, curve relation and scalar replay
are checked in the Rust tests and independent benchmark verifier.
