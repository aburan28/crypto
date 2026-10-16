# Dynamic one-root table proof map

[`Root11DynamicCapacity.lean`](Root11DynamicCapacity.lean) checks the exact
arithmetic of the three dynamic candidates under the source traversal count
`S = 83*K*(K-1)/2 + 42*K` and request `ceil(4*S/3)`. The constructor doubles
that request and rounds up to a power of two.

| Result | Lean declaration | Hypotheses | Status |
| --- | --- | --- | --- |
| K1024, K1271 and K1536 state counts | `k1024_state_count`, `k1271_state_count`, `k1536_state_count` | Source symmetric traversal retains every state | Checked arithmetic; source bridge asserted and tested |
| Requests are 58,021,888; 89,388,583; 130,548,736 | `k1024_request`, `k1271_request`, `k1536_request` | Exact state count and integer ceiling formula | Checked |
| Slot buckets are `2^27`, `2^28`, `2^28` | `k1024_slot_bucket`, `k1271_slot_bucket`, `k1536_slot_bucket` | Twofold request followed by power-of-two rounding | Checked arithmetic; allocator bridge by source inspection |
| State identifiers fit 29 bits through K1536 | `k1536_state_ids_fit` | `stateCount + 1 < 2^29` | Checked arithmetic; encoder assertion remains active |
| Exact first-witness equality across capacities | Open | Same `(key,state,selector)` insertion stream, free slot, and full canonical-root equality on collision and lookup | Existing six-representative source test checks equality; complete source-level Lean model is pending |

The pinned toolchain is `leanprover/lean4:v4.34.1`. Build with
`lean Root11DynamicCapacity.lean` from this directory and audit its
`#print axioms` output in `proof_check.log`. The standard Lean foundations
may appear; the file uses no `sorry`, `admit`, custom claim axioms, or native
evaluation. The measured run must separately check actual allocation,
resource use, factor-base relation and scalar replay.
`lean4checker --fresh` was unavailable on the development host.
The `proof-check` CI job rebuilds these declarations with the pinned
toolchain and compares the `#print axioms` output with `proof_check.log`.
