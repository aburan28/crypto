# State-ID capacity proof map

The [Lean module](StateId30Capacity.lean) uses the pinned Lean 4.34.1
toolchain. `lean StateId30Capacity.lean` completed with the recorded
[`#print axioms` output](proof_check.log). The two concrete boundary
declarations have no axioms. The structural and slot theorems use only
Lean's standard `propext`, `Quot.sound` and `Classical.choice` foundations
through `omega`; no custom
axiom, admitted proof or unverified native evaluation is used.
`lean4checker --fresh` was unavailable locally.

| Mathematical result | Lean declaration | Exact scope |
| --- | --- | --- |
| K5086 has fewer than 2³⁰ symmetric N83 states | `k5086_fits` | Exact natural-number formula `83*K*(K-1)/2 + 42*K` at K5086. |
| K5087 reaches the state-ID limit | `k5087_exceeds` | Same exact formula at K5087. |
| Every K ≥ 5087 exceeds the 30-bit state-ID field | `excluded_above_5086` | Uses monotonicity of both factors `K` and `K-1`; no premise about timing, memory, or relation yield. |
| The 30-bit ID, selector and nine-bit fingerprint decode exactly | `slot_roundtrip` | For every natural state with `state+1 < 2³⁰` and selector `< 2`; fingerprint is unrestricted for this arithmetic identity. |
| The encoded word occupies at most five bytes | `slot_fits_five_bytes` | Adds the source's fingerprint `< 2⁹` bound and proves the packed integer is `< 2⁴⁰`. |

The source bridge is the opt-in encoder assertion
`state_index+1 < 2^30`, the same symmetric traversal's exact
`max_states` formula, and the flag check requiring one-root N83 states.
Tests cover the encoder/decoder boundary and full-root equality after
fingerprint matches. A future other-format index is outside this bound.
