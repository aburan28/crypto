# Proof map: exact four-byte root slots

Lean toolchain: `leanprover/lean4:v4.34.1` in `lean-toolchain`.
Run `lean RootFour.lean` from this directory and inspect `proof_check.log`.
The `#print axioms` audit reports only Lean's standard `propext` and
`Quot.sound` for the arithmetic decode lemmas; the closed capacity
inequalities use no axioms.

| Result | Lean declaration | Exact assumptions | Status |
| --- | --- | --- | --- |
| Selector bit survives encoding | `decode_selector` | Natural state and fingerprint, selector less than two | Checked |
| State ordinal survives encoding | `decode_state` | `state+1 < 2^29`, selector less than two | Checked |
| Fingerprint bits survive encoding | `decode_fingerprint` | `state+1 < 2^29`, selector less than two | Checked |
| Encoded occupied word differs from empty | `encode_nonzero` | Natural state, selector and fingerprint | Checked |
| Encoded word fits 32 bits | `encode_fits_u32` | `state+1 < 2^29`, selector less than two and fingerprint less than four | Checked |
| K3550 complete two-root bound fits `2^31` slots | `k3550_complete_key_bound_fits_fixed_slots` | `stateBound` as defined | Checked |
| K3550 state IDs fit 29 bits | `k3550_state_id_fits_29_bits` | `stateBound` as defined | Checked |
| K3550 arrays have payload below 14 GiB | `k3550_one_root_and_slot_payload_below_14_gibibytes` | 11 bytes per state and four bytes per slot | Checked |
| K3596 state IDs fit 29 bits | `k3596_state_id_fits_29_bits` | `stateBound` as defined | Checked |
| K3597 state ID exceeds 29 bits | `k3597_state_id_exceeds_29_bits` | `stateBound` as defined | Checked |
| K3596 arrays have payload below 14 GiB | `k3596_one_root_and_slot_payload_below_14_gibibytes` | 11 bytes per state and four bytes per slot | Checked |
| Full root-key equality preserves exact lookup and first witness | Open | Rust open-addressing insertion order and GF(2^83) canonical-root reconstruction | Formalization partial |

The checked decode arithmetic is a model of the slot bit layout. The source
uses `u32` bit shifts and masks, and compares the full canonical root after
the one-bit fingerprint. Rust tests compare the four-byte table with the
existing exact five-byte table on complete small indexes and probes.
The missing Lean bridge is a formal model of the Rust table transitions,
the source's canonical-root function, and the GF(2^83) representation.
The characteristic-two root reconstruction identity is recorded in the
parent experiment's `RootSum.lean` and `PROOF_MAP.md`.
