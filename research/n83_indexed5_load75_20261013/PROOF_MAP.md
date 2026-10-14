# Indexed five-byte table density proof map

The optional `KIC_INDEXED5_LOAD75=1` request allocates the next power of two
above `2*ceil(2S/3)` slots for `S` retained symmetric-S3 states. The table
stores at most one exact canonical key per slot and checks the complete key
against the retained root before accepting a fingerprint match. It inserts
states in the same source order as the regular request.

| Result | Lean declaration | Exact hypotheses | Status |
| --- | --- | --- | --- |
| K3114 request fits the 2^29-slot bucket | `N83Indexed5Load75.k3114_request_fits_2pow29_slots` | Source-loop bound `S=83*K*(K-1)/2+42*K`; request `ceil(2S/3)`; table rounds twice request up to a power of two | Checked |
| K3115 crosses that bucket | `N83Indexed5Load75.k3115_request_crosses_2pow29_slots` | Same bound and request | Checked |
| K3114 state indices fit the slot's 29-bit field | `N83Indexed5Load75.k3114_state_id_fits_29_bits` | State count at most the source-loop bound | Checked, conditional on source-count bridge |
| K3114 state and slot array payload below 13 GiB | `N83Indexed5Load75.k3114_state_and_slot_payload_under_13_gibibytes` | 28 bytes/state and 5 bytes/slot; no allocator or other arrays | Checked, conditional on source-count and allocator bridges |
| Exact key and first-witness equality under a changed slot count | Open | Source insertion visits the same `(key,state,selector)` stream; probing terminates with a free slot; both lookups recompute the full canonical key | Source test and cloud replay gate; Lean source-level bridge partial |

The arithmetic declarations build with Lean 4.34.1 and `#print axioms` reports
no axioms. `lean4checker --fresh` was unavailable on the development host.
The same-field source test compares every canonical key and first witness
for a six-representative N83 index at both table densities. A larger exact
run must establish root-key count, target relation, scalar, allocation, RSS,
and timing. The 13-GiB bound covers only state and slot array payloads.
