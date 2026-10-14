# N83 guarded root-only proof map

| Result | Lean declaration | Hypotheses and status |
| --- | --- | --- |
| Diagonal and off-diagonal offset decoding | `N83RootOrdinal.diagonal_offset_roundtrip`, `off_diagonal_offset_roundtrip` | Given the correct left-prefix ordinal, `relative<42` on the diagonal and `<83` off it; checked. |
| K3114 fits the `2^30`-slot request, K3115 crosses it | `N83RootOrdinal.k3114_request_fits_2pow30_slots`, `k3115_request_crosses_2pow30_slots` | Exact source-loop bound and `ceil(4*S/3)` request; checked. |
| K3114 state IDs fit 29 bits | `N83RootOrdinal.k3114_state_id_fits_29_bits` | Every state is retained, so count equals loop bound; checked conditional arithmetic. |
| K3114 two-array payload is below 13.3 GiB | `N83RootOrdinal.k3114_root_and_slot_payload_below_13point3_gibibytes` | 22 bytes/state and five bytes/slot; checked. Other memory is measured. |
| Every source pair has S3 roots and the binary search finds its exact left label | Open | The producer asserts root presence and total count; source tests enumerate six real N83 representatives and compare exact keys and first witnesses. A faithful Lean model of curve addition, S3 roots, and the source traversal/binary search is still needed. |

`RootOrdinal.lean` uses pinned Lean 4.34.1. `proof_check.log` retains its
successful build and `#print axioms` audit. The offset lemmas use only
standard `propext`, `Classical.choice`, and `Quot.sound`; numeric bounds use
no axioms. There are no `sorry`, `admit`, claim axioms, or native checks.
