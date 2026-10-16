# Indexed five-byte table density proof map

The K2543 Tier 1 receipt recorded 536,737,856 exact root keys for
268,375,505 retained symmetric-S3 states. Thus the quotient keeps almost
**two distinct keys per state**. The optional `KIC_INDEXED5_LOAD75=1` mode
requests `ceil(4*S/3)` entries for `S` retained states. The constructor
allocates the next power of two above twice that request. The source test
checks exact first-witness semantics on six N83 representatives. The K2543
cloud control compares base hashes, root-key counts, and recovered scalars;
the large index's global witness map remains a source-level proof obligation.

| Result | Lean declaration | Exact hypotheses | Status |
| --- | --- | --- | --- |
| K2830 request fits the 2^30-slot bucket | `N83Indexed5Load75.k2830_request_fits_2pow30_slots` | Source-loop bound `S=83*K*(K-1)/2+42*K`; request `ceil(4*S/3)`; table rounds twice request up to a power of two | Checked |
| K3114 fits and K3115 crosses the 2^30-slot bucket | `N83Indexed5Load75.k3114_request_fits_2pow30_slots`, `k3115_request_crosses_2pow30_slots` | Same bound and request; these are exact adjacent allocation checks, not an RSS feasibility claim | Checked |
| K2830 state indices fit the slot's 29-bit field | `N83Indexed5Load75.k2830_state_id_fits_29_bits` | State count at most the source-loop bound | Checked, conditional on source-count bridge |
| K2830 state and slot payload is below 13.7 GiB | `N83Indexed5Load75.k2830_state_and_slot_payload_under_13point7_gibibytes` | 28 bytes/state and 5 bytes/slot; no allocator or other arrays | Checked, conditional on source-count and allocator bridges |
| Exact key and first-witness equality under a changed slot count | Open | Source insertion visits the same `(key,state,selector)` stream; probing terminates with a free slot; both lookups recompute the full canonical key | Source test checks the six-representative fixture; cloud replay checks counts and valid target relations; Lean source-level bridge partial |

Command: `lean DenseCapacity.lean` from this directory with the pinned Lean
4.34.1 toolchain; `proof_check.log` retains its output. `#print axioms`
reports no axioms. `lean4checker --fresh` was unavailable on the development
host. The same-field source test compares every canonical key and first
witness for a six-representative N83 index at both densities, plus 2,047
additional lookup probes. The larger run must establish root-key count,
target relation, scalar, allocation, RSS, and timing. The 13.7-GiB bound
covers only the two array payloads. The K2543 measured overhead above those
arrays was about 0.297 GiB; K2830's total RSS is still a hypothesis to test.
