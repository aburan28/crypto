# Capacity proof map

The [Lean module](OnlineRangeCapacity.lean) uses pinned Lean 4.34.1 and
the natural-number state and slot formulas exposed by the pinned Rust
source. `lean OnlineRangeCapacity.lean` completed with the recorded
[`#print axioms` output](proof_check.log); every declaration below reports
no axioms. `lean4checker --fresh` was unavailable locally.

| Result | Lean declaration | Hypotheses and scope |
| --- | --- | --- |
| Exact symmetric-state counts for K3450, K3584, K4096 | `k3450_state_count`, `k3584_state_count`, `k4096_state_count` | Arithmetic definition `83·K(K−1)/2 + 42K`; equality for those three integers. |
| K3450 and K4096 dynamic requests require the 2³¹ slot bucket | `k3450_dynamic_bucket`, `k4096_dynamic_bucket` | Request is `(4·states+2)/3` and the constructor allocates the next power of two at least twice the request. The Lean inequalities pin that bucket. |
| Dynamic K3450 and K4096 exceed the guard with a 1 GiB working reserve | `k3450_dynamic_exceeds_guard`, `k4096_dynamic_exceeds_guard` | Five bytes per slot, eleven bytes per state, 15.25 GiB guard, 1 GiB policy reserve. The theorem is an exact integer inequality for each K. |
| Fixed K3450 and K3584 fit the capacity model | `k3450_fixed_fits_capacity_model`, `k3584_fixed_fits_capacity_model` | Fixed 2³⁰ slots, same state bytes, guard and reserve. This establishes the model inequality; measured RSS must still pass the run guard. |
| Both fixed candidates satisfy the source's maximum-root assertion | `k3450_fixed_root_count_fits`, `k3584_fixed_root_count_fits` | `2·states < 2³⁰` for each K, matching the implementation's fixed-index assertion. |

The Rust bridge is the frozen source's symmetric `max_states` count,
`KIC_INDEXED5_LOAD75` request, `KIC_INDEXED5_FIXED30` fixed request, and
five-byte indexed table constructor. The proof assumes those implementation
branches are selected by the recorded flags; the workflow verifies the
source hash, flags and resulting table receipts. The 1 GiB working reserve
is an admission policy, not an empirical RSS theorem. Final memory
compatibility remains an empirical gate with process RSS and OOM receipts.
