# K3114 dynamic-table boundary proof map

The [Lean module](Capacity.lean) is pinned to Lean 4.34.1 by
[`lean-toolchain`](lean-toolchain). Running `lean Capacity.lean` produced the
versioned [`#print axioms` output](proof_check.log). The seven concrete
arithmetic declarations have no axioms; the all-K implication uses Lean's
standard `propext` and `Quot.sound` foundations through `omega`. There are
no admitted proofs or research axioms. `lean4checker --fresh` is unavailable
on this host.

| Result | Declaration | Exact scope |
| --- | --- | --- |
| Retained-state and request counts at the boundary | `k3114_states`, `k3115_states`, `k3114_request`, `k3115_request` | Natural-number formula `83·K(K−1)/2 + 42K` and integer request `(4·states+2)/3`, for K3114 and K3115. |
| K3114 uses 2³⁰ five-byte slots | `k3114_uses_two_to_30_slots` | Twice the request lies in `(2²⁹, 2³⁰]`, matching the source's next-power-of-two allocation. |
| K3115 uses 2³¹ five-byte slots | `k3115_uses_two_to_31_slots` | Twice the request lies in `(2³⁰, 2³¹]`. |
| K3114 fits the current 29-bit state ID | `k3114_state_id_fits` | The 402,426,891 retained states are fewer than 2²⁹. |
| Every larger K crosses the 2³⁰-slot bucket | `every_k_above_3114_crosses_bucket` | For every natural K ≥ 3115, the exact state formula and integer request are monotone, hence twice the request exceeds 2³⁰. |

The readable proof of monotonicity is: for K ≥ 3115, both K and K−1 are
at least 3115 and 3114, so their product is at least 3115·3114. Multiplying
by 83, dividing by two, and adding 42K preserves the inequality for the
state count. Multiplying by four, adding two, and integer-dividing by three
preserves it for the request. Since K3115's doubled request already exceeds
2³⁰, every greater K does too.

The implementation bridge is the pinned Rust source SHA-256
`3ee6bc11399b382e72789b152043e5a67abc5b783c4e998ad8c534d60baebc03`:
its symmetric-state count uses the same N83 formula; `KIC_INDEXED5_LOAD75`
requests `(4·states+2)/3`; and its five-byte constructor allocates the
next power of two at least twice the request. Five bytes times 2³⁰ slots is
5 GiB; five bytes times 2³¹ is 10 GiB. This calculation fixes table-slot
capacity. A candidate still needs an actual base dump, memory receipt,
independent target replay and matched A/A-controlled timing.
