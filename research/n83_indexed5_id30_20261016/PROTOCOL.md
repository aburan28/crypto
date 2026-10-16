# N83 opt-in 30-bit exact root-state IDs

The five-byte indexed root slot currently uses 29 bits for the state ID,
one root-selector bit and a ten-bit fingerprint. This candidate moves one
fingerprint bit to the state ID under `KIC_INDEXED5_ID30=1`. It retains the
same five-byte slot size, the same zero-empty marker, and full canonical-root
equality after every fingerprint match. The flag requires N83 symmetric
one-root states, indexed five-byte slots and the dynamic load-75 request.
With the flag absent, the prior encoding and source behavior remain selected.

The N83 symmetric state count is `83*K*(K-1)/2 + 42*K`. The new slot admits
`state_index+1 < 2^30`, so K5086 is the largest admissible folded-column
count for this exact format. The [Lean module](StateId30Capacity.lean) proves
the boundary and its implication for every K ≥ 5087. This is a state-ID
encoding bound; available RAM, table capacity, full relation collection,
and verified online timing are separate empirical gates.

Test old and new slot round trips at their boundaries, exact first-witness
behavior under fingerprint collisions, and equal recovered scalars and base
hashes on the same exposed public point for K1024. The correctness control
uses the public point with SHA-256
`2c068d2f659588edcdbaa685fc73bb146e2bfee3be123d377d537658f8021c6d`,
rank seed 29, two fresh processes in legacy-then-ID30 order, a 900-second
cap per process, and the source SHA-256
`f14d9747cd24627d39d7cae4c9bef69ed2819544125ef15f31306ce5d75f3fc7`.
Each producer receives the public point only. The verifier independently
replays the factor base, target relation and scalar and requires equal base
hashes and exact-index entry counts. Shared-host timing is descriptive.
Then run a bounded
K4096 feasibility diagnostic on a host with enough RAM, preserving source,
binary, flags, base, relation, scalar and resource receipts. Any speed claim
requires a paired Tier 1 online IC/rho panel with same-point replay and
the declared cloud-instance envelope.
