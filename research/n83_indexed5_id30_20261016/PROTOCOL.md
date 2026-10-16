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
The executed pre-lint control and first K4096 scripts are preserved byte for
byte in `evidence_v1/` because the current helper now pins the later
lint-clean source. The original source is retained at commit
`5a88a28299de2b1ae82a3a1a6a6a1cd6ed2dc076`; the archive receipts
record the executed script and binary hashes. The historical materializers
verify those snapshots rather than silently reading the current helper.
Then run a bounded
K4096 feasibility diagnostic on a host with enough RAM, preserving source,
binary, flags, base, relation, scalar and resource receipts. Any speed claim
requires a paired Tier 1 online IC/rho panel with same-point replay and
the declared cloud-instance envelope.

The first K4096 diagnostic used source SHA-256
`f14d9747cd24627d39d7cae4c9bef69ed2819544125ef15f31306ce5d75f3fc7`.
The current-head confirmation uses source SHA-256
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`
after adding a lint-only dead-code attribute. Both use `run_k4096.py` with
the same public point and rank seed, `root11_dynamic`, and
`KIC_INDEXED5_ID30=1`. Each runs one
fresh process for at most 900 seconds. The state formula gives 696,256,512
retained states, 679,936 expected usable base points, a 10 GiB exact root
table and about 7.13 GiB of retained root states before other working
memory. Admit the local host only with at least 40 GiB physical RAM and 50%
reported free memory. Sample process RSS and memory pressure each second;
stop at 32 GiB RSS or below 15% reported free memory. Save the exact
source/binary/input hashes, flags, stdout, stderr, base dump, target record,
outer wall and sampled resource receipt, including timeout and resource
stops. `verify_k4096.py` independently replays the full base, target
relation and scalar after a completed producer run. This local diagnostic
tests feasibility and correctness; its wall time is descriptive.
