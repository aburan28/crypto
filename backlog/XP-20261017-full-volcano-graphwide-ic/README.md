# Full graph-wide IC campaign

Run the full matrix in [the protocol](PROTOCOL.md), preserving every vertex,
conductor, independent one-target workload, factor base and method, including
failure, timeout and unsupported cases. This card is a draft, not a launched
or leased job. [Certified input preparation](../../research/full_volcano_graphwide_ic_20261017/README.md)
covers all 24,823 vertices of the three additional completed graphs. It does
not discharge the original fixed A--E campaign.

Before ready: register exact per-vertex subgroup/generator representations;
validate subgroup/quotient-aware IC relations and counted factor-base setup;
freeze and audit a native ecbench control/pilot; derive total core, memory,
storage and time needs from that pilot; select/provision a worker and durable
checkpoint/output storage; implement and test the exact full-matrix driver.
The existing k=2 supersingular MOV teaching implementation is not finite-field
IC for the actual large embedding degrees in the three classes.

Acceptance: a per-vertex ledger matches every planned row with verified
recovery or a specific retained failure status; full exclusive cold/online
costs, raw relation yields/ranks and min/max/mean distributions are auditable;
all deterministic records replay and every recovered scalar is independently
checked. No timing is qualified without the required isolation/controls, no
generic fallback is called IC, and no sample is called full graph coverage.
No resource estimate is recorded yet; the 128 GiB disk entry is a proposed
worker envelope, not an allocation or a measured requirement.
