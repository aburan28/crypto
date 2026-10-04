# n37 target-only instruction attribution

This directory holds the [preregistered protocol](PROTOCOL.md) and
[frozen ecbench spec](SPEC.json) for a fresh K8/K16 target-only
instruction comparison. The protocol and spec were pushed in `a2f9edc6`
before target planning or measurement. The subsequent [plan](PLAN.json)
contains 16 distinct public one-target workloads and 384 executions. Its
SHA-256 is
`8777d3ac7d853c57611294700d497c88bc7cd2b6b2e1f6d859fecf4922a85c72`;
the spec SHA-256 is
`655b05ceda0ba07226d04cf3bb21140b0e9f36b49f7237ddeb2d5333ce40821c`.
Its workload IDs and exact encoded target points have empty intersections
with both the preceding 16-target panel and the older eight-target sweep.

Native execution, independent replay, Callgrind profiling, and the decision
are pending. The
[preceding cold-solve result](../ecbench_n37_k8_k16_20261004/RESULT.md)
selected K8 in complete-solve Ir but left the primary online base choice
open.
