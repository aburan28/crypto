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

The [sealed native session](sessions/mac_arm64_l0_01) is
`ECBS1he4e70185e6e2`, measured by binary SHA-256
`1a21947f0501585e4e3784bc43e3ac635cb78abe8b15a4dacbecf9e5c5166922`.
It contains 384/384 verified recoveries, including 320 measured executions.
Its records SHA-256 is
`60838cb661b6df9ba801cca089975ce0d6bca259e1bb18f148461340fbe0f169`.
The [local audit](AUDIT.json) replayed all 320 measured executions exactly;
its SHA-256 is
`cdee74be352be7b3939d6efd3575386c4b2ef7afe944f58836325c3230b3bbb4`.
The [frozen first-round jobs](JOBS.json) cover sequences 64–127 and have
SHA-256
`a2d63e9063c98955e7fdf49b137d1c9baafc9912b903f8aa79cd69a43035a32f`.
Every Mac run earned L0, so wall ratios from this session are descriptive.

Independent Linux replay, Callgrind profiling, and the decision are pending.
The
[preceding cold-solve result](../ecbench_n37_k8_k16_20261004/RESULT.md)
selected K8 in complete-solve Ir but left the primary online base choice
open.
