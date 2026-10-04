# Untouched K8/K16 confirmation

The [preregistered protocol](PROTOCOL.md) and [frozen ecbench spec](SPEC.json)
were pushed at `61f3946c` before the new target plan or any solve. The
[native plan](PLAN.json), generated afterward, contains 16 unique public
targets and 384 executions. Its 16 workload IDs and exact target points have
empty intersections with the preceding eight-target panel. The plan file has
SHA-256 `af6f80d8b3a9fdafdd3f78de8c062d7a4ab2268017e93fbdaac3abf3209afbf4`;
the spec file has SHA-256
`96dbd9b4b77ae876b057c5ae5e6bc121095aaa1a303a1d7513699e43c215fdca`.

The [sealed native session](sessions/mac_arm64_l0_01) is
`ECBS1hdebde0f96835`, measured with binary SHA-256
`356f0cfc7f53763bf9b4e4e0bf898f4a4fe860bbd267d9fca80d5ffe4cbd01e3`.
It contains 384/384 verified records, including 320 measured records. The
[local full audit](AUDIT.json) reproduced all 320 measured records exactly;
its SHA-256 is
`8fc88942c5aeb47022b325b690715e7e53b0900f5d481a3d7c8984ff512a6403`.
The sealed records have SHA-256
`d639e9a9ff8c5d152291bcae0683c070cd7a130b4359a27f59a25c5df06ca284`.
The [frozen first-round jobs](JOBS.json) have SHA-256
`36b15b44eb464957d3a1fd820a4522f314de61c05e1522605a912b5f72368b79`.

The saved counted-operation comparisons report K16/K8 = 1.0175
[0.973, 1.059], K8/rho = 4.5501 [4.110, 5.040], and K16/rho = 4.6298
[4.209, 5.095] on 80 pairs per arm. These are **incomplete cold counted-cost
diagnostics**. Every macOS run earned L0, so its wall measurements are
descriptive only. Independent Linux replay, whole-solve Callgrind profiling,
the predeclared K8/K16 decision, and the canonical scoreboard update remain
pending; no base is selected from this local result.
