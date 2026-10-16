# First retained-N83 mixed-orbit partial-search gate, revised baseline

The prior local experiments are conservatively charged at **12,600 seconds**.
The [900-second build/import result](../mixed-orbit-import-900s-20261012/RESULT.md)
accepted the exact nested K=1,182/K=2,048 S3-backed bases but did not start a
relation search. This protocol supersedes the earlier source-only protocol's
11,700-second accounting. It is **prepared, not authorized or launched**.

The current supervisor requires a new, run-specific
`n83.additional-local-budget/v2` record with purpose
`mixed_orbit_two_residual_first_gate`,
`prior_total_charged_seconds: "12600"`,
`prior_import_audit_sha256:
"76a8bf74e46feb0ef4fbe650395f39733f4ff95d31c35f8829cb506a9ad04c7d"`,
at least **2,400 newly authorized seconds**, exact fresh build and run
directories, and a user approval reference. The prior 900-second authorization
cannot start this gate. Both build and run recheck the committed, byte-pinned
import audit before Docker. The old v1 budget record is rejected.

The frozen first arm uses the primary `a=0` Koblitz curve, seed `2026100801`,
`public_x_hash` K=1,182 as the small signed-Frobenius base and the exact
K=2,048 envelope, leaving 866 residual signed orbits. The public fixture-zero
target has canonical corpus BLAKE3
`831eaa7be2387634fd610f6706103b2af9b4e274f0db59923c0280af3cbed74b`.
Search one deterministic equation with three small and two residual summands,
ordered first small pair, at most 150,000 variables, 12,000,000 domain clauses,
32 models and 5,000 conflicts. Independently derived first-trial coefficients
are `a=1833078656608906180176350` and `b=163089659398503559150141` modulo
the pinned subgroup order; the outer audit rejects a changed equation.

The small and envelope compressed S3 objects must have SHA-256 values
`eb342b574e244e1e67e8e78c805f33ee88edcdb4f3a14d4114e6f413beede380`
and `aacb62784dee0c9f5e06328dbb85a4b6a1c7243e778328da82900e7284f3a74a`.
Their exact S3 URIs are in the [original object table](../mixed-orbit-supervisor-source-20261012/PROTOCOL.md).
The image is pinned to Linux arm64
`sha256:7c4ae649a84014c467d79319bbf17ce2632ae8b8be123ac2fb2ea5be46823f31`.
Build is offline with four CPUs, 8 GiB and a 450-second outer wall, charged
600 seconds. Import has one CPU, 4 GiB, zero swap and a 360-second wall.
Only an accepted import receipt and verified container cleanup permit search.
The one-trial search has one CPU, 4 GiB, zero swap, a 1,200-second outer wall,
and a 1,150-second worker watchdog. Both runtime containers have no network.
A 240-second reserve covers cleanup, frozen-input checks, independent group
replay and receipts. Full conservative charge is 2,400 seconds; each launched
stage and the applicable reserve are charged even if the run fails early.
There is no automatic retry.

The entry point is
`research/koblitz_n83_factor_base_sweep_20261008/mixed_orbit_supervisor.py`.
Its `build` mode requires `--authorization`, `--build-dir` and `--output-dir`;
`run` also requires `--small-object` and `--envelope-object`. The same v2
authorization binds both modes. Preserve the authorization bytes, clean source
commit, source hashes, Linux ELF and binary hash, object hashes, image ID,
stage logs, cgroup and container states, worker config/result/events,
independent group replay, and outer and budget receipts. Keep timeouts, OOMs,
invalid models, malformed receipts and cleanup failures as separate outcomes.

This is a capacity and partial-yield gate. Its result cannot rank complete
cold runtime. A factor-base winner requires comparable candidates through
graph completion, independent rank growth, sparse linear algebra, an
individual logarithm, one complete cold single-target runtime and matched
rho with calibration and host-isolation evidence. The prior 900-second import
charge remains part of the historical experiment ledger and cannot be
silently omitted from any claimed cost of this line of work.
