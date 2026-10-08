# Shared F6-IC pair index: operation cut did not reduce online time

The bounded follow-up in [PROTOCOL.md](PROTOCOL.md) completed on October 6,
2026. All 16 fresh-process executions exited zero, recovered and replayed
the same scalar for each paired public point, and had five exclusive online
phases whose sum equals the reported online interval. The target T1 point
is `(73407,129763)` with scalar 4785 and workload `146a1e9ee3c8`; T7 is
`(98625,98119)` with scalar 2391 and workload `ced1677f0976`. Both use
the exact 62-point usable n17 base, 29 folded columns, one prepared log
state, one target per process, and one Rayon thread. The same worker binary
`a6996e355c8fa52202845404f98cea426ca2d8c601baac52ea1c7ede81ab9afb`
ran all arms. [FREEZE.tsv](FREEZE.tsv), `candidates/`, and
[measurements.jsonl](measurements.jsonl) bind the exact candidate, workload,
run, source, input, and binary identities. All raw outputs, errors, exit
statuses, timestamps, and hashes are in `runs/`.

| Target | Rep | Inherited F4 online ms | Original F6-IC online ms | Attempt-local pair online ms | Shared pair online ms | F4 / shared | F6 / shared |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| T1 | 1 | 3.994 | 3.176 | 6.843 | 3.143 | 1.271× | 1.010× |
| T1 | 2 | 7.210 | 3.379 | 3.317 | 3.286 | 2.194× | 1.028× |
| T7 | 1 | 227.897 | 165.924 | 174.224 | 179.205 | 1.272× | 0.926× |
| T7 | 2 | 250.282 | 163.537 | 164.290 | 167.062 | 1.498× | 0.979× |

| Frozen target | PDP attempts per arm | Original F6 geometric additions | Attempt-local pair additions / builds / lookups | Shared pair additions / builds / lookups | F6 Boolean reductions in each pair arm | Shared pair peak RSS bytes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| T1 | 1 | 633 | 633 / 0 / 0 | 633 / 0 / 0 | 9 | 6,717,440; 6,553,600 |
| T7 | 11 | 79,635 | 62,295 / 10 / 300 | 6,776 / 1 / 599 | 691 | 8,208,384; 7,979,008 |

Sharing cut T7's counted geometric additions by **91.49%** against
original F6-IC, and by 89.12% against the attempt-local pair variant.
The first failed attempt built the 1,953-pair table after the registered
threshold; subsequent attempts used it immediately. The complete online
wall interval did not improve: both F6/shared ratios are below 1, and
both F4/shared ratios are below the registered 2× gate. The identical F6
A/A ratios across repetitions are 1.064× on T1 and 1.015× on T7; the T1
F4 and attempt-local pair repetitions vary much more. These two
repetitions and this unisolated Mac host cannot establish a statistical
wall-time bound. The operation-count result is deterministic for these
frozen inputs and is **not** a wall-time speedup claim.

**Decision:** The correctness, one-build, T1 operation-count, and T7
≥40%-addition-cut gates pass. The T7 F4/shared ≥2× and F6/shared ≥1.1×
wall gates fail in both repetitions. Stop this candidate at the pilot;
do not expand it to eight targets or route F6-IC to it by default. The
stage evidence indicates that counted geometric additions were a poor
proxy for this complete target's dominant runtime. A useful next
optimization needs a profile of the Boolean search and exact oracle
costs under a controlled host, rather than another pair-table variant.

No same-point rho reference was run here, and no host-level CPU isolation
receipt exists. The one-target IC-versus-rho speedup, normalized total
operation ratio, n83 performance, and 2× F6 claim remain unknown. Peak
RSS is process-wide `getrusage(RUSAGE_SELF)` and is outside the online
timing interval; `/usr/bin/time -l` failed in this sandbox before timing,
as recorded in [HOST_MEMORY_PREFLIGHT.txt](HOST_MEMORY_PREFLIGHT.txt).
The native shared-cache and original pair-closure tests passed in release
mode. Reproduce the 16-row phase and correctness check with `sh derive.sh`;
the result is [DERIVATION_CHECK.json](DERIVATION_CHECK.json).
