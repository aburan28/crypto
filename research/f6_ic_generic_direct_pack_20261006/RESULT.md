# Generic cached-layout direct packing did not reach F6-IC

**Decision: reject.** The corrected frozen comparison completed 8/8
one-target solves and independently replayed every scalar, but the F6-IC
inherited engine had **zero cached-layout hits and zero misses** in both
arms on both targets. The opt-in direct branch in the generic builder was
never entered. The predeclared requirement for at least one T7 F6 hit
failed, as did the consistent build and complete-call gain gates. The
default remains unchanged.

The [protocol](PROTOCOL.md) was committed as `e6f497177` before code or
timing. The corrected worker implementation was committed as `5e4549c43`
and frozen at `ff56ecbe4` before the successful eight-run comparison.
The worker binary SHA-256 is in [BUILD_IDENTITY.txt](BUILD_IDENTITY.txt),
with exact [IC1 manifests](candidates/), [input and candidate hashes](FREEZE.tsv),
[source hashes](SOURCE_SHA256SUMS), and [raw runs](runs/). Two earlier
frozen runs failed the worker's declared-environment preflight before a
target solve. Their raw rows, old frozen identities, causes, and checksums
are retained in [PREFLIGHT_FAILURE.md](PREFLIGHT_FAILURE.md) and its
adjacent directories. Neither failure is counted as a successful or
timed one-target run.

The corrected workload fixes the same prepared n17 Koblitz curve,
62 actual usable factor-base points, 29 folded columns, T1/T7 public
targets, one Rayon thread, three summands, degree three, and imported
certified logs. Both arms explicitly declare and set
`KIC_F4_ACTIVE_MULTIPLIERS=1`; only the direct-packing flag differs.
The online interval begins after reusable setup and includes query,
every target PDP attempt, relation check, descent, and scalar replay.
All five exclusive online phases summed exactly to each reported wall
time. No paired rho solve was run in this stage study.

| Target | Rep | Sorted online ms | Direct online ms | Sorted/direct | Sorted build ms | Direct build ms | Layout hits, both |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| T1 | 1 | 3.639 | 3.575 | 1.018× | 2.763 | 2.973 | 0 |
| T1 | 2 | 6.507 | 6.979 | 0.932× | 4.698 | 4.287 | 0 |
| T7 | 1 | 272.039 | 181.706 | 1.497× | 183.501 | 136.308 | 0 |
| T7 | 2 | 188.530 | 188.570 | 1.000× | 141.706 | 144.443 | 0 |

The T7 first-repetition movement cannot be credited to direct packing:
the branch did not execute, and the second repetition was equal within
0.03%. All sorted/direct pairs had identical attempts, reductions,
geometric additions, matrix rows/columns, word-XOR counts, and layout
counters. The first direct parity unit, new generic-builder cache and
stale-layout unit, existing fused-packing equivalence, F6 geometric
closure, and the worker environment-declaration unit passed. Release
build logs and source/binary hashes are retained. The host was an
unisolated macOS host whose physical CPU model was unavailable through
the sandbox; CPU ratios are exploratory and memory peak is unknown.

The root cause is architectural. F6-IC selects `SolverEngine::InheritedF4`.
That engine enables support-local multipliers by default, and its
`ReducedBasis::from_system_closed_with` path calls
`build_inherited_macaulay_support_local`. This builds rows separately
for each generator and then computes columns and packs them. It bypasses
both the inherited cached-layout and generic cached-layout paths tested
here. `KIC_F4_ACTIVE_MULTIPLIERS` changes a different decisive
`matrix_f4_f2_counted_impl` path, so setting it did not reroute this F6
workload. The next F6 matrix-build change must act on support-local row
generation and preserve its exact row space and final one-target result.

Run `sh research/f6_ic_generic_direct_pack_20261006/derive.sh` to rebuild
[measurement rows](measurements.jsonl), [paired correctness check](DERIVATION_CHECK.json),
and [derived hashes](DERIVED_SHA256SUMS) from the raw successful runs.
There is no n83 ordinary-relation, end-to-end IC-versus-rho, or 2×
speedup claim from this pilot.
