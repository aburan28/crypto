# F6-IC direct parity packing at the generic cached-layout builder

Registered before implementation or timing. The preceding [direct-packing
pilot](../f6_ic_direct_fused_pack_20261006/RESULT.md) was exact but reached
zero F6 layout hits on the frozen n17 workloads: its opt-in path was wired
only into the inherited F4 builder. This pilot wires that same parity
mechanism into the generic F4/F6 cached-layout builder, where F6 can use
it. It remains opt-in and leaves the generic uncached and F5-criterion
paths unchanged. A repeated Boolean product monomial must cancel before
row/column support is counted.

Freeze the prepared n17 `icv1-f2m17-tm101-00378d4e` curve, actual 62-point
usable standard base, 29 folded columns, archived public targets T1 and
T7, three summands, degree three, 8,192 node budget, 32 target trials,
one Rayon thread, and imported certified logs. Force
`KIC_F4_ACTIVE_MULTIPLIERS=1` in **both** arms so this small-system
control exercises the wider-system layout-reuse policy. The arms are
original F6-IC with the sorted generic builder and F6-IC with opt-in
direct parity packing. Use one new release worker binary, distinct IC1
candidate IDs that record both algorithmic flags, same archived workload
IDs, and two fresh-process repetitions per arm and target in alternating
order. Source, binary, input, candidate, environment, and run hashes
must be frozen before timing. Preserve every raw failure, timeout, or
successful run, along with five exclusive online phases and scalar
replay. These are seen public points and an unisolated Mac; wall ratios
remain exploratory.

Before timing, add a native test that the generic sorted/direct cached
matrix builders produce the same exact columns and packed rows, including
colliding products, zero rows, a valid cache hit, and stale cached layout
fallback. Run the existing F6 closure and fused-packing controls. In the
paired target runs, require same verified scalar, attempts, reductions,
geometric additions, matrix rows/columns, word-XOR counts, and layout
hit/miss counts. Record online wall, target PDP, F4 build/reduce/readback,
and layout counters.

Reject if F6 has zero cached hits on T7, if exactness differs, if direct
packing fails to lower F6 Macaulay build by at least 20% in **both** T7
repetitions, or if it fails to lower complete F6 online time by at least
15% in **both** T7 repetitions. Require no more than 5% regression on
either T1 run. The 2× complete-call target is a further goal, never
inferred from a build-only result. A retained candidate still needs an
isolated-host replay and n83/full-IC evidence before any controlled or
end-to-end speedup claim.
