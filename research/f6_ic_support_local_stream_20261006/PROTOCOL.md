# Stream F6 support-local Macaulay rows with shared scratch

Registered before implementation or timing. The previous cached-layout
pilots [#1460](../f6_ic_direct_fused_pack_20261006/RESULT.md) and
[#1461](../f6_ic_generic_direct_pack_20261006/RESULT.md) reached zero
F6 layout hits. Source tracing shows F6-IC's inherited engine uses
`build_inherited_macaulay_support_local`, which materializes one temporary
row collection per generator before extending the global rows. This
candidate keeps one pair of product/parity scratch vectors across all
generators and writes surviving rows directly to the global collection.
It also reuses each support-mask/degree-gap multiplier schedule within
the build. Keep the exact monomial sort and odd-multiplicity cancellation,
global row cap, column order, packed matrix, solver decisions, and default
path unchanged. The new path is opt-in.

Use the same prepared n17 `icv1-f2m17-tm101-00378d4e` curve, actual
62-point usable standard base, 29 folded columns, archived public T1/T7
targets, certified logs, three summands, degree three, 8,192 node
budget, 32 target trials, one Rayon thread, and default algorithm
environment. Compare original F6-IC and opt-in support-local streaming
from one release worker. Freeze distinct IC1 candidate IDs, the same
archived workload IDs, source/binary/input hashes and two fresh-process
repetitions per arm and target before timing. Preserve all raw statuses,
failures and timestamps. Both arms must match verified target scalar,
attempts, reductions, geometric additions, matrix rows/columns, word
operations, and the five exclusive online phases must sum to wall time.

Before timing, test exact column and packed-row equality on colliding
Boolean products, zero rows, mixed support masks and at the row cap.
Run existing F6 closure and fused-packing controls. Record target PDP,
F4 build/reduce/readback, complete online time, layout counters and
memory peak when available. The host is unisolated, so CPU ratios remain
exploratory. The public targets are previously seen controls.

Retain only if T7 F6 Macaulay build falls at least 20% and complete
online wall at least 15% in **both** paired repetitions, with no more
than 5% complete-call regression on T1. The requested 2× complete
end-to-end IC target is separate and needs an isolated-host replay plus
ordinary n83 relation yield and a same-target rho comparison.
