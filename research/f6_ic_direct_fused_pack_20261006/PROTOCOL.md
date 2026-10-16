# F6-IC direct parity packing into a cached Macaulay layout

Registered before implementation or timing. The frozen n17 F6-IC T7
profile in PR #1459 attributed 83.1–83.9% of target PDP to Macaulay
work, and 98.5% of that to matrix construction. Its geometric oracle
took only 1.6%. A previous Boolean F4 experiment that replaced each
sorted product with a fresh parity hash set regressed complete calls
([receipt](../boolean_f4_unsorted_products_20260926/RESULT.md)). This
candidate uses a different exact cancellation mechanism: when a cached
column layout is present, XOR each product term directly into its packed
row bit, so repeated monomials cancel there. It allocates no per-product
hash set and does no product sort in that cached-layout path. Only bits
that remain set after cancellation count as observed columns. The
ordinary sorted path handles cache misses, and the default remains it.

Freeze the same prepared n17 `icv1-f2m17-tm101-00378d4e` curve, actual
62-point usable standard base, 29 folded columns, dimension 6, certified
logs, T1 and T7 public targets, seed `20261004039`, three summands,
degree-three inherited F4, 8,192 node budget, 32 target trials, and one
Rayon thread. Compare four arms from one new binary: inherited F4 sorted,
inherited F4 direct, original F6-IC sorted, and F6-IC direct. Freeze
distinct IC1 candidate IDs for each algorithm-affecting option and the
same archived workload IDs. Run two fresh-process repetitions per arm
and target in alternating order. Preserve every status, raw output,
stderr, source/binary/input hashes, and five exclusive online phases.
This is a bounded remeasurement on seen public points, not an independent
speedup validation.

Before timing, add native tests that construct products with colliding
Boolean monomials and prove byte-for-byte equality of sorted and direct
packed rows, including zero products, missing cached columns, and the
exact-support condition. Run the existing Boolean F4 and F6 closure
correctness controls. Keep all solver decisions, scalar replay,
attempts, reductions, matrix dimensions, and word-XOR counts equal
between each sorted/direct pair. Record cached-layout hit and miss
counts, F4 build/reduce/readback timers, online PDP, complete online wall,
and peak RSS when available. The new direct mode must not change the
generic flat F4/F5 matrix builder.

Retention requires exact verified outputs, at least 20% less F6 Macaulay
build time on T7 in both repetitions, at least 15% lower complete F6
online wall time on T7 in both repetitions, and no more than 5% higher
T1 F6 online time. The user's broader **2× complete F4/F6 target** is
an additional gate, not inferred from a build-only result. If either
retention gate fails, preserve the negative result and keep direct
packing opt-in. The host is an unisolated Mac; all CPU timing ratios
remain exploratory until an auditable isolated-host replay. This n17
three-summand result cannot establish n83 ordinary relation yield or
end-to-end IC-versus-rho speedup.
