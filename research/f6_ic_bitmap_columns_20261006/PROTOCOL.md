# Exact ranked bitmap for F6 support-local Macaulay columns

Registered before implementation or timing. The F6 inherited engine
builds support-local Boolean Macaulay rows. Its current column collector
flattens all row-term occurrences, sorts them, deduplicates and then sorts
the unique masks into degree reverse lexicographic order. This candidate
maps each degree-at-most-four monomial mask to its combinatorial rank,
marks a compact bitmap, retains a mask only on its first observation,
and sorts only the unique masks into the unchanged column order. It sees
only already parity-cancelled row terms, so colliding Boolean products
must remain cancelled before they enter the bitmap. Unsupported widths
or degrees fall back to the original collector. The path is opt-in and
the default remains unchanged.

Use the frozen prepared n17 `icv1-f2m17-tm101-00378d4e` curve, 62 actual
usable standard-base points, 29 folded columns, archived T1/T7 public
targets, imported certified logs, three summands, degree three, 8,192
node budget, 32 target trials, one Rayon thread, and default algorithm
environment. Compare original F6-IC against F6-IC with bitmap columns,
one release binary and distinct IC1 candidate IDs. Freeze binary,
source, input and candidate hashes before timing. Run two fresh-process
repetitions per arm and target in alternating order. Preserve all raw
statuses, failures, timestamps and five exclusive online phases.

Before timing, verify exact columns and packed matrices on mixed masks,
duplicates, zero rows, degree bounds and a deterministic family of
small systems; run F6 closure and existing F4 packing controls. In every
paired target solve require the same verified scalar, attempts,
reductions, geometric additions, matrix dimensions and word-XOR counts.
Record target PDP, matrix build/reduce/readback, complete online wall,
layout counters and memory peak where available. The host is unisolated,
so wall ratios are exploratory. Seen planted/public controls do not
establish ordinary n83 yield.

Retain only if T7 matrix build falls at least 20% and complete online
time at least 15% in **both** repetitions, with no more than 5% T1
regression. A controlled 2× complete F6 or IC-versus-rho claim requires
an isolated replay and verified end-to-end single-target evidence.
