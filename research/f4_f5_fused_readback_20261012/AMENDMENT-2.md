# Second candidate, before its measurement

The first candidate `07151ff0` passed exact-output/accounting and memory
controls, but **failed the frozen two-arm instruction acceptance gate**.
The native [first comparison](records/comparison.json) records F4 as neutral
(1.4363% reduction, below 2%), F5 as fewer instructions (11.1305%), and zero
memory regressions. F4 allocation requests decreased 18.4494%. Keep all
first-attempt profiles, controls, build logs and source/binary hashes unchanged.
Do not relabel its instruction verdict as accepted.

The fresh baseline profile records 20,596,539 F4 self instructions in column
lookup over 33 complete calls. The next hypothesis is a compact direct index
for small observed monomial-mask ranges. Represent present columns as `u16`
indices and missing masks by `u16::MAX`; use it only when the table's byte
count is no larger than the minimum `(u64, usize)` entry storage of the current
hash map. Reject index overflow and too-sparse/high-bit ranges before allocating.
Preserve the current mixed-mask hash map fallback and the standard-hash control.
`KIC_F4_HASH_COLUMN_INDEX=1` forces the former mixed-mask index for diagnostics.

This changes lookup representation only. Preserve column lists, row order,
packed matrices, caps, pruning and XOR accounting. Retain the first candidate's
ordered F4 readback and fused F5 builder. Add independent sparse-hole, bit-63,
overflow and sentinel-boundary lookup checks against a standard map; keep
existing exact-layout, stale-cache, flat/nested and solver checks.

All eleven inputs, seeds, one-worker coverage, three primary repetitions,
five rotated A/B/A′ pairs, bootstrap parameters and **both complete-call
instruction/no-memory-regression gates remain unchanged**. Do not expand or
replace workloads to earn acceptance. Stop after this second candidate; if its
gate fails, retain the rejected result and deliver only independently accepted
changes, explicitly separating any allocation-only improvement from an
instruction improvement.

Save final-attempt data under `records/attempt-2/`. Reuse the unchanged original
baseline binary and its primary/A/A profiles (copy byte for byte); run fresh
candidate primary and all five paired rounds. Retain the first candidate binary
in scratch as well as its committed source revision and hashes. The final
compatibility control sets **both** `KIC_F5_FUSED_BUILD=0` and
`KIC_F4_HASH_COLUMN_INDEX=1` in a fresh process, then compares all outputs and
accounting to the baseline. Retain its explicit environment file. Each attempt
contains the fixed 792-call core dataset and a 33-call compatibility dataset;
the two baseline processes are shared, not additional independent samples.

The final report/PDF includes the rejected first attempt. Native instruction
history may plot the previous hashing pair and the final pair as disconnected
segments, with the native cross-study output/accounting check cited. No ratio
is compounded across intervening unrelated merges. Old study plots remain
frozen. Required library/example checks and publication rules are unchanged.
