# F6-IC two-word geometric closure: n=83 arithmetic gate

Frozen before the first timed run on 2026-10-04. This is an arithmetic and
correctness diagnostic, not an IC candidate or an end-to-end DLP comparison.

Use the public K_1 curve and subgroup generator from
`research/ic_tool_program/conformance/v2-b3b/params/C082-kic-two-word-n83.json`.
This is the 53-bit subgroup on the confidence gate's field, **not** the
81-bit K_0 confidence-gate subgroup. Generate the first 64 nonzero multiples
of that generator outside timing. For each prefix size 16, 32 and 64, form
all unordered pair sums including repeated points in lexicographic pair
order. Compare the reference `point_add` loop to `batch_add_fixed`, byte for
byte. Time each method three times, alternating reference/batch order by
round, with one warmup for each, using the release build. Retain raw
durations and every failure.

Build one capped pair index for each prefix, check a planted four-summand
target and a target outside the reachable scalar interval. These targets
are correctness controls; they do not estimate ordinary relation yield.
Record the index's pair count. The primary observation is the paired
pair-formation time, including output allocation. Also time a full
no-witness target query through the same index using scalar residual
additions versus the batched residual path, with the same three-round
alternating order and one warmup. No result here can be
reported as an F6/F4/F5 or IC/rho speedup. CPU ratios on this host remain
exploratory unless the repository's isolation receipt gate passes.
