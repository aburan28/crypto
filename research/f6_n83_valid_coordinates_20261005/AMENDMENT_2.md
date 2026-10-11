# Paired complete planted-block search

Registered after the root pairs and before the following measurement.
The root coordinate gate passed with identical outcomes, and the
candidate 128-node planted probe returned `[0,2,4]` with exact source
and subgroup-projection replay. Pair the previous lazy-child solver
against the indexed-coordinate solver on that **complete planted n83
three-summand block** to check that the root saving survives a real
witness search.

Use the frozen previous witness binary SHA-256
`8d79eebbb93634a7a4faa3f58ef3cb1827797d079abe22ef35aedb1dd3cc7cd2`
and current candidate witness binary SHA-256
`55a25c8ccf25e5c9184990160a586ebe41ec3e1f5b835f5dd5715ca097f23eeb`.
Keep the exact K0 curve, source base, planted indices `[0,2,4]`,
`m=3`, node budget 128, `RAYON_NUM_THREADS=1`, and release profile.
Run five alternating baseline/candidate pairs with a 120-second cap
per process. Preserve all JSON, status and stderr files. Require both
arms to find and replay a witness on every pair, with matching search
counters and matrix dimensions. A local diagnostic gain requires at
least 2× lower median query interval and no more than 1.25× median
peak RSS. Preserve a failed gate or regression; do not change the
runtime code after freezing the binaries without a new comparison.

This is still a planted correctness control on an unisolated Mac. It
does not estimate ordinary-query yield, eight-summand cost, an IC
candidate's one-target online time, or a same-point rho ratio.
