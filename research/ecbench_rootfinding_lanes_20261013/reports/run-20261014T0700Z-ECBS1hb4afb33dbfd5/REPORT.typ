#set page(margin: 18mm)

#set text(size: 10pt)

= Prime rootfinding fixed-stride coverage train

Session `ECBS1hb4afb33dbfd5`. Source commit `9036e0654613eb48296005b244f4790d7e3027da`; binary SHA-256 `b24171038803d4873d63dbdf5b5c95d5958b0549a7c40266f559fb3331405105`.

Status `complete`; 160/160 execution records. Host class `ECBENV2hd82681268e96`.

All local wall values are L0 diagnostics; the figure shows counted cold work and table relation yield. Native field operation counts are absent from this ecbench record.

#image("counts.svg", width: 100%)

The bars are arithmetic means over measured executions; whiskers span their observed minimum and maximum. No confidence interval or wall advantage is inferred from the figure.

- `b48-stride-1`: 24/24 verified; mean 7734.042 cold GAE; 8.516201985935899 relations per 1,000 lookups; peak 31120 KiB RSS.

- `b48-stride-23`: 24/24 verified; mean 7584.375 cold GAE; 8.179858934169278 relations per 1,000 lookups; peak 31120 KiB RSS.

- `b48-stride-7`: 24/24 verified; mean 8408.667 cold GAE; 8.029034285598426 relations per 1,000 lookups; peak 31184 KiB RSS.

- `rho-neg`: 24/24 verified; mean 1131.333 cold GAE; N/A relations per 1,000 lookups; peak 27728 KiB RSS.

- `rho-neg-aa`: 24/24 verified; mean 1131.333 cold GAE; N/A relations per 1,000 lookups; peak 27696 KiB RSS.

Raw records: `research/ecbench_rootfinding_lanes_20261013/sessions/coverage-train-local-20261013/records.jsonl`. Plot data: `summary.json`. Audit receipt is linked from the study README.
