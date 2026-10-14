#set page(margin: 18mm)

#set text(size: 10pt)

= Prime rootfinding fixed-stride coverage holdout

Session `ECBS1h1b585c89ebd1`. Source commit `e443d306a1b79a157ef10a96a4eed72e9496814b`; binary SHA-256 `b24171038803d4873d63dbdf5b5c95d5958b0549a7c40266f559fb3331405105`.

Status `complete`; 160/160 execution records. Host class `ECBENV2hd82681268e96`.

All local wall values are L0 diagnostics; the figure shows counted cold work and table relation yield. Native field operation counts are absent from this ecbench record.

#image("counts.svg", width: 100%)

The bars are arithmetic means over measured executions; whiskers span their observed minimum and maximum. No confidence interval or wall advantage is inferred from the figure.

- `b48-stride-1`: 24/24 verified; mean 8097.833 cold GAE; 7.84288166156119 relations per 1,000 lookups; peak 31184 KiB RSS.

- `b48-stride-23`: 24/24 verified; mean 7802.125 cold GAE; 8.29856335456132 relations per 1,000 lookups; peak 31200 KiB RSS.

- `b48-stride-7`: 24/24 verified; mean 8238.417 cold GAE; 7.865774366456279 relations per 1,000 lookups; peak 31152 KiB RSS.

- `rho-neg`: 24/24 verified; mean 886.333 cold GAE; N/A relations per 1,000 lookups; peak 27664 KiB RSS.

- `rho-neg-aa`: 24/24 verified; mean 886.333 cold GAE; N/A relations per 1,000 lookups; peak 27648 KiB RSS.

Raw records: `research/ecbench_rootfinding_lanes_20261013/sessions/coverage-holdout-local-20261013/records.jsonl`. Plot data: `summary.json`. Audit receipt is linked from the study README.
