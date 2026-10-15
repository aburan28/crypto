#set page(margin: 18mm)

#set text(size: 10pt)

= Prime rootfinding B retune train with batch-8

Session `ECBS1hbdfb5bc7d8a8`. Source commit `193f50982720da9608d295f2e090020f063d12fd`; binary SHA-256 `b24171038803d4873d63dbdf5b5c95d5958b0549a7c40266f559fb3331405105`.

Status `complete`; 256/256 execution records. Host class `ECBENV2hd82681268e96`.

All local wall values are L0 diagnostics; the figure shows counted cold work and table relation yield. Native field operation counts are absent from this ecbench record.

#image("counts.svg", width: 100%)

The bars are arithmetic means over measured executions; whiskers span their observed minimum and maximum. No confidence interval or wall advantage is inferred from the figure.

- `b40-batch-8`: 24/24 verified; mean 8225.208 cold GAE; 5.645351902692679 relations per 1,000 lookups; peak 31248 KiB RSS.

- `b44-batch-8`: 24/24 verified; mean 8262.875 cold GAE; 6.841947144504637 relations per 1,000 lookups; peak 31232 KiB RSS.

- `b48-batch-8`: 24/24 verified; mean 7957.417 cold GAE; 8.018646978787947 relations per 1,000 lookups; peak 31184 KiB RSS.

- `b52-batch-8`: 24/24 verified; mean 8192.167 cold GAE; 9.172080742104715 relations per 1,000 lookups; peak 31088 KiB RSS.

- `b56-batch-8`: 24/24 verified; mean 8289.333 cold GAE; 10.63267565935172 relations per 1,000 lookups; peak 31200 KiB RSS.

- `b64-batch-8`: 24/24 verified; mean 8706.292 cold GAE; 14.132731113863189 relations per 1,000 lookups; peak 31296 KiB RSS.

- `rho-neg`: 24/24 verified; mean 1128.958 cold GAE; N/A relations per 1,000 lookups; peak 27728 KiB RSS.

- `rho-neg-aa`: 24/24 verified; mean 1128.958 cold GAE; N/A relations per 1,000 lookups; peak 27680 KiB RSS.

Raw records: `research/ecbench_rootfinding_lanes_20261013/sessions/retune-train-local-20261013/records.jsonl`. Plot data: `summary.json`. Audit receipt is linked from the study README.
