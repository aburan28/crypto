#set page(margin: 18mm)

#set text(size: 10pt)

= Prime rootfinding B retune holdout with batch-8

Session `ECBS1hadb0689b12a0`. Source commit `0bd3f1e53ad6e214bb74376f29a2499c128a0be4`; binary SHA-256 `b24171038803d4873d63dbdf5b5c95d5958b0549a7c40266f559fb3331405105`.

Status `complete`; 256/256 execution records. Host class `ECBENV2hd82681268e96`.

All local wall values are L0 diagnostics; the figure shows counted cold work and table relation yield. Native field operation counts are absent from this ecbench record.

#image("counts.svg", width: 100%)

The bars are arithmetic means over measured executions; whiskers span their observed minimum and maximum. No confidence interval or wall advantage is inferred from the figure.

- `b40-batch-8`: 24/24 verified; mean 8356.750 cold GAE; 5.664739884393064 relations per 1,000 lookups; peak 31232 KiB RSS.

- `b44-batch-8`: 24/24 verified; mean 8290.208 cold GAE; 6.839522392673596 relations per 1,000 lookups; peak 31152 KiB RSS.

- `b48-batch-8`: 24/24 verified; mean 8328.292 cold GAE; 7.973655869439152 relations per 1,000 lookups; peak 31232 KiB RSS.

- `b52-batch-8`: 24/24 verified; mean 8283.750 cold GAE; 9.515432592355669 relations per 1,000 lookups; peak 31216 KiB RSS.

- `b56-batch-8`: 24/24 verified; mean 8151.708 cold GAE; 11.355478692162109 relations per 1,000 lookups; peak 31200 KiB RSS.

- `b64-batch-8`: 24/24 verified; mean 8764.125 cold GAE; 14.889122030439344 relations per 1,000 lookups; peak 31392 KiB RSS.

- `rho-neg`: 24/24 verified; mean 1046.708 cold GAE; N/A relations per 1,000 lookups; peak 27680 KiB RSS.

- `rho-neg-aa`: 24/24 verified; mean 1046.708 cold GAE; N/A relations per 1,000 lookups; peak 27696 KiB RSS.

Raw records: `research/ecbench_rootfinding_lanes_20261013/sessions/retune-holdout-local-20261013/records.jsonl`. Plot data: `summary.json`. Audit receipt is linked from the study README.
