#set page(margin: 18mm)

#set text(size: 10pt)

= Prime rootfinding synchronized lane batch inversion

Session `ECBS1hdb2c745b39ed`. Source commit `84373b34c481c3dbd7b9bc067983f7241df247e0`; binary SHA-256 `b24171038803d4873d63dbdf5b5c95d5958b0549a7c40266f559fb3331405105`.

Status `complete`; 256/256 execution records. Host class `ECBENV2hd82681268e96`.

All local wall values are L0 diagnostics; the figure shows counted cold work and table relation yield. Native field operation counts are absent from this ecbench record.

#image("counts.svg", width: 100%)

The bars are arithmetic means over measured executions; whiskers span their observed minimum and maximum. No confidence interval or wall advantage is inferred from the figure.

- `b48-batch-16`: 24/24 verified; mean 7807.292 cold GAE; 8.262035499936152 relations per 1,000 lookups; peak 31216 KiB RSS.

- `b48-batch-4`: 24/24 verified; mean 7814.000 cold GAE; 8.265189483865296 relations per 1,000 lookups; peak 31152 KiB RSS.

- `b48-batch-8`: 24/24 verified; mean 7466.417 cold GAE; 8.411948705806115 relations per 1,000 lookups; peak 31072 KiB RSS.

- `b48-serial-16`: 24/24 verified; mean 7807.292 cold GAE; 8.262035499936152 relations per 1,000 lookups; peak 31264 KiB RSS.

- `b48-serial-4`: 24/24 verified; mean 7814.000 cold GAE; 8.265189483865296 relations per 1,000 lookups; peak 31232 KiB RSS.

- `b48-serial-8`: 24/24 verified; mean 7466.417 cold GAE; 8.411948705806115 relations per 1,000 lookups; peak 31232 KiB RSS.

- `rho-neg`: 24/24 verified; mean 1000.417 cold GAE; N/A relations per 1,000 lookups; peak 27664 KiB RSS.

- `rho-neg-aa`: 24/24 verified; mean 1000.417 cold GAE; N/A relations per 1,000 lookups; peak 27696 KiB RSS.

Raw records: `research/ecbench_rootfinding_lanes_20261013/sessions/lanes-local-20261013/records.jsonl`. Plot data: `summary.json`. Audit receipt is linked from the study README.
