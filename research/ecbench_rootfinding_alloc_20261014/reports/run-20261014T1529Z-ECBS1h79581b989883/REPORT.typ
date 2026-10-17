#set page(margin: 18mm)

#set text(size: 10pt)

= Prime lane relation collector buffer reuse, baseline versus candidate

Session `ECBS1h79581b989883`. Source commit `431d0b3167027092443705e21d363489b2637a22 (prior postmerge-build-check.json)`; binary SHA-256 `ab5764921ad616ae8177ca57c125ec8cd7ba09896f1e0644a6999d1fb33b7159`.

Status `complete`; 128/128 execution records. Host class `ECBENV2hd82681268e96`.

All local wall values are L0 diagnostics; the figure shows counted cold work and table relation yield. Native field operation counts are absent from this ecbench record.

#image("counts.svg", width: 100%)

The bars are arithmetic means over measured executions; whiskers span their observed minimum and maximum. No confidence interval or wall advantage is inferred from the figure.

- `b48-batch-8`: 24/24 verified; mean 7952.583 cold GAE; 7.739420935412027 relations per 1,000 lookups; peak 31200 KiB RSS.

- `b48-serial-8`: 24/24 verified; mean 7952.583 cold GAE; 7.739420935412027 relations per 1,000 lookups; peak 31360 KiB RSS.

- `rho-neg`: 24/24 verified; mean 1156.792 cold GAE; N/A relations per 1,000 lookups; peak 27744 KiB RSS.

- `rho-neg-aa`: 24/24 verified; mean 1156.792 cold GAE; N/A relations per 1,000 lookups; peak 27760 KiB RSS.

Raw records: `research/ecbench_rootfinding_alloc_20261014/sessions/baseline-a0/records.jsonl`. Plot data: `summary.json`. Audit receipt is linked from the study README.
