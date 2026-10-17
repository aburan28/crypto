#set page(margin: 18mm)

#set text(size: 10pt)

= Prime lane relation collector buffer reuse, baseline versus candidate

Session `ECBS1hd3d1c462e777`. Source commit `316cc8386611d9b6ba7a57c86442b4d59ba9ee58`; binary SHA-256 `db9cf52602e8b3bc19a885c3eae99c899bd6acc541958246f3e540e58e8de743`.

Status `complete`; 128/128 execution records. Host class `ECBENV2hd82681268e96`.

All local wall values are L0 diagnostics; the figure shows counted cold work and table relation yield. Native field operation counts are absent from this ecbench record.

#image("counts.svg", width: 100%)

The bars are arithmetic means over measured executions; whiskers span their observed minimum and maximum. No confidence interval or wall advantage is inferred from the figure.

- `b48-batch-8`: 24/24 verified; mean 7952.583 cold GAE; 7.739420935412027 relations per 1,000 lookups; peak 31376 KiB RSS.

- `b48-serial-8`: 24/24 verified; mean 7952.583 cold GAE; 7.739420935412027 relations per 1,000 lookups; peak 31424 KiB RSS.

- `rho-neg`: 24/24 verified; mean 1156.792 cold GAE; N/A relations per 1,000 lookups; peak 27680 KiB RSS.

- `rho-neg-aa`: 24/24 verified; mean 1156.792 cold GAE; N/A relations per 1,000 lookups; peak 27776 KiB RSS.

Raw records: `research/ecbench_rootfinding_alloc_20261014/sessions/candidate-b/records.jsonl`. Plot data: `summary.json`. Audit receipt is linked from the study README.
