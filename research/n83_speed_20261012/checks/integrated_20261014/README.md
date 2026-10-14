# Integrated validation, 14 October 2026

The clean integrated revision `2ea7fd32086dfaf4189b1f95071739ae3f1d8ffe`
passes the required release checks:

| Command | Passed | Failed | Ignored |
| --- | ---: | ---: | ---: |
| `cargo test --release --locked --lib` | 4,089 | 0 | 117 |
| `cargo test --release --locked --example n83_rho_speed_eval` | 7 | 0 | 0 |

Both Cargo commands exited zero. The new worker source remains SHA-256
`bd2e5b2ed4ca4229f187b3f1c5c613633d9cb773200c961fdd8cbe8fb1d65241`.
The integrated library includes upstream diagnostic and unrelated changes;
these checks establish correctness regression coverage for that revision.
They are **not** new performance measurements. The frozen October 12 measured
worker, sources, raw records and final PDFs remain unchanged, and the historical
results retain their original executable and source bindings.

`commands.txt` records the exact compiler/cache settings and commands.
`receipt.json` binds the revision, source and log hashes, result counts and
frozen artifact hashes. The source hash recheck passes for all seven bound
files. Library tests ran with network permission for the repository's local
socket checks; no credential values or process-environment dump was recorded.

Compilation and all tests succeeded on their first integrated attempt. A
post-run console summary command could not find `rg` under the minimal build
PATH. Its exit-127 reproduction is retained separately. Extracting the final
result lines under the ordinary tool PATH succeeded; neither Cargo result
was affected, and no tests were skipped or rerun to alter the result.

Conductor was attempted before writing this check directory and remains
unavailable. No source code, benchmark records or frozen reports were edited.
