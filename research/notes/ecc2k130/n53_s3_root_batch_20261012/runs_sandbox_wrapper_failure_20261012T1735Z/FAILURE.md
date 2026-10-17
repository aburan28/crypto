# First execution attempt: measurement-wrapper failure

UTC start: 2026-10-12 17:35. This attempt used the frozen source, inputs,
schedule, and binaries in `../FROZEN.json`, but ran inside the filesystem
sandbox. macOS `/usr/bin/time -l` returned exit code 1 and wrote
`time: sysctl kern.clockrate: Operation not permitted` after the child
process exited. The wrapper therefore supplied no maximum-resident-set-size
line, and the native runner classified every completed cell as `INCOMPLETE`.
The run was stopped after eight completed cells and while the ninth was
starting. Its raw intents, outputs, stderr, and status records are retained
here. These cells are excluded from the frozen comparison; no timing or
correctness verdict is inferred from them. A one-command wrapper check outside
the sandbox returned exit code 0 with maximum RSS, so the same frozen binary
will be rerun there from the beginning. The new run writes to `../runs/`.
