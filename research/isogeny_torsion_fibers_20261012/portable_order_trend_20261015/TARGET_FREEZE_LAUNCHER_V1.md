# Target-freeze launcher failure retained

The first three-cell shell launcher ran the `p=409,r=433` freeze child to
completion and wrote its public/checker JSON and complete stdout, but then
exited 1 on `zsh:8: read-only variable: status`. `status` is a reserved zsh
parameter. That launcher had no child-exit receipt for the first cell; its
child stdout reports the source/input hashes, and independent `shasum -a 256`
of both files agrees exactly. The first cell was not rerun or overwritten.

The corrected thin launcher used `run_status` and froze only the two
remaining cells. Both child exit receipts are zero and their stdout hashes
match the files. The shell failure is an orchestration error, not a curve or
factor-base rejection. All three target files are frozen in
`BATCH_RUN_SPEC_V1.json` before any measured child is launched.
