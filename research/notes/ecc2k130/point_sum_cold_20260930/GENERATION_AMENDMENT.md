# Pre-output fixture-generator repairs

The first generator invocation stopped before calling rho: this sparse
worktree had not materialized the already-committed `PROTOCOL.md`.
After materializing it, the second invocation generated and independently
checked the six streams in memory but stopped at `out.mkdir` because the
protocol directory already existed. Neither invocation wrote a point,
label or `FROZEN.json` file; no timed arm was run or inspected. The
generator now uses `exist_ok=True` for the existing protocol directory
while still refusing an existing `fixtures/` or `FROZEN.json` output.
The committed seed, corpus names, source hashes, cell list, decision rule
and group-law checks are unchanged. The next invocation is the first
one permitted to publish the frozen files.
