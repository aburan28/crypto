# Amendment 01: streamed second-pass accumulators

Date frozen: 2026-10-06

This amendment is frozen after the initial protocol commit and before source
implementation or execution.

The protocol's candidate-state sentence incorrectly included `16S` bytes for
an accumulator vector.  That vector is unnecessary to the candidate:
second-pass replay visits one start's eight requests contiguously, so it can
construct, test, digest, and emit that start's 16-byte accumulator before
advancing to the next start.

The corrected algorithmic materialisation is:

```text
ceil(U / 64) * 8-byte bitmap words
+ 16U bytes of pair values
+ one 16-byte streaming accumulator
```

At the full 34,562,148,612-entry signed-pair universe this is 557,480,271,372
bytes (about 557.48 GB), excluding allocator metadata and any downstream
output queue.

For complete checked cells, the harness may retain the reference's
accumulator vector so that every candidate output can be compared.  That
reference/harness memory must be recorded separately and must not be presented
as candidate algorithmic state.  Timed candidate runs must stream a digest and
must not allocate an `S`-entry accumulator vector.

All other protocol requirements, including two complete request-generation
passes, exact replay, traffic accounting, scaled-table timing limits, and
promotion gates, remain unchanged.
