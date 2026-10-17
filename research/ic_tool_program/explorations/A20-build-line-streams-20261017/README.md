# A20: the folded build's streams written a line at a time (smoke only)

**What this is.** This is a candidate, its tests and a smoke check, on 2026-10-17. It is not a
timed exploration: the smoke showed nothing worth one.

**Why.** [A19](../A19-filter-pages-build-parts-20261017/README.md)'s build
timers put the append into the partition streams at 12–13 ns of a stored
pair's 31–35 ns. The usual remedy is software write-combining:
- **the pending line:** each stream keeps a 64-byte line in L1;
- **the stores:** a full line goes out with non-temporal stores, so the
  stream's lines are neither read before they are written nor kept in the
  cache.

**The outcome: no visible gain, so it is not pursued.** In four isolated processes per size, A20's
build reads 119 and 114 ms against A19's 128 and 108 at 2^44.3, and 238 and
220 against 231 and 229 at 2^47.2 ([`smoke.tsv`](smoke.tsv)). Both rows'
logarithms are unchanged.
- **Probable cause:** the append is not bound by line traffic. The likelier
  cost is the first touch of the streams' fresh memory: about 12 bytes a
  stored pair, 45–75 MB a build, or 11–19k page faults. Non-temporal stores
  do not avoid those faults. This is not measured.
- **What it rules out:** the streams' write traffic, as a lever on the build.

## The candidate

[`line-streams-on-a19.patch`](line-streams-on-a19.patch) is one commit (`21564f85`) on
[A19](../A19-filter-pages-build-parts-20261017/)'s candidate. It adds `LineStream`:
- **The line:** a pending 64-byte line, flushed with four `_mm_stream_si128` stores.
- **The storage:** an over-reserved `Vec` whose elements start on a 64-byte boundary.
- **The fallback:** ordinary appends if a stream outgrows its room.
- **The table:** unchanged, since push order is kept.

The tests passed: `cargo test --release --locked --lib -- koblitz_index_calculus
koblitz_fast`, 207 passed, none failed, 3 ignored. They include the checks
that the partitioned and scalar builds give the same table.

The binaries are hashed in [`binaries.sha256`](binaries.sha256).
