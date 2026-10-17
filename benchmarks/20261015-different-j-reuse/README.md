# Different-j full-image setup reuse and process cost

The [source study](../../research/isogeny_torsion_fibers_20261012/different_j_reuse_20261015/RESULTS.md)
holds the frozen protocol and correction, public and checker inputs,
source/binary hashes, raw V2 failure, V3/V4 arm receipts, separate
arithmetic replay, native summaries and figure sources.

This browser copy contains the [two-page run report](run.pdf),
[one-page results report](results.pdf), and the V3
[prefix SVG](prefix.svg)/[150 dpi PNG](prefix.png) and V4
[process SVG](process.svg)/[150 dpi PNG](process.png). Each file is
byte-identical to its source-study counterpart and was generated from
the same native Rust plotting source and Typst renderer.

All 60 distinct order-61 targets were solved by seven arms in three
rounds; V3 and V4 each passed 1,260/1,260 independent scalar replay.
At `N=60`, V4 median whole-child CPU was 15,351 µs full-image,
6,148 µs direct BSGS and 6,908 µs rho. The paired full-image/BSGS
and full-image/rho ratio medians were 2.483 and 2.222. These are
shared-host L0 measurements, with portable counted `S` and isolated
performance replay still open.

The required newest-first row for the shared
[`benchmarks/README.md`](../README.md) is prepared in
[INDEX_ENTRY.md](INDEX_ENTRY.md). Conductor task T-15 currently holds
the index path exclusively, so indexing remains pending; this mirror
does not claim that requirement complete.
