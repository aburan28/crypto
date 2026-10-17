# Run: fifth-order missing-price sensitivity

**Start:** 2026-10-16 01:53 UTC, Darwin arm64 shared host. The
[frozen protocol](../../sensitivity_five/PROTOCOL_V1.md),
[source/input/binary hashes](../../sensitivity_five/RUN_SPEC_V1.json),
[exact command and receipt](../../sensitivity_five/RUN.md), and
[native Rust source](../../../../../examples/different_j_price_sensitivity_five.rs)
define this deterministic derivation. It reads completed, independently
replayed solver records; no new DLP solver run or target selection occurs.

The executable exited zero. Its [60 paired rows](../../sensitivity_five/SENSITIVITY_FIVE_V1.json)
cover orders 229, 421, 433, 457 and 1,201, four target prefixes and
three rounds. It checked five run hashes, all 420 child records and
8,925 scalar outputs against their independent replay receipts, and
all recomputed group totals against the native five-order trend.
The first 48 rows and four N64 summaries match the earlier V1 output
byte for byte after canonical JSON serialization.

The plotted unit is **extra omitted rho group-addition equivalents per
observed rho step** needed to tie the full-image counted group lower
bound, assigning zero cost to all omitted full-image work. N64 median
values are 112.973, 303.573, 334.844, 288.595 and **1,390.865** in
ascending order. The new fifth cell's three rounds are 1,266.482–
1,575.602. This threshold is a sensitivity calculation, not a measured
primitive price, a complete portable `S`, or a wall-time result.

![Run 2026-10-16 01:53 UTC: N64 five-order threshold median and three-round range](GRAPH.svg)

The [results report](../../sensitivity_five/RESULTS.md) gives the full
table and interpretation. The host and order variation do not support
a scaling exponent. The complete-price gate still needs disjoint native
operation counts and pinned conversion; the graph adds no canonical
frontier point.
