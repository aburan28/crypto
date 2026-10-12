# Wide-S4 CryptoMiniSat correctness controls on disclosed n19 queries

The [exact-yield audit](../generic-exact-yield-audit/RESULTS.md) found two
natural n19a0 queries with real three-point relations that the registered
native-XOR and CNF SAT paths missed under their conflict cap. This
[panel](panel.json) freezes those **already disclosed, outcome-selected**
queries as correctness controls, along with the first two exactly infeasible
queries. They are not a new natural-yield sample. The panel pins the exact
archive and exact-yield result digests, four `[a]G` points, source exporter,
external binary, parser, process meter, order, limits and outcome policy
before any export or solver call. The curve, field representation, geometric
factor base and query points must independently replay from the parent raw
evidence before execution. No sealed or fresh target is used.

**Question.** Can the existing *wide symmetrised S4 circuit* exporter and
CryptoMiniSat's XOR-DIMACS path find a real point decomposition on the two
known-positive queries? The preceding in-process SAT candidate used a
chained-S3 Boolean system with degree-2 Macaulay consequences; this is an
explicit encoding/backend change, not a speed comparison of equivalent
instances. The two negative controls test for invalid point witnesses.

Build `koblitz_pdp_export` from the same pinned Rust commit
`765c3c5f19032bd852163805f257c56babef2040`, with the declared Cargo lock,
offline dependencies and compiler receipt. Verify its file hash and compare
the complete library/dependency source manifest to the parent build before
and after compilation. Each exporter run receives only the registered public
point coordinates and a public nonce, with `--export-only` so it cannot
solve or reveal the known class. Preserve the manifest and all three source
exports, including the exact XOR-DIMACS bytes. Check the exported field,
curve, basis and target against the parent fixture and require the manifest's
source mode to be wide S4. The original geometric factor-base point count
must agree.

Run the pinned local arm64 CryptoMiniSat 5.14.7 executable **once** per fixed
export with one thread, random seed 1, at most one returned model, one million
conflicts and a 120-second process watchdog. `process_meter.py` records
process wall, CPU and peak RSS; the 60-second exporter watchdog is separate.
Preserve stdout/stderr, status, binary/source hashes and exact commands for
every control. A timeout, solver error, missing/partial model or invalid
XOR-DIMACS assignment is inconclusive or invalid, never UNSAT. For a valid
SAT source assignment, decode the first 18 coordinate bits against the
declared six-dimensional basis; independently enumerate rational curve
lifts and group-check all sign choices against the supplied point and the
verified geometric base. A source SAT model without a valid point witness is
reported as nonlifting. A point witness on an exact-negative query or an
UNSAT verdict on an exact-positive query fails scientific consistency.

This is **binary-bound exploratory stage evidence**, because the installed
CryptoMiniSat binary has a retained hash and version but lacks a complete
source/dependency build receipt. Even a valid point witness here does not
make a source-bound, complete one-target SAT IC pipeline. This pilot reports
no SAT natural-yield rate, no paired solver cost, no F5/incumbent/rho speedup,
and no promotion. A positive result would justify integrating a rebuilt,
source-bound SAT backend into the full relation-collection/LA/descent worker
under a new frozen protocol. A negative or timed-out result is retained.
