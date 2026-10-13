This follows merged PR #1606 and extends the standalone `isogeny` binary for the requested campaign of **at least 2000 additional P192/P224 maps**. It raises the walk request bound and adds native accounting, exact replay controls, source/binary archives and durable recovery evidence. The complete campaign result and canonical model registration remain open.

Two fresh 1050-step walks use P192 degree 2341 and P224 degree 1117. The first 50 maps per class reproduce the preceding campaign and are excluded from additional yield. The native accounting utility compares source ICV1, prime degree and canonical monic kernel coefficients, so differing receipt metadata cannot count an old map again. A later archival check rejects alternate coefficient encodings and unreduced values.

The CLI now admits 1–10000 maps, using the same bound for construction and full-chain replay. Every step still requires exact independent kernel, codomain and rational curve-equation checks, subgroup order and 20 public scalar-transport checks before advancing. Coverage stays PARTIAL at one Frobenius eigenline per step. Map count and distinct degree count are separate; these are consecutive edges at two prime degrees.

The running CLI was built clean at `ecc2b2250`; the native supervisor was frozen at `316150cba`. Both binary hashes, protocol hash and preceding manifest hash are recorded in the launch receipt. The supervisor retains child outputs/exit receipts and writes a final unique-count summary automatically when both walks end. Raw ongoing output uses the system disk because the research SSD has limited headroom. The initial committed report/PDF records 25 repeated steps and **0 additional maps**; live counts are separate.

Research records: [protocol](research/p192_p224_2000_more_large_prime_isogenies_20261010/PROTOCOL.md), [running state](research/p192_p224_2000_more_large_prime_isogenies_20261010/RUNNING.md), [native commands](research/p192_p224_2000_more_large_prime_isogenies_20261010/COMMANDS.md), [initial report](research/p192_p224_2000_more_large_prime_isogenies_20261010/snapshots/initial/RESULTS.md), [artifact inventory](research/p192_p224_2000_more_large_prime_isogenies_20261010/FILES.json), and [previous hundred-map closeout](research/p192_p224_hundred_large_prime_isogenies_20261010/FOLLOWUP_VALIDATION.md).

Validation for this extension:

- CLI: 63 tests plus 3 independent field-reference tests pass.
- Native running supervisor: 14 tests pass. Later archival utility: 16 tests and strict Clippy pass. Changed source passes rustfmt 1.98.
- A saved two-step chain exactly replays under a 1050-step request while reporting `complete: false`.
- Modified source/model and mapped generator with recomputed hashes, and an empty chain, are rejected. Positive replay passes.
- The archived standalone binary runs outside the checkout and passes complete two-map P384/19 construction and independent replay; these control maps do not count toward the campaign.
- The required root release-library gate passed 3815 tests, with 0 failures and 110 existing ignored tests, in 2847.54 seconds. Its durable receipt is retained in the October 12 recovery evidence.

The earlier 100-map raw evidence, report and canonical-at-run audit remain frozen. The 449-model inventory covers the pinned standards sources with explicit operation capabilities and unresolved adapters; the IC performance ledger is unchanged. Timings are operational shared-host measurements, with constructor/replay/invocation costs separate and repeated prefix costs retained.

The durable October 11 run stopped on storage exhaustion after sealing 153 additional maps: 95 P192 and 58 P224. These certificates survive. The October 12 continuation uses recovered SSD capacity and fresh output, with union-based deduplication across attempts; repeated prefix work does not inflate yield. See RECOVERY_20261012.md. The original scratch-loss observation and initial zero-yield report remain historical records.

PR #1606 is merged. The latest published release checked has no isogeny-named asset, so the committed standalone archive provides the verified local runtime. Current-head CI and subsequent ReleaseMe publication remain separate gates. Final new-map archival, canonical registry/graph updates and the 2000-map acceptance check remain outstanding.
