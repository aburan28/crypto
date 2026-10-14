# Delivery integration and artifact boundary

The native139-launch campaign ran on Oct12, under source/protocol commit
9da9ff2c and preflight pins adddc355. Its results, all-answer replay and proofs
were committed as0b401aee; the inspected report and durable PDF inputs as6edb9f437.
All measurements and accepted PDFs retain their original hashes.

Delivery resumed on Oct14. Main had advanced fromc62555b7 to275beb93f.
Before integration, a dry merge found only a benchmark-index insertion conflict.
The unpublished branch aligned that metadata index to already-landed main,
merged main cleanly at2ea7fd320, then restored the new dated benchmark row in
chronological order. No shared branch was rebased or force-pushed and no frozen
raw record was selected as a conflict resolution.

Upstream adds diagnostic DP/abandoned-trail vectors to StrongRho, changes public
visibility on one-word Gf2.folds and changes unrelated IC paths. New release
library/example validation is recorded under checks/integrated_20261014.
It passed4089 library tests with117 ignored and all seven worker tests.
That validation concerns the integrated source tree, not a replacement benchmark
or a new runtime measurement. The measured source files are copied byte-for-byte
under archives/measured_sources, with the full tree still pinned in Git.

The Oct12 report is an as-measured record. It is not a survey of measurements
made after its frozen campaign, a new IC/rho ratio, or a fastest-setting claim
over every later implementation. The final comprehensive70-page PDF prepends
the new10 pages to the earlier60, with text and144dpi pixels independently
verified for every historical page. Original accepted SHA-256 values:

- Addendum: ceed3b30b7ee40c68d3854f865f72fa9dabb28b73f73437c421673cbb93796b4
- Comprehensive:856edd55dd4172374f0f60b23e8605230ce7e63ee41d025c2bc7e9e6d0a8fb0c

The new comprehensive PDF and its immutable prior60-page input are committed
alongside the report source, with a separate additional-delivery checksum file.
The original46-entry report and23-entry benchmark-copy manifests remain unchanged.
Study-wide MANIFEST.sha256 is generated after final delivery metadata is staged.
The archive packages existing evidence only; no result is recomputed by packaging.

Owner workflow guidance changed on Oct14 from stopping at a ready PR to following
the authored PR until it lands, using activity/check-in cadence and applicable
CI/review gates. This delivery follows the newer guidance. Earlier frozen protocol
text describes the standing workflow when the experiment was registered.
