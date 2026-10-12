# Native isogeny setup and receipt-hashing continuation

Dated 2026-10-09. Continue the owner's request for an optimized dedicated binary,
comprehensive curve coverage, a PR, and ReleaseMe publication. Prior evidence remains
frozen; this round retains every coverage and release gate still outstanding.

## Hypothesis and reference

The reference is commit `5cbc3a2dd61028aa5fcbcd6cba7def1818d8ce8e`, with the 349-model
canonical catalogue. The recorded preceding full P192/73 profile used 12,477,942,149
instructions in round 1, of which 832,131,218 were in the parent and 11,645,810,931
in the constructor child. Parent scalar SHA-256 accounted for 647,091,865 instructions
(77.76% of that parent, 5.1858% of the whole command). This historical profile used
345 models; it locates a cost but is not the matched reference for this new experiment.

Test a build-generated compact execution index, a build-time checksum of the exact
embedded catalogue, and the pinned native RustCrypto SHA-256 implementation for
receipt hashing. Preserve the public SHA-256 values, full inventory, exact curve
parameters, generators, cofactors, ambiguity rules, map output and independent replay.
These are engineering changes to command setup and hashing, with no new mathematical
algorithm or ECDLP performance claim. Construction/verification arithmetic is unchanged.

The primary success condition is at least 1% fewer median instructions for the entire
P192/73 search than the matched 349-model reference, with identical certified maps.
Screening and verification are secondary command diagnostics. Retain regressions and
inconclusive outcomes; do not discard a panel because a different panel improves.

## Frozen inputs and accounting

- Exact 349-model registry at the reference revision; record its SHA-256 and byte count.
- Existing public P224/1471 construction output, independent certificates and seed 1.
- Entire screening command: P224, prime degrees 1010–5000, maximum order 8.
- Entire independent verification command: P224/1471; construction excluded.
- Entire P192/73 search: constructor child, setup, file hashing, output and replay included.
- Five interleaved reference/candidate rounds per panel, same Rust 1.93.1 compiler,
  portable release flags and ARM Linux Docker image; Callgrind Ir summed over all
  traced children. No native macOS wall-time claim.
- Each invocation is limited to 1800 seconds, one VM CPU and 4 GiB container memory.
  Preserve stdout, stderr, profiles, exit status, source archives, executable hashes,
  toolchain/environment receipts and timeouts. Do not profile during builds or tests.

Compare full inventory JSON, normalized screening candidates, independent certificate
records, and byte-for-byte P192/73 map outputs. Independently compare accelerated
SHA-256 to the repository's unchanged implementation across block/padding boundaries
and large input. Confirm public subgroup and coefficient-mutation rejection remain.
Run dedicated CLI tests, field reference tests, standalone algorithm tests and the
required root library gate before pushing. Root failures remain explicit blockers.

## Requirement status before execution

| Requirement | Current status and continuation |
| --- | --- |
| Dedicated installed binary | Verified at the reference; repackage the passing candidate |
| Pinned curve inventory | 349 models; preserve all aliases, identities and operation capabilities |
| Broader prime construction/replay | 174 construction models, 146 replay models; retain existing support |
| Montgomery/Edwards and binary/extension maps | Partial: construction/replay adapters remain open requirements |
| Improved complete-search cost | Not established at the reference; execute the matched full-command panel |
| Root compile and hosted CI | Root has 643 recorded compiler errors; hosted jobs remain queued without an assigned runner |
| Merge and ReleaseMe publication | Blocked until applicable build, review and CI gates pass; do not bypass them |

Deliver a native-generated report, diagram, PDF, compact raw evidence and manifest in
this directory; update PR 1606. Leave the original study reports and catalogue graphs
unchanged when there is no new curve or map finding, and explain that in the report.
