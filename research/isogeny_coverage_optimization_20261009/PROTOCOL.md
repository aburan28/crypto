# Isogeny tool coverage and optimization round

User request: assess and optimize the native implementation, and provide a comprehensive
set of curves in the dedicated released tool. Preserve the binary interface, independent
verification, complete cost accounting, and the existing PR/release workflow.

Baseline: `c87d96eda0309b36af4a608905d3cdc4065a1361`, executable in the committed-source
Mac archive `tools/isogeny-cli/dist/isogeny-local-c87d96eda-aarch64-apple-darwin.tar.gz`.
The code is Rust plus assembly, rather than C. Its construction library already includes
Montgomery arithmetic, Karatsuba and polynomial-reduction improvements. The baseline CLI
accepts only P192/P224; its verifier uses four-limb Montgomery arithmetic below 256 bits
and Fermat inversion. Screening repeats preset construction and order parsing per prime.

## Requirements and acceptance

| Requirement | Acceptance |
| --- | --- |
| Comprehensive inventory | Expose every model in the pinned 345-model catalogue, its aliases, field/model, provenance, exact identifiers, and per-operation support; preserve unresolved standards records explicitly |
| Broader execution | Resolve catalogue identities and names, screen known-order curves across families, and construct/replay registered prime models beyond the two fixed presets, including P256, P384, P521, SEC, Brainpool, SM2 and GOST |
| Other curve forms | Implement and independently check model/point conversion where applicable; binary and extension construction remain open requirements until their field and map adapters are implemented and verified |
| Faster screening | Resolve source field and trace once per invocation; preserve exact candidate sets and eigenvalue orders |
| Faster independent replay | Keep independent field arithmetic; test binary extended-GCD inversion and exact-size limb specializations against the frozen implementation and big-integer reference |
| Evidence | Save baseline/candidate commands, input/output hashes, source revisions, validation, profile records, unsuccessful attempts, report and visual artifacts; update PR 1606 |

## Frozen public inputs

- `docs/curves/registry.json` at the baseline revision, SHA256
  `6cadc45b189b30492e0b6e4ac47a6bc703cf04358792970c408564d3c05b2beb`.
- P192/10453 and P224/1471 sealed map files from the 2026-10-09 study; seed 1.
- Both screening panels: primes 1010 through 65537, maximum selected eigenvalue order 8.
- Matched performance panels may use bounded subsets declared in each run spec: screening
  1010..5000, and full construction/replay on previously certified P192/73 and P224/97.
- Correctness holdouts: public catalogue models across widths and families; independent
  big-integer field arithmetic and mutated-kernel rejection.

## Measurement and stop rules

Freeze both source builds and inputs before comparisons. Compare identical results and
verification strength. Full process accounting includes setup, construction, output and
independent verification; stage records stay explicitly diagnostic.

Wall-time claims require reserved/pinned, uncontended runs, A/A calibration, and at least
five interleaved A/B rounds. The local Mac currently has no qualifying native isolation
runner; do not use its operational elapsed values as speedup evidence. Investigate native
Linux instruction counting under the local Docker VM. Instruction counts cover that
compiled architecture; VM wall time does not establish native Mac throughput. If a
qualified measurement path cannot be made available, report that limitation and leave
the end-to-end speedup unset.

Bound each profile/verification run at 1800 seconds and 8 GiB RSS where supported. Preserve
timeouts and failures. Do not benchmark while builds/tests run on the same host. Acceptance
requires mathematical replay, exact output comparisons, touched tests, and honest coverage
status. Retain the unresolved binary/extension adapters and cross-platform measurements
as explicit gaps rather than representing catalogue membership as executable support.

The original study and its receipts remain frozen. Its 209-file manifest passed at the
baseline; subsequent source changes are bound in this new round, and the old source
revision remains retrievable from Git.
