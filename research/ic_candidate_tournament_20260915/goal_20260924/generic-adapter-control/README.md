# Retained mixed-adapter tournament control

The [registered control](../generic-qualification-adapter/PROTOCOL.md) passed in
[workflow 36280983990](https://github.com/aburan28/crypto/actions/runs/36280983990),
at commit `717879287fa389c1f2eb3817316446dc86664f97` in
[PR 862](https://github.com/aburan28/crypto/pull/862). The existing tournament
completed 90 verified native/profile pairs and retained all 12 intentionally
incomplete profiles. Its frozen evaluator reconstructed all 102 trial receipts
and 969 source files, both in Linux CI and after artifact transport to macOS.
No worker was executed during either replay.

This is an accounting and integration result on one toy cell, **not a reference
selection, observer qualification, improvement round, or speedup claim**. The
accepted incumbent remains unchanged. One improvement round is complete without
a qualifying winner; two remain.

## Frozen workload and outcomes

The field is F_(2^13), curve parameter a=0, with subgroup order **2003**; field
degree is not subgroup bit length. Seed 2026092662 generated five supplied public
points: one A/A point, one smoke point and three development points. There are
three process repetitions, one target per job, one pinned CPU and Rayon worker,
8 GiB, and a 180-second process limit. Linux amd64, Rust 1.94.1 and Valgrind 3.22.0
are recorded in the frozen contract. No selection, confirmation or replay-stage
targets were generated. Here, *frozen replay* means auditing saved evidence.

The requested subgroup-orbit recipe was seed43, points78. Each actual IC base
contains **182 distinct usable points and seven folded columns**. Both adapters
independently check that support. The one-attempt budget arm cannot complete the
full-rank preparation, so its native solve is never launched.

| Arm | Verified pairs / planned slots | Retained incomplete profiles | Actual base / folded columns | Performance claim |
| --- | ---: | ---: | ---: | --- |
| Prepared scaled incumbent, including A/A | 15 / 15 | 0 | 182 / 7 | Unqualified control |
| Prepared A/A copy | 3 / 3 | 0 | 182 / 7 | Unqualified control |
| Generic pair-table, dense LA | 12 / 12 | 0 | 182 / 7 | Unqualified control |
| Generic pair-table, sparse LA | 12 / 12 | 0 | 182 / 7 | Unqualified control |
| Generic dense, one-attempt budget | 0 / 12 | 12 | 182 / 7 | Incomplete; costs unknown |
| Prepared rho, width 1 | 12 / 12 | 0 | Not applicable | Unqualified control |
| Prepared rho, width 4 | 12 / 12 | 0 | Not applicable | Unqualified control |
| Generic rho, width 1 | 12 / 12 | 0 | Not applicable | Unqualified control |
| Generic rho, width 4 | 12 / 12 | 0 | Not applicable | Unqualified control |

[RESULTS.json](RESULTS.json) maps aliases to canonical candidate or reference IDs.
All paired arms share one workload ID per point. There are 162 distinct stored
record IDs: prepared records describe a pair, while generic records distinguish
native and profiled processes. Twelve generic native IDs reserve work marked
`NOT_RUN`; **162 is not an executed-process count**. All twelve incomplete rows
retain identity and failure evidence, with null verified costs and native timing.

Native one-target online time remains the primary research metric. Complete cold
native time and Callgrind Ir are supplementary. The raw control records preserve
their clocks, closed phases, instructions, certificates and development statistics;
none of those wiring diagnostics selects a research reference or supplies a
qualified speedup or normalized-S claim. All observer costs remain charged.

## Durable evidence and replay

[EVIDENCE.json](EVIDENCE.json) lists the archive SHA-256, every member's hash/size,
original artifact identity, source commit and workflow. The archive contains
8,374 files and 48,688,480 expanded bytes, including source/build identities,
executables, fixtures, admissions, failed reports, unchanged compressed profiles,
process receipts, summaries and the frozen evaluator. Only disposable build and
bytecode caches and the transient lock are omitted from the control subtree.
The original 143,283,824-byte CI ZIP additionally contains other driver controls
and a build cache; those are not silently presented as part of this control.

The durable gzip archive is 27,761,146 bytes with SHA-256
`52a3f197f47cbb55558bb67497518bb4fe1440c9e78ce515822f82019a28b48e`.
It preserves measured bytes; archiving does not rerun an experiment.

With Python 3.12 and a new output path:

```sh
python3.12 research/ic_candidate_tournament_20260915/goal_20260924/generic-adapter-control/audit_archive.py \
  --out /tmp/ic-generic-adapter-audit.json
```

The script checks the archive and complete member set, extracts into a temporary
directory, runs the **archived** evaluator, checks failure/identity/base/workload
invariants, and reproduces the result export. Linux CI performs strict frozen
statistics replay; the portable archive/export test runs on both platforms.
The earlier macOS transport replay is retained in
[transport-replay.json](transport-replay.json).

[fixtures.json](fixtures.json) is an exact copy of the five exposed points. Pass
it to the existing `--exposed-fixtures` mechanism before any later confirmation,
together with all other supplemental exposures and prior round history.

The next research gate is a registered five-cell generic/reference comparison
and an explicit observer study. Keep generic arms out of the two remaining
competitive rounds until those gates are audited. Do not redispatch round one.
