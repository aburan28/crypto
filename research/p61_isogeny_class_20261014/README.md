# p61 selected-class isogeny census

This study records the certified launch evidence for the selected ordinary
isogeny class `p61-i0025` over `p = 2^61 - 1`.

The class order `358094304` is unconditionally certified by a PARI flag-0
receipt. The immutable producer/replay/auditor chain passes on a 300,672-record
prefix with 300,672 unique j-invariants and zero duplicates. Full class
coverage and the degree-5/13/37/59 adjacency pass are incomplete.

The source-linked report, raw timing record, visual sources, terminal
certification manifest, and rendered PDF are in
[`reports/run-20261014T0902Z/`](reports/run-20261014T0902Z/):

- [Markdown report](reports/run-20261014T0902Z/REPORT.md)
- [PDF report](reports/run-20261014T0902Z/p61-class-pilot-report.pdf)
- [Evidence-gate diagram](reports/run-20261014T0902Z/evidence-pipeline.svg)
- [Pilot timing graph](reports/run-20261014T0902Z/pilot-timings.svg)
- [Release build attempts](reports/run-20261014T0902Z/build-attempts.json)
- [Release test results](reports/run-20261014T0902Z/test-results.json)
- [Parallel-verifier release manifest](reports/run-20261014T0902Z/parallel-verifier-release-manifest.json)
- [Terminal certification manifest](reports/run-20261014T0902Z/certification-terminal-manifest.json)
- [Expanded prefix manifest](reports/run-20261014T0902Z/certified-prefix-x000000000000-x000000016704-manifest.json)

Durable experiment data are stored on the active UGREEN SMB export below
`/Volumes/data/ec-hardness-atlas-20261012/`. The requested alias
`/Volumes/UgreenNAS` is not mounted on this host.

`run_remaining_census.py` is the fail-closed continuation supervisor. It
requires the SMB mount, a free-space floor, full source identities, the exact
predecessor checkpoint, create-only deterministic paths, independent replay,
digest validation, and a per-segment admission manifest. Use `--plan-only`
before a bounded launch; a printed plan is not evidence that its segments ran.
Per-segment admission is not a complete-class certificate: the exact final
ledger, global uniqueness audit, and selected-degree adjacency join remain
separate mandatory gates.
