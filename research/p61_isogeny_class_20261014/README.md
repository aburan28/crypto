# p61 selected-class isogeny census

This study records the certified launch evidence for the selected ordinary
isogeny class `p61-i0025` over `p = 2^61 - 1`.

The class order `358094304` is unconditionally certified by a PARI flag-0
receipt. The active producer/replayer campaign has per-segment admission through
`x=147776`, or 2,659,968 vertex records. The last sealed global uniqueness
ledger remains the 595,584-record `x=0..33088` prefix, so those two evidence
levels are kept distinct. Full class coverage remains incomplete.

The selected-degree adjacency producer and its distinct kernel-first verifier
now pass one-, 64-, and 4,096-vertex windows. The largest create-only admission
contains 28,672 directed edges and 471,040 kernel coefficients. Its observed
representation density projects to approximately 3.632 TB logical JSONL and
0.748 TB gzip for the complete four-degree graph. Production adjacency coverage
is 4,096 of 358,094,304 source vertices; the dense target-index join and
reciprocity audit remain future gates.

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

The adjacency calibration report, visual sources, machine-readable metrics,
failure records, and rendered PDF are in
[`reports/run-20261015T0436Z-adjacency-pilot/`](reports/run-20261015T0436Z-adjacency-pilot/):

- [Markdown report](reports/run-20261015T0436Z-adjacency-pilot/REPORT.md)
- [PDF report](reports/run-20261015T0436Z-adjacency-pilot/p61-adjacency-pilot-report.pdf)
- [Calibration graph](reports/run-20261015T0436Z-adjacency-pilot/adjacency-calibration.svg)
- [Evidence pipeline](reports/run-20261015T0436Z-adjacency-pilot/evidence-pipeline.svg)
- [Recorded metrics](reports/run-20261015T0436Z-adjacency-pilot/metrics.json)
- [Build attempts](reports/run-20261015T0436Z-adjacency-pilot/build-attempts.json)
- [Test and replay results](reports/run-20261015T0436Z-adjacency-pilot/test-results.json)
- [Render attempts](reports/run-20261015T0436Z-adjacency-pilot/render-attempts.json)

`run_adjacency_campaign.py` is the fail-closed adjacency follower. It accepts
only independently admitted vertex segments, creates nonoverlapping windows,
replays source and output hashes, requires distinct producer/verifier releases,
checks exact degree and kernel totals, and admits a window only after replay and
a host sync barrier. The active contract reserves at least 1.5 TB free and caps
logical adjacency output at 6.0 TB.
