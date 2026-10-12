# Static-source wide-S4 SAT controls: four verified point witnesses

The six [registered disclosed controls](panel.json) completed once each.
The source-receipted static CryptoMiniSat executable passed `--version` from
its final copied run path before any instance was exported. All six pinned
Rust exports were valid wide-S4 XOR-DIMACS source systems. The solver reported
UNSAT on both exact-negative controls and returned a complete SAT assignment
on each of the four exact-positive controls. An independent replay checks
**every clause and XOR constraint** of each positive assignment and re-adds
the decoded factor-base points to the supplied public point. All four are
valid group witnesses. The two negative solver reports are consistent with
the separate exhaustive group oracle; no independently checkable UNSAT proof
trace was emitted.

| Trial | Parent exact relation | SAT status | Source assignment and group replay | SAT process wall |
| --- | --- | --- | --- | ---: |
| 0 | absent | reported UNSAT | no model | 21.262 s |
| 3 | absent | reported UNSAT | no model | 30.598 s |
| 4 | present, bounded F5 missed | SAT | valid point witness | 2.348 s |
| 67 | present, bounded F5 missed | SAT | valid point witness | 8.369 s |
| 71 | present, bounded F5 missed | SAT | valid point witness | 12.797 s |
| 75 | present, bounded F5 missed | SAT | valid point witness | 18.732 s |

These process times are **stage diagnostics on outcome-selected controls**.
They include each external SAT child process only; source export, setup,
relation collection, final matrix/rank/LA, target descent, and recovery are
separate or absent. The geometric base contains 63 points, of which 62 are
subgroup-usable; its 29 folded relation columns are properties of the parent
F5 full worker, not an end-to-end result of this SAT pilot. There is no
natural-yield estimate, candidate `IC1` run, one-target IC online interval,
paired incumbent/rho measurement, or speedup claim. The four positives were
selected *because* bounded F5 had missed known relations; their 4/4 result
must not be reported as a natural success rate.

The complete [raw evidence archive](evidence.tar.gz) is SHA-256
`3eb9c693bcc00519b4f94b39cab4b96adc5e8e67711a34e670e90f27e3106e97`.
It retains 81 files: six exact input exports and manifests, process commands,
metrics, SAT streams, the exporter source/build receipt, static solver
executable and Phase-B source/dependency receipt, preflight and host data,
registered code and all rows. `summary.json` is SHA-256
`13a4b0203115756547236cea57b3abf760720109423af63d98ad07b89eff60a6`.
The [independent archive audit](../../test_static_cms_s4_evidence.py)
checks the sealed bytes, source/build identities, full formulas and four
point sums. The Phase-B terminal bundle verifier itself passes on the
unchanged source-receipted static CryptoMiniSat build.

The next scientific gate is a **new frozen full-worker protocol** with this
encoding/backend wired into ordinary natural relation collection and target
handling. It must charge every failed PDP attempt, construct and solve the
actual relation matrix, recover and independently verify one unseen target,
and pair F5/SAT/incumbent/rho on the same new public points and resources.
These six controls are exhausted and will not be retried or treated as a
performance sample.
