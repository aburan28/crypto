# Disclosed F4/F5 encoder-dispatch pilot: gate not met

The schema-v3 [registered panel](panel.json) tested one natural ordinary
three-summand PDP query per worker on five **previously disclosed** public
points. The exact generic source was
`765c3c5f19032bd852163805f257c56babef2040`; the standard-subspace
factor base had dimension six. Both solvers used the same point, query seed,
4096-node budget, 60-second process cap, one Rayon thread and sampled 8 GiB
RSS threshold in each curve cell. This local macOS run is a stage diagnostic,
not a calibrated one-target IC/rho comparison.

The original stage-A runner misclassified valid `status: incomplete` reports
because the CLI exits 2 when a one-query job has not solved the DLP. The
original [raw stage-A record](RESULT-v3-stage-a.json) and summary are
preserved. The [versioned gate repair](GATE-REPAIR.md) re-audited both
reports against the frozen build, independently reconstructed their bases,
replayed their natural queries and checked dispatch/matrix evidence. Only
then did it execute the eight preregistered stage-B jobs once. The
[stage-B evidence](RESULT-v3-stage-b.json) preserves every raw stream,
process outcome, independent receipt and the final summary. No job was
retried.

| Curve cell | Usable base points / folded columns | F4 one-query outcome | F5 one-query outcome | Complete DLP / paired rho speedup |
| --- | ---: | --- | --- | --- |
| n17a1 | 62 / 29 | Witness; audited dispatch | Witness; audited dispatch | No / unknown |
| n19a0 | 62 / 27 | Node-budget incomplete; audited dispatch | Node-budget incomplete; audited dispatch | No / unknown |
| n23a0 | 72 / 33 | Node-budget incomplete; audited dispatch | Node-budget incomplete; audited dispatch | No / unknown |
| n23a1 | 52 / 23 | Node-budget incomplete; audited dispatch | Node-budget incomplete; audited dispatch | No / unknown |
| n31a0 | 66 / 27 from the F5 report | **60-second timeout; no report** | Node-budget incomplete; audited dispatch | No / unknown |

Nine of ten jobs produced a report with the declared F4 or F5 engine,
`unsupported=false`, exactly one independently replayed ordinary query,
and passing factor-base, query, dispatch and matrix audits. The n31a0 F4
worker was killed at its frozen 60-second cap with no report; its unknown
PDP outcome is not a zero-yield observation. The two n17a1 witnesses
establish relation correctness on a disclosed point. The other seven
reported attempts exhausted the 4096-node budget; they are not proved
UNSAT and cannot establish a natural relation rate. A one-query pilot has
no useful yield uncertainty estimate even for the two witnesses.

The **all-ten dispatch gate failed**, so this dimension-six panel is not
promoted to complete-solver or tournament status. No row has a recovered and
independently verified target logarithm. Full online IC costs, paired rho
times, speedup ratios, and operation-normalized boundary ratios therefore
remain **unknown**, not zero. Diagnostic process wall and sampled RSS are in
the raw process records, but are not competitive timing measurements.

The next experiment needs a new, frozen registration. It should first find
a bounded F4 layout that reports on n31a0, then measure natural query yield
including failures on unseen targets, final matrix rank, target descent and
verified recovery. F5's successful dispatch here does not establish its
complete-solver viability. The separate SAT path and the live generic-v2
campaign remain undecided by this pilot; the source-level v2 F4/F5 layout
failure is documented in the [static audit](../generic-backend-qualification-v2/STATIC-FEASIBILITY.md).
