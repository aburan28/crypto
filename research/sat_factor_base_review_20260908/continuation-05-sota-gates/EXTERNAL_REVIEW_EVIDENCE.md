# Evidence accompanying an external review

Use this guide with [the review template](EXTERNAL_NOVELTY_REVIEW_TEMPLATE.json).
A source review, an independent execution, and a replay of retained outputs
establish different properties. Identify which work the reviewer performed;
an invitation, a completed template, or a `CONCUR` verdict alone does not
establish that the seven gates passed.

The added template fields describe a human review packet. They are not an
automatic admission schema, and they do not change archived execution records.
Keep unassessed fields null and evidence lists empty until evidence exists.
Preserve the existing Magma artifact bindings when reviewing that lane.

## Evidence records and coverage

Add one object to `reproduction_evidence` for each independently identifiable
backend/instance observation. Use these fields; record unavailable information
as null with an explanation, rather than zero or an inferred success.

| Field | Required meaning |
|---|---|
| `evidence_id` | Unique identifier within this return packet. |
| `evidence_kind` | `external_execution`, `retained_output_replay`, or `source_review`. A replay does not become an independent execution by using another host. |
| `backend` and `instance_id` | Actual backend and archived instance covered; a review spanning several instances must list them explicitly. |
| `parameters` | n, ell, m, curve coefficients, field representation, factor-base construction and actual size, target-population and sampling rule. |
| `source_binding` | Repository commit, relevant clean/dirty state, instance/source-system identity, factor-base/predicate identity, and input artifact hashes. |
| `tool_and_host_binding` | Tool versions and executable hashes, OS/CPU, requested and observed worker/thread counts, with evidence for observations. |
| `outcome` | Observed terminal result and its scope: source-system result, validated point decomposition, inconclusive cap/timeout, operational failure, or not executed. |
| `validation_evidence` | References to source-assignment validation and rational-point/signed-sum validation, assessed separately. A proper ideal without a point witness does not certify a PDP solution. |
| `resource_evidence` | References to raw process wall/user/system/total CPU, peak RSS and its scope, conflicts when exposed, measurement method, failed attempts, and enclosing receipts. |
| `artifacts` | Paths or immutable archive URLs, byte sizes, SHA-256 hashes, and the role of each raw output, metric, witness or review file. |
| `limitations` | Missing fields, contradictions, caps, unpriced resources, provenance uncertainty and other limits on this observation. |

List each omitted requested backend/instance combination and its reason in
`coverage_gaps`. Link every non-null gate assessment to the relevant records
through `gate_evidence_ids`. Record the scope of a partial assessment in the
assessment text. Five replicates in one cell do not establish coverage of a
different cell. An assessment-only report is useful even with no execution
records; it leaves reproduction requirements open.

Do not include credentials, license contents, private data or environment
dumps. Public professional identity, affiliation, conflicts, and the reviewer's
independence statement belong in the existing reviewer fields. Project-authored
CI and internal agent review remain identified as such, even when their outputs
agree on separate machines.

## Matching the literature

Fill `literature_matching` using exact sections/pages and the source copies
identified in `primary_sources`. The existing
[literature comparison](LITERATURE_MATCHING_20260910.md) records three material
differences for reviewers to assess:

- The campaign's n31 standard and GGMP cells use different curve coefficients
  and actual factor-base sizes: a=1 with 31 points versus a=0 with 63 points.
  Their timing ratio does not isolate the factor-base construction's effect.
- The campaign's n41 ell5 m3 cell differs from the prominent WDSat n41 ell20 m2
  experiment. Matching n alone does not reproduce that benchmark. GGMP's n41
  fraction examples also must not be called Frobenius-invariant constructions.
- A balanced SAT/UNSAT sample measures conditional behavior. Its imposed 50/50
  allocation does not estimate natural relation yield; planted multiplicities,
  the target population and P/-P source-system dependence need separate treatment.

The novelty assessment includes the 2019 WDSat, 2020 GGMP and 2025 SATIC sources.
Fill retrieval/version, digest and section fields from the source actually
reviewed. A citation supplied by the project is not evidence that the reviewer
retrieved or assessed it. Record unavailable sources as a review limitation.

## Costs and the supported claim

The historical meter sets both `single_core_seconds` and `total_core_seconds`
to user plus system CPU. They are aliases for one observation; do not add them
or treat the former as independently measured single-core elapsed time.
Requested thread settings do not prove observed single-core execution. State
the actual observation and its limits, keeping elapsed time separate from CPU.

Identify inclusive and nested receipts before summing costs. Sum disjoint
process charges across all workers/hosts once, and report the inclusive wall
interval separately. Retain preparation, discovery, predicate construction,
builds, collection, linear algebra, validation, failed attempts and external
resources with their accounting boundaries. Do not add an enclosing CPU timer
to its child CPU charges. A maximum per-process RSS does not establish the
simultaneous memory peak of a parallel process tree or multiple hosts.
Missing conflict counters remain unavailable, not zero.

An external Magma run can supply missing Magma evidence. A source review can
support correctness or prior-art classification. Neither alone establishes
the full comparison, larger-regime index-calculus scaling, or a SOTA result.
The final verdict must name the narrowest supported claim and unresolved gaps.
