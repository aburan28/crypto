# Research problems and obstructions

Track the major named blockers that prevent a research claim, proof step,
implementation, experiment, or audit conclusion from being established.
The [agent rule](../AGENTS.md#name-the-obstruction) requires each record to
explain the gap and what would let the work move past it.

## Index

Keep resolved and superseded records in the index.

| Problem | Status | What it blocks | What would unblock it |
| --- | --- | --- | --- |
| [Root-finding pipeline complete-cost calibration](rootfinding-complete-cost-calibration.md) | open | qualified full-cost or wall-time comparison for the assessed root-finding IC pipeline; 4/8/16-lane, B train/holdout and coverage panels are verified and replayed; a fresh-target equal-prefix component control exists, while native field counts are absent from sealed records and wall data are L0 | serialize and price native work; repeat matched public targets under L2 isolation and on an independent host |
| [EC Bench GPU matched workload and validation](ecbench-gpu-matched-workload-contract.md) | open | admitted CPU/GPU comparison of the same ECDLP workload | frozen CUDA input/output, exact costs, independent answer replay and GPU contention controls |
| [N83 recorded-query sampler contract](n83-recorded-query-sampler-contract-mismatch.md) | open | green whole-library validation on the pinned baseline | traced and versioned query law, point checks and passing regression |
| [GHS → Smith characteristic mismatch](ghs-to-smith-pipeline-characteristic-mismatch.md) | open | the GHS → Smith → Diem pipeline | a char-2 hyperelliptic→quartic isogeny, or an odd-char descent onto genus-3 hyperelliptic curves |
| [Smith 18.57% walk extension](smith-kernel-isogeny-class-walk.md) | open | costing walks that lift Smith's coverage toward 50.16% | native Smith Prop. 2/4 tests and per-class success counts |

## Add or update a problem

1. Check this index and the affected work for an existing record. Reuse it when
   the obstruction is the same; link the canonical record across repositories
   rather than duplicating its evidence.
2. Copy [TEMPLATE.md](TEMPLATE.md) to a stable, descriptive filename such as
   `subgroup-preserving-transfer-not-established.md`. Name the mechanism or
   missing obligation, not a vague "hard problem". The filename above is an
   illustration, not a recorded finding.
3. State the exact blocked task or claim and its parameters, assumptions, and
   acceptance criteria. Separate evidence from conjecture and cite durable
   artifacts, proofs, reproductions, or source records.
4. Specify the evidence or capability needed to move past the obstruction and
   one concrete next experiment or proof obligation that can distinguish
   explanations. Link an existing execution protocol where applicable; this
   record does not replace it or impose a new campaign spending cap.
5. Add the index row and link the problem from the affected note, plan, audit
   finding, or PR. Update these pointers and the record together as work lands.

Use this directory for consequential or recurring blockers. Routine defects
and TODOs can stay in their existing trackers. Record a blocker even when its
cause is unknown; explicit unknowns are useful, invented explanations are not.

## Status and evidence

- `open`: the obstruction remains and its next action is recorded.
- `investigating`: work on that action is underway; link the task or run.
- `resolved`: linked evidence meets the stated exit criteria for this specific
  obstruction. State any remaining limitations.
- `superseded`: a linked successor now tracks this obstruction; this is not a
  claim that it was solved.

These statuses describe the problem record, not the scientific status of a
hypothesis or the completion of a campaign. Stopping work does not resolve
a problem or authorize abandoning the user's requirement.

Keep observations, calculations, hypotheses, resource/access/assistance gaps,
and proved limits distinct. Bound a negative result by its exact inputs and
search limits. A resource ceiling, missing implementation, or failed search
does not establish a mathematical obstruction or a security theorem.

Preserve historical evidence and append dated updates explaining corrections.
A resolution cites its evidence and the criteria met; a supersession cites its
successor and the reason for replacing the record. Follow the repository's
existing review, archival, and claim-promotion rules.
