# Open-PR integration, 12 October 2026

This batch combines 37 eligible source heads captured at
`2026-10-12T02:31:34Z`. [heads.json](heads.json) records every exact commit,
the excluded draft #1690, and the later main snapshot reconciled during
authoring. The captured commits remain ancestors of the integration branch.

[validation.json](validation.json) records the local checks, their source
revisions, earlier failures, corrections, and remaining external gates.
The checks validate integration and reproduction; they do not create new
performance measurements or complete open research requirements.

- [Workflow reconciliation](workflow-reconciliation.json) describes the
  sealed Lean budget exceptions and the isolated historical n53 lock view.
- [Proof verification](proof-validation.json),
  [catalogue verification](catalogue-validation.json), and
  [Eigen-Séta verification](eigen-seta-validation.json) retain their specific
  commands, inputs, and limits.
- [Final metadata checks](metadata-final.json),
  [benchmark copy coverage](benchmark-coverage.json), and
  [main reconciliation](final-main-reconciliation.json) track the generated
  views, source-linked visual copies, and preserved F5 study snapshots.
- [Queue cleanup](queue-cleanup.json) records cancellation of 192 obsolete
  conventional QA runs. [Pages diagnosis](pages-diagnosis.json) is a dated
  observation, with the concrete post-merge publication checks still owed.

The `receipts/` directory contains portable logs and input hashes. Earlier
failed runs are retained and distinguished from their successful replacements.
GitHub checks, required reviews, the guarded merge commit, and live Pages
verification belong to the separate CI follow-up after this ready PR opens.
