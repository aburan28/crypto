# Linux v2 single-edge DAG-to-DIMACS gate

**HELD; measured attempt 2 has not run.** The merged #804 first cold receipt
is `PRODUCER_FAILURE / LAUNCH_ERROR` before its toy producer executed. The
merged #831 preparation is `PASS_HARMLESS_PREPARATION_ONLY`, with no toy, n13
or n131 measured child. This PR contains the proposed second release gate,
wrapper and independent replay, with release fields null and no outcome
archive. It reports no SAT/UNSAT, CNF capacity, cost, PDP, ECDLP or rho verdict.

Next action: review this exact held diff and hash-only CI. A separately
reviewed release commit must freeze this PR number and main head; only then
may the unique label initiate one bounded cold toy → n13 → n131 attempt. The
raw success or failure must be committed and independently replayed before
this result can change.

Held validation on the isolated merged-main checkout: inherited #831
`verify_preparation.py` PASS (18 raw files and both Linux build hashes);
`ci_replay.py` PASS with no outcome archive; both gate-only entrypoints
returned `HASH_ONLY_NO_MEASURED_CHILD`; all four no-network one-shot/refusal
controls passed; `actionlint` and Python compilation passed. A local attempt
to enter the measured supervisor under this held freeze produced zero phase
children and a replayable `ARCHIVED_PRE_DISPATCH_REFUSAL`. This refusal
control is outside the PR evidence namespace and is not attempt 2.
