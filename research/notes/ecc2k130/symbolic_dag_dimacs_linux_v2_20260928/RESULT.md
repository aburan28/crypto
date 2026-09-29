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
returned `HASH_ONLY_NO_MEASURED_CHILD`; all nine no-network one-shot/refusal
and relation-mutation controls passed; `actionlint` and Python compilation passed. A local attempt
to enter the measured supervisor under this held freeze produced zero phase
children and a replayable `ARCHIVED_PRE_DISPATCH_REFUSAL`. This refusal
control is outside the PR evidence namespace and is not attempt 2.

This held revision strengthens replay of a future n131 PASS: it reconstructs
the pinned relation and checks each canonical gate clause plus the output
assertion. A dimension-preserving gate mutation, an output-unit mutation, and
the formerly accepted synthetic `p cnf 1 1` are regression controls. The
measure workflow now archives checkout/setup/static preflight refusals before
the supervisor. Hosted harmless control on the target image showed that
`git merge-base` exits 128 under the inherited toy 512-MiB address-space cap:
Git cannot map the checkout packfile. Under the same cap, the Go-based `gh`
CLI exits 2 before reading a PR because it cannot reserve page-summary memory.
The full Git/GitHub/preparation gate
therefore stays in the uncapped supervisor before dispatch; capped children
verify its sealed dispatch and local byte identities. CI probes both paths
without starting a measured phase. The hosted control receipt is evidence of
this host limitation, not an attempt-2 outcome.
