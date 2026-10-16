# Frozen N83 known-witness stage control

This run probes whether the pinned native macOS arm64 worker reaches ordered
SAT search on the full retained K1182/K2048 domains with a known-answer
encoding. It is a stage-localization control, not a natural relation-yield
measurement or a complete index-calculus runtime comparison.

The [authorization](authorization.json) applies the user's 2026-10-14 grant
conservatively as **900 new local seconds** after 22,200 previously charged
seconds. The immutable [launcher](../known-witness-stage-launch-preflight-20261014/launch.py)
pins the worker, source, checker and domains by SHA-256, sets a **840-second**
process-group watchdog and a **5,000-conflict** cap, and charges the full 900
seconds even if the worker exits early. No second invocation is covered by
this receipt.

Run only after the local release-test processes finish. Launch into a fresh
directory on `/private/tmp`, where there is sufficient free space; copy the
complete raw launch, exit, stdout, stderr and audit receipts to `evidence-01/`
here without rewriting them. Recheck the copied evidence with the launcher's
`--check` mode and record the local disk and competing-load observations.
Retain any timeout, producer error or malformed trace as such. The audit's
`natural_relations`, `total_index_calculus_runtime_ms` and
`selected_best_total_runtime` fields must remain unset for this control.
