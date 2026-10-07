<!-- conductor:begin -->
Before making code changes, obtain or attach to a Conductor task. Run
`conductor check --summary "…" --scope path:…` first — if someone already holds
those files, it will tell you who and what to do about it. Read `.conductor/WORKFLOW.md`
and the active task card. Report scope expansion before editing outside the reserved paths.
Do not publish chat transcripts or secrets as task metadata.
<!-- conductor:end -->

Conductor is not preinstalled everywhere; a cloud session has it only when its
environment's setup script installs it. [`docs/CONDUCTOR.md`](docs/CONDUCTOR.md)
says how to install it and connect a session to a control plane. Until then the
hooks warn and let edits through.
