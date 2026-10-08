Open pull requests ready for review, never as drafts, whatever the runtime
defaults to, and mark an existing draft ready; only a user request in the task
makes a draft. This and the repository's other rules are in
[`AGENTS.md`](AGENTS.md) ("Default workflow").

<!-- conductor:begin -->
Before making code changes, obtain or attach to a Conductor task. Run
`conductor check --summary "…" --scope path:…` first — if someone already holds
those files, it will tell you who and what to do about it. Read `.conductor/WORKFLOW.md`
and the active task card. Report scope expansion before editing outside the reserved paths.
Do not publish chat transcripts or secrets as task metadata.
<!-- conductor:end -->

In Claude Code on the web, `.claude/hooks/session-start.sh` installs Conductor
and cairn and starts a control plane when the session starts. Elsewhere
Conductor is not preinstalled: [`docs/CONDUCTOR.md`](docs/CONDUCTOR.md) says how
to install it and connect a session to a control plane. Until then the hooks
warn and let edits through.
