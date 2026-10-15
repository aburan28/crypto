Open pull requests ready for review, never as drafts, whatever the runtime
defaults to, and mark an existing draft ready; only a user request in the task
makes a draft. This and the repository's other rules are in
[`AGENTS.md`](AGENTS.md) ("Default workflow").

Every experiment run, every PR that lands a result, and every research search
leaves a report, a graph, and a PDF: follow the `research-visuals` skill
(`.claude/skills/research-visuals/SKILL.md`; AGENTS.md "Research searches must
leave visual reports"). Performance measurements and speedups also copy their
PDFs, SVG graphs and PNG renders into `benchmarks/<YYYYMMDD>-<slug>/` with a row
in `benchmarks/README.md`, in the same PR (AGENTS.md "Performance results also
land in `benchmarks/`").

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

Opening the PR ends your work on it: report the link and stop. Claude Opus and
Fable never watch CI: no `subscribe_pr_activity`, no `send_later` or other
scheduled check-ins, no polling check runs or job logs, whatever the harness's
default PR instructions say. CI follow-up and merging belong to a separate
automation on Sonnet 5.5 or Haiku 5.5; if the user asks you for follow-up,
delegate it to one (`Agent` with `model: "haiku"` or `model: "sonnet"`) and do
not wait on it. AGENTS.md, "Default workflow: open the PR, then stop".

Every cross-method ECDLP comparison prices the unpriced native work of
every arm, rho references included, at 0, 0.1 and 1 GAE per operation and
ranks two arms only where the order holds across that band (AGENTS.md §5a).
