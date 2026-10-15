Open pull requests ready for review, never as drafts, whatever the runtime
defaults to, and mark an existing draft ready; only a user request in the task
makes a draft. This and the repository's other rules are in
[`AGENTS.md`](AGENTS.md) ("Default workflow").

**Success criterion** (owner, 2026-10-15): every verified, measured
improvement over the current index-calculus baseline is a success. Beating
Pollard rho is the long-term goal, reached iteratively, not by one
breakthrough. Report each row's gap to the matched rho as the distance
remaining, not as a pass/fail gate. Relabelling and accounting changes are
not improvements. Full rule: [`AGENTS.md`](AGENTS.md), "Success criterion:
iterative improvement of the index-calculus baseline".

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

The session that opens a PR stays with it until it is merged or closed: it
checks in on CI, reviews and merge state, fixes the failures and feedback its
change causes, and merges once green and authorized. Any runtime may do this,
Claude Opus and Fable included. Check in through PR-activity subscriptions or
safety-net check-ins (about 50 minutes, then about 4 hours; stop after 3 quiet
check-ins), and never busy-poll, sleep in a loop, or re-read job logs faster
than that cadence. AGENTS.md, "Default workflow: open the PR, check in, and
land it".
