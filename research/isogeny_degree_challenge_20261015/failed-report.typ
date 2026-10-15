#set page(paper: "us-letter", margin: 0.7in, numbering: "1")
#set text(font: "DejaVu Sans", size: 10pt)
#let run = sys.inputs.at("run")
#let freeze = json(run + "/freeze.json")
#let entry = json(run + "/controls/p192/RESULT.json")
#let child = json(run + "/controls/p192/judgement/verification.json")
#let policy = json(run + "/policy.json")
#align(center)[#text(size: 21pt, weight: "bold")[Prime Degree Frontier]\
Failed admission-control run · 2026-10-15T23:25Z]

= Requested and executed scope
Run the frozen native executable's all-standard-curve policy and candidate
ladder, then native admission controls. The first positive P192/13 judgement
was censored, so this run stopped. P521, the eight mutation cases, forged-score
leaderboard, and scheduling stress were not reached in this run. Their earlier
passing packet and later corrected rerun are separate evidence.

#table(columns: (1fr, 3fr), inset: 6pt,
  [Source revision], [#text(size: 8pt)[#freeze.source_commit]],
  [Executable SHA256], [#text(size: 8pt)[#freeze.binary_sha256]],
  [Host], [#freeze.host],
  [Roster], [#policy.curve_count standard prime models, all listed separately.],
  [Judge result], [#entry.status; score absent.],
  [Resource receipt], [#entry.verification.status],
  [Replay child output], [#child.status],
)

= Observations and cause
The child emitted a passing independent certificate, but the supervising parent
reported unavailable resource monitoring and withheld the score. Its receipt
retains a positive sampled RSS value, five unavailable samples and exit code 0.
This is an admission/supervision failure, not a mathematical rejection.

The source treated every zero-RSS read as unavailable even after the startup
handshake had succeeded. On the observed macOS shutdown transition that made
the parent censor an already completed child. The correction keeps the positive
sample requirement before the handshake, then treats a successful zero-RSS
read as a valid sample afterward. Genuine unavailable monitoring still fails
closed. The follow-up run must pass the full integration controls and twenty
fresh short-map judgements with positive recorded RSS before publication.

#figure(image("failure-flow.svg", width: 100%), caption: [Observed control flow: a passing child was censored by its resource supervisor.])

= Reproduction and retained evidence
The exact executable/source identities are in freeze.json and binary.sha256.
policy.json and plan.json retain the complete roster and candidate screen.
controls.log retains the failed assertion. The raw map, constructor receipt,
exact independent verification certificate, stderr, resource receipt and
censored entry are under controls/p192/. No score is assigned to this run.
The full-package failure is also retained in the validation archive.

Commands were run through the standalone installed binary from an installation
directory. Replay and construction had the published separate budgets. The
host was shared; operational receipts establish no isolated speed result.
The canonical IC/rho dashboard, leaderboard and curve-edge charts remain
unchanged because no new graphable mathematical or performance fact was found.
