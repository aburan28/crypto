# Native prime-field base/oracle tuning

Status: **draft**. No experiment was run or scheduled. This card records the
next bounded engineering experiment after historical evidence PR #1898.
See the [evidence-linked handoff](../../research/factor_base_evaluation_20261015/followup_20261016/README.md),
its diagram and PDF for the decision sequence and sources.

## Question and decision value

Can joint selection of factor-base size and an exact native decomposition
oracle reduce complete one-target IC cost against the current verified
incumbent on fresh prime-order prime-field instances?

A robust holdout improvement becomes the next baseline, even above rho.
A null or pricing-sensitive outcome retains the incumbent. Low-x is a
historical control, not an assumed winner of the current native pipeline.
Coset/tower and random bases are matched-size controls where admissible.
This is a separate bounded engineering task; it does not replace the open
PKM reduced-presentation degree check or claim NIST-size recovery.

## Work before ready

- Audit existing native ecbench implementations and select the best verified
  per-curve pipeline; avoid duplicating already-landed oracle optimizations.
- Freeze exact registered curves, subgroup/target identities, seeds, base-size
  grid, arity, train/holdout split, repetitions and resource caps in PROTOCOL.md.
- Use archived inputs as regression fixtures only. Use fresh public holdout
  targets and prohibit known scalar labels from base selection.
- Pin a merged implementation commit, resumable command and native independent
  checker/replay command. The draft argv is a placeholder, not a command to run.
- Emit all unpriced counters on every arm, including rho, and compare at common
  weights 0, 0.1 and 1 GAE per native operation. Freeze the success/stop rules.
- Run `backlog validate` and record its receipt before changing status to ready.

## Resources and cost

The one-core, 16 GB RAM, 10 GB disk and two-hour timeout fields are proposed
screening ceilings, not predictions, measured costs or authorization to provision
a host. Runtime estimates are deliberately null until a protocol specifies the
grid and a correctness smoke run supplies a cost basis. Overflowing the cap is
a recorded incomplete result, never a reason to silently shrink the holdout.

Funding is proposed as zero while draft: no current queue priority changes
and no reward is posted. A reviewer can set funding relative to the existing
131,072-1,000,000 seed cards when a runnable protocol and cost are available.

## Completion

Land the frozen protocol, source/input hashes, all outcomes and exclusive phase
costs, independent relation/rank checks, verified scalar recovery and replay
receipts in a fresh `research/prime_base_oracle_tuning_<run-date>/` directory.
Include baseline/candidate and matched-rho ratios, the common native-price band,
paired uncertainty, memory, limits and failures. Add reports, benchmark copies
and affected dashboard updates in the result PR. Draft planning alone does not
complete this experiment.

Local availability: Conductor, Cargo and `backlog validate` returned 127;
`.conductor/WORKFLOW.md` returned 404. JSON-schema validation is a structural
check only; it does not replace native backlog validation or freeze a protocol.
