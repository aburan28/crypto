# Repeat the D=-43 endpoint comparison with a generic degree-11 recurrence

Status: blocked. This card is a proposal, not a prepared or completed run. No cairn task, lease or objective was created.

## Question and decision

Does the ten-gap D=-43 query-count signal repeat on fresh public targets after adding the generic degree-11 recurrence control?

A fresh-target repeat distinguishes the current small-family signal from target/query noise. A repeated phase-specific effect merits a larger preregistered population; parity or censoring does not support a conductor-hardness claim.

## Baseline and protocol to freeze

Source: [the completed study](../../research/conductor_gap_cm_orbits_20261015/REPORT.md), its native JSON summaries and independently certified models. All completed large-gap timing comparisons use 53-bit fields; the new 96-bit result is a certified construction.

Use the existing ten certified D=-43 pairs and two new public targets per endpoint. On each crater compare CM, direct scalar and recurrence; on each floor compare direct scalar and recurrence; add rho on every target. Freeze the exact seeds, 8192-point bases, query caps and balanced endpoint/method order. Bootstrap whole gaps, retaining both fresh targets. Treat ten curve clusters as ten clusters.

## Prerequisites

PR #1888 must land. Implement and independently replay the q=11 recurrence arm and a blinded fresh-target cohort driver. The proposed driver subcommand and independent checker named in card.json are not implemented by this card. Before promotion, implement them, freeze PROTOCOL.md and its blob, pin a merged code commit, add resumable checkpoints, choose the cairn lab and record the actual environment recipe. Follow backlog/README.md for leases and run receipts.

## Expected cost and funding

Planning estimate only: 100 IC attempts and 40 matched rho attempts for ten gaps and two fresh targets per endpoint, allowing several minutes per cold IC attempt, replay and setup. Six core-hours and four wall-hours are provisional estimates, not bounds; the driver has a four-hour wall cap.

These planning values use [session summaries](../../research/conductor_gap_cm_orbits_20261015/results/session-summaries.json), [replication metrics](../../research/conductor_gap_cm_orbits_20261015/results/replication-metrics.json) and, for the upper construction, [the timed receipt](../../research/conductor_gap_cm_orbits_20261015/results/v4-96-time.txt). Every extrapolation is provisional and must be replaced by a measured probe where stated.

Proposed priority: 200000 funding units. Below the existing 250000-unit frozen-cell card: this is a narrower, cheaper confirmation with an existing ten-pair baseline. Review and merge approve the proposal; this card performs no remote funding action.

## What done means

Recoveries and every retained relation pass independent replay. Same-endpoint scalar/recurrence and crater CM/scalar records agree exactly. Report queries, cold CPU and phases, with all failures/censors. Replicated query yield would motivate more gaps; parity would retire the pilot timing lead. Neither outcome alone is an intrinsic-hardness theorem.

Land source, frozen inputs, every result/censor, independent replay, resource receipts and a source-linked report in the card's research destination. Stop at the recorded cap and preserve the checkpoint or bounded failure. Do not reinterpret an unfinished measurement as a successful solve.
