# Certify fast point evaluation and dual maps on the new 96-bit pair

Status: blocked. This card is a proposal, not a prepared or completed run. No cairn task, lease or objective was created.

## Question and decision

Can the certified 96-bit codomain construction be extended to independently verified point transport and dual composition?

A verified transport permits comparisons on the same hidden scalar and prices the actual bridge between strata. A codomain without point transport cannot support that end-to-end comparison.

## Baseline and protocol to freeze

Source: [the completed study](../../research/conductor_gap_cm_orbits_20261015/REPORT.md), its native JSON summaries and independently certified models. All completed large-gap timing comparisons use 53-bit fields; the new 96-bit result is a certified construction.

First compare point maps against direct Velu on tiny controls and a frozen subset of the 53-bit references. Then use the certified p=78125000378223625457769953627 pair and ell=125000000302579. Transport G and a public unknown-log Q. Construct the dual kernel and check composition on several independent points. Charge kernel preparation and verification separately and preserve all caps.

## Prerequisites

PR #1888 must land. Implement square-root Velu point evaluation and the independent dual-map checker; the current constructor outputs codomains only. The proposed driver subcommand and independent checker named in card.json are not implemented by this card. Before promotion, implement them, freeze PROTOCOL.md and its blob, pin a merged code commit, add resumable checkpoints, choose the cairn lab and record the actual environment recipe. Follow backlog/README.md for leases and run receipts.

## Expected cost and funding

Planning estimate only: the recorded codomain process took 2697.60 wall seconds and 9362772 KiB maximum RSS. Point and dual evaluation have not been measured. Eight core/wall hours is a provisional allowance for small-reference checks and sequential large-map evaluations, not a bound.

These planning values use [session summaries](../../research/conductor_gap_cm_orbits_20261015/results/session-summaries.json), [replication metrics](../../research/conductor_gap_cm_orbits_20261015/results/replication-metrics.json) and, for the upper construction, [the timed receipt](../../research/conductor_gap_cm_orbits_20261015/results/v4-96-time.txt). Every extrapolation is provisional and must be replaced by a measured probe where stated.

Proposed priority: 131072 funding units. At the lower seed-card funding level: the largest-gap endpoints exist, but point evaluation remains new implementation work. Review and merge approve the proposal; this card performs no remote funding action.

## What done means

Independent arithmetic verifies on-curve images, homomorphism controls and dual composition [ell]P. Reject corrupted coefficients/images. Keep codomain-only and point-map claims distinct. Publish cold/warm setup and evaluation receipts; do not claim a 96-bit ECDLP solve.

Land source, frozen inputs, every result/censor, independent replay, resource receipts and a source-linked report in the card's research destination. Stop at the recorded cap and preserve the checkpoint or bounded failure. Do not reinterpret an unfinished measurement as a successful solve.
