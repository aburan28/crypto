# Transport all four geometric seeds in the coordinate control

Status: blocked. This card is a proposal, not a prepared or completed run. No cairn task, lease or objective was created.

## Question and decision

Do isomorphic models produce identical folded relation work when G, Q and all four geometric seed points are transported?

This removes the known seed-selection confound in the completed coordinate experiment and determines whether any remaining difference belongs to arithmetic or to point-set selection.

## Baseline and protocol to freeze

Source: [the completed study](../../research/conductor_gap_cm_orbits_20261015/REPORT.md), its native JSON summaries and independently certified models. All completed large-gap timing comparisons use 53-bit fields; the new 96-bit result is a certified construction.

Reuse the frozen ten-gap coordinate cohort. Transport G, Q and the four public seed points by the documented isomorphism. Keep scalar coefficients, query seeds and orbit lengths fixed. Normalize signs after pullback before comparing records. Do not rerun prefix seed selection in the new model.

## Prerequisites

PR #1888 must land. Implement explicit public seed inputs and an independent pullback/sign-normalization audit. The proposed driver subcommand and independent checker named in card.json are not implemented by this card. Before promotion, implement them, freeze PROTOCOL.md and its blob, pin a merged code commit, add resumable checkpoints, choose the cairn lab and record the actual environment recipe. Follow backlog/README.md for leases and run receipts.

## Expected cost and funding

Planning estimate only: ten gaps, both endpoints and two coordinate models give 40 IC attempts plus 40 rho controls. Three core-hours is a provisional allowance based on the existing isomorphic cohort, including cold setup and replay; rho remains diagnostic because its coordinate-hash partition changes.

These planning values use [session summaries](../../research/conductor_gap_cm_orbits_20261015/results/session-summaries.json), [replication metrics](../../research/conductor_gap_cm_orbits_20261015/results/replication-metrics.json) and, for the upper construction, [the timed receipt](../../research/conductor_gap_cm_orbits_20261015/results/v4-96-time.txt). Every extrapolation is provisional and must be replaced by a measured probe where stated.

Proposed priority: 160000 funding units. Below D=-43 replication and the existing 250000-unit card; this resolves a specific interpretation gap with a small cohort. Review and merge approve the proposal; this card performs no remote funding action.

## What done means

The independent pullback audit matches all base points, coefficients and retained relations. Reject a deliberately corrupted transport. Compare IC query counts and per-query/cold CPU; preserve every failure. A residual CPU effect isolates representation arithmetic more closely than the completed prefix-reselection control.

Land source, frozen inputs, every result/censor, independent replay, resource receipts and a source-linked report in the card's research destination. Stop at the recorded cap and preserve the checkpoint or bounded failure. Do not reinterpret an unfinished measurement as a successful solve.
