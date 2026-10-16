# Protocol: wider implications of the lopsided refutation for this program

Frozen 2026-10-06. **Audit and analysis round**: no performance claim, no
`S` accounting, no `ecbench` sessions, no dashboard edits. This protocol
authorises (a) a systematic audit of what in this repository depends on the
refuted hypotheses, (b) minimal precision repairs where language is now
wrong, and (c) recorded method implications. `S` and speedup remain `null`.

Paper: Alman–Vassilevska Williams, `arXiv:2610.06783v1` (local digest
`sha256:d91a79c5…cc9e76005`; see the
[round-1 protocol](../lopsided_thin_product_20261006/README.md)).
Rounds 1–2: thin-product interface + four insertion points (merged PR
#1485); twelve-target survey with the screening gate (open PR #1493).

Conductor note: no `conductor` binary exists in this environment, so no
scope reservation could be registered (as in round 2). Edited paths are a
new study directory, one research-note sentence, and one encoding repair of
a round-1 figure; `git status` shows no other in-flight work on them.

## Questions

1. Which statements in this repository rest on the refuted 3SUM/APSP/Exact
   Triangle hypotheses, and which rest on unconditional ground?
2. Does any measured verdict, floor, or reference change?
3. What does the discovery method (machine-found identity, Lean-verified)
   imply for how this program should handle algebraic-identity claims?

## Scope and method

- Grep the whole tree for the refuted-hypothesis vocabulary
  (`3SUM`, `APSP`, `Exact Triangle`, `Zero-Weight`, `hinted OMv`,
  fine-grained `SETH`/`OV`) plus the neighbouring false friends
  (`matrix-vector`, `Sethi-Ullman`, `b3sum`, `conditional`), and classify
  every hit with quoted context as: refuted-dependence, already-correct,
  historical measurement, or false positive. Commands and counts are
  receipts in the report.
- One minimal precision repair is pre-authorised: the RR solver panel's
  "3SUM hardness is a conjecture" sentence, which names the conjecture
  itself. The repair preserves the note's conclusion and adds the
  refutation as a sharpening. No frozen variant folder, manifest, or
  measurement is touched.
- One mechanical repair is pre-authorised: the round-1 `figure.svg`
  encoding (legacy bytes and control characters from authoring) to clean
  UTF-8 with identical rendered text.

Out of scope: re-measuring anything, re-matching rho, touching dashboards,
the registry, or frozen evidence.

## Reference and boundaries

Unchanged: generic-group floor (Shoup 1997, unconditional) + matched rho in
unit `S`. The audit's expected result — stated here before running it — is
that no floor, reference, or verdict moves, because the repository never
cited the refuted hypotheses as load-bearing. A hit that contradicts this
expectation would be the finding of the round.

## Success and stop conditions

- Success: every hit classified with quoted context; any required precision
  repair landed; method implications recorded as proposals with exact
  next actions.
- Stop: protocol, dated report with the audit table, dependency diagram
  plus sources, and PDF committed; PR opened. `S`/speedup `null`.

## Deliverables

1. This protocol (`README.md`), the dated audit (`REPORT.md`).
2. `figure.mmd` (editable dependency-graph source) + `figure.svg`.
3. `report.html` + `report.pdf` (report with visual included).
4. The one-sentence RR-panel precision repair + the `figure.svg` encoding
   repair, disclosed in the report and the PR body.
