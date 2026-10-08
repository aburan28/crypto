# Protocol: the general lopsided regime and the rest of Section 6 (round 4)

Frozen 2026-10-06. **Audit and analysis round**: no performance claim, no
`S` accounting, no `ecbench` sessions, no dashboard edits. This protocol
closes round 2's open check 3 ("the paper's general regime
(epsilon<0.1204, size-dependent gamma) is not encoded in the gate") and
works through the remaining Section 6 directions for transfer targets the
earlier rounds did not cover.

Paper: Alman-Vassilevska Williams, `arXiv:2610.06783v1` (local digest
`sha256:d91a79c5…cc9e76005`; cached at `/tmp/lopsided.pdf`).
Rounds 1–3: thin-product interface + four insertion points (merged PR
#1485); twelve-target survey with the screening gate (open PR #1493);
refutation audit with one sentence repaired (open PR #1494).

Conductor note: no `conductor` binary exists in this environment, so no
scope reservation could be registered (as in rounds 2–3).

## Questions

1. Does the general envelope (Theorems 24–25, Table 2, paper p.40: every
   eps < eps* = 0.1204 admits some gamma > 0) admit any of the 16
   already-screened targets that the concrete 1/18 line excludes?
2. Does the technique transfer to sparse matrix multiplication — i.e. to
   our sparse linear algebra (block Wiedemann over `Z/rZ`)?
3. Do the Section 6 reduction-loss remarks (reductions designed for
   hardness keep only a fraction of the saving; SVX27 black-box
   optimality) or the gray-box problems change any screening verdict?

## Scope and method

- Freeze Table 2's Theorem-25 side as `table2.json`, each cell verified
  against a rendered image of paper p.40 (one extraction error corrected
  in provenance).
- Re-screen all 16 targets (4 insertion points + 12 survey candidates)
  against the *widest* proven envelope (D <= N^0.114 row, kappa up to 1).
  Record per target which gate blocks it and whether the envelope moves
  that gate. Gates for shape (shared thin integer middle with known sparse
  W), ring (small-integer products over Z), and scope are orthogonal to
  the (eps, kappa, gamma) envelope by construction; size is the only gate
  the envelope can move.
- Confirmatory greps: `girth`, sparse-MM hardness claims, balanced-sparse-
  triangle vocabulary outside our own study files.
- No `src/` changes on this branch: the gate-parameterization this
  envelope enables is recorded as a precise next action for after PR
  #1493 merges, with `table2.json` as its frozen input. No dependency on
  unmerged branches is taken.

## Reference and boundaries

Unchanged: generic-group floor (Shoup 1997, unconditional) + matched rho
in unit `S`. Pre-registered expectation: zero of the 16 targets moves,
because every recorded block sits at shape, ring, or scope — gates the
envelope does not touch. A target whose *only* block is size under the
concrete line but which passes under the general envelope would overturn
the headline and open round 2's promotion protocol.

## Success and stop conditions

- Success: envelope frozen with provenance; all 16 targets re-screened
  with the blocking gate named; sparse-MM and reduction-loss questions
  answered with paper citations; method implications recorded as
  proposals with exact next actions.
- Stop: protocol, dated report with the re-screen table, envelope figure
  plus sources, and PDF committed; PR opened. `S`/speedup `null`.

## Deliverables

1. This protocol (`README.md`), the frozen envelope (`table2.json`).
2. The dated analysis (`REPORT.md`).
3. `figure.mmd` (editable source) + `figure.svg`: the proven (eps, gamma)
   envelope at fixed kappa with the concrete line, the eps* ceiling, and
   the beyond-ceiling region marked as needing a new base identity.
4. `report.html` + `report.pdf` (report with visual included).
