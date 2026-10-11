# Tree-split cut-out screen, round 2: results

[PREREGISTRATION.md](PREREGISTRATION.md) is unchanged since registration commit `8a7109c4`,
which was pushed before any registered cell ran. The rule was also posted to cairn before
the run.

- Run report, graph and PDF: [reports/run-20261011T2224Z/](reports/run-20261011T2224Z/).
- Records: [runs/registered/](runs/registered/), run commit `6506a1e2`, binary
  `28fdd513…cfdb9`.
- Lean: [lean/](lean/).
- cairn log: [cairn/log.jsonl](cairn/log.jsonl).

## Registered verdict: gate **PASS**; all six groups **CLOSED**

| group | cells (uncensored) | gaps `D_s(curve) − min(controls)` | slopes | verdict |
|:--|--:|:--|:--|:--|
| G1 `T2`/`frob2`, n = 21 | 16 (16) | all 0 | 0.6–1.5 | **CLOSED** |
| G2 `T2`/`nbw`, n = 17–23 | 40 (40) | 39 × 0, 1 × −1 | 0.8–1.1 | **CLOSED** |
| G3 `T2`/`poly`, n = 23, 25 | 20 (20) | all 0 | 0.6–0.8 | **CLOSED** |
| G4 `Z`/`poly` | 8 (7) | all 0 | no clause | **CLOSED** |
| G5 `Z`/`nbw` | 8 (7) | all 0 | no clause | **CLOSED** |
| G6 `Z`/`frob2`, n = 15 | 4 (4) | all 0 | no clause | **CLOSED** |
| `frob2`, n = 35 (supplementary) | 4 (4) | all 0 | — | not in the verdict |

- **Run status.** 102 of 102 cells completed, none failed.
- **Censoring.** Two cells were censored, as registered:
  - `Z`/`nbw` with n = 15, a = 0, ℓ = 5, where `|Z|` = 111,628;
  - `Z`/`poly` with n = 15, a = 1, ℓ = 5, where `|Z|` = 132,870.

  For each, degree 5 needs 174,437 monomials, which is over the 90,000 cap.
- **The single −1** is `nbw` with n = 21, a = 1, ℓ = 4, a degenerate cell with
  `|F| = |T2| = 4`:
  - the curve set reads `D_s = 1`;
  - its random 4-point control drew one of its 2,000 samples from its 8-point affine
    closure, a chance of about 0.4%, and so reads 2.

  That is noise in the control, not structure in `T2`.
- **What closes.** Removing every subfield-contained `V` (`frob2`) removes all six of
  round 1's −1 cells. Every `frob2` gap is 0, so round 1's post-hoc reading (subfield
  leakage) holds.

**What this closes.** Subfield-free Frobenius kernels, normal-basis windows and the
polynomial window up to n = 25 give no low-degree description of `T2`. Neither does `Z` on
structured `V`, at the sizes run. Together with round 1 and the 2026-09-30 study, no
factor-base family tried so far gives `T2` or `Z` a low-degree shortcut.

## Lean and cairn

- **`Verdict.lean`, kernel `decide`.** It proves
  `report cells = ⟨gate := true, all six := closed⟩` over the 102 recorded cells, using the
  registered `Rule.lean`, sha256 `ec513363…b1e3`. Changing `t2Nbw` to `.alive` fails the
  build.
- **cairn.**
  - The registration objective `sha256:49b79374…` was accepted before the run.
  - The result objective `sha256:84dde9a6…` was accepted by the kernel.
  - `CAIRN_EPOCH_SECONDS=5 cairn --log cairn/log.jsonl audit` reports the chain intact.
- **No independent replay this round, as registered.** Round 1's `Replay.lean`
  cross-checked the shared code.

## Scope

- **Covered.** Koblitz curves; `T2` at n ≤ 25 (35 supplementary) and ℓ ≤ 8; `Z` at
  n ≤ 15 and ℓ ≤ 5. The metric is a sampled screen.
- **Not covered.** Larger `Z`, other bases, and solving-degree audits.
- **Claims.** No exponent, speedup or deployed-curve claim.
- **Canonical graphs** were checked and left unchanged: the scoreboard, the progress
  timeline and the leaderboard. This study prices no attack.
