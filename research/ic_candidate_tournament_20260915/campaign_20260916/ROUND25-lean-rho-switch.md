# Round 0025 — a lean rho, and a collector that picks pairs or triples

Pre-registered before any round-0025 stage ran. Everything here that reads
measurements reads round 0024, which is closed, or the pre-round probe in §4,
which draws no fixture from round 0025's seed.

## 1. Why this round exists

Round 0024 ([RESULTS.md](RESULTS.md#round-0024-a-matched-rho-and-the-triple-table))
measured the pair and triple collectors against a *matched* rho, one that
names Frobenius classes the way the IC arm does. It left three directions
open, and this round is where two of them are scored. The third is §6.

* **A leaner rho.** The matched rho still ran its walk in the library's
  general arithmetic, with Fermat inversion and double-and-add setup, while
  the IC arm used its own single-word arithmetic.
  [round25-rho-lean-notes.md](round25-rho-lean-notes.md) moves rho onto the IC
  arm's own `Arith`. The iteration function, class names, sign rule, cycle
  escape, restarts, draw order and report are unchanged. On 72 of round 0024's
  fixtures the walk is **identical step for step**, with the same iterations,
  walk additions, restarts and logarithm, at 0.28–0.58× the rho-phase
  instructions.
* **A pair/triple switch.** [round25-switch-notes.md](round25-switch-notes.md)
  prices both collectors with their own sizing models (W₂, W₃) before any base
  point is drawn, and runs the triple table if and only if W₃ < W₂. Doing so
  needs an additive checker amendment, which lets a report declare fewer
  summands than its arm is configured with. It lives in its own evaluator copy,
  [evaluator-r25](evaluator-r25/).

**This round is not expected to beat rho anywhere, and says so here rather
than after.** §4 measures both IC arms above the lean rho at every cell.

## 2. The round

| | |
|:--|:--|
| seed | `2026092525` |
| cells | `--cells 13a0,17a1,19a0,23a0,23a1,31a0,37a0,43a1` (rounds 0023–0024's) |
| holdouts | `--holdout-cells 19a1,29a1,59a0,61a1` |
| allocation | `--confirmation-cases n23a1=40,n37a0=40,n43a1=40`; 228 confirmation cases |
| profile | pilot, one target, objective `rho`, `--require-native-progress` |
| comparison | `--comparison-kind factor-base-policy` |
| limits | `--timeout 300`, `--max-processes 6000` |
| arms | `incumbent` (pair, `--source-root round25-sources/lean`), `switch` (`pair_or_triple`, `summands: 4`); `round-0025-candidates.json`, trees from `round25_candidates.py` |
| rho | the incumbent build's rho mode: matched walk, lean arithmetic |
| evaluator | `evaluator-r25/tournament.py` (sha256 `1e5b27743bf3…`, unchanged from round 0023's) and `evaluator-r25/oracle.py` (`ae8f9b2a21b5…`, round 0023's plus `round25-oracle-mixed-summands.patch`) |
| trials | 4,656: aa 48, smoke 72, development 216, selection 216, confirmation 2,052, replay 2,052 |

**Why the amended checker is safe to score with.** When a report's declared
summand count equals its configuration, the amended checker is identical to
the old one, down to the certificate and every refusal message. All 34,668
stored profile and native reports of rounds 0014–0020 verify identically under
both ([round25-switch-check.txt](round25-switch-check.txt)). A report may
declare fewer summands than configured, down to three, uniform within the
report, and every relation is still re-added in the group.

## 3. Before the round: tests

[round25-tests.txt](round25-tests.txt). Each tree was built and tested in its
**own** cargo target directory.

* `lean`: `koblitz_tiny_ic` 10/10. `--lib rho` 41/42; the one failure is
  `single_word_rho_walk_matches_the_reference_step_for_step`, which pins the
  old walk's class representative and fails by construction, as on round
  0024's matched tree. The 41 include the lean patch's two walk-identity
  tests.
* `switch`: `koblitz_tiny_ic` 16/16, including the three switch tests.
  `--lib rho` 41/42, with the same single failure.

**Correction to commit `1722bad1`.** That commit and PR #925's first
description reported "switch 16/16" and "lean 16/16" from a build in which
both trees shared one cargo target directory. Cargo treated the second tree
as up to date, so both workers and both test runs were one tree's code. The
switch worker was byte-identical to the lean worker, which is how the probe
caught it. The results above were re-run in separate directories. The lean
worker rebuilt byte-identical to the one the probe used.

## 4. Before the round: the probe

[round25_probe.py](round25_probe.py), raw output in
[round25-probe.txt](round25-probe.txt). It used four of round 0024's closed
fixtures a cell, and checked every report with `evaluator-r25`. **Every arm
finished and all four agreed on every logarithm.** The lean and matched rho
reported identical walks on every fixture. Instructions, geometric means over
the four fixtures:

| cell | switch picks | pair / lean rho | switch / lean rho | switch / pair | lean / matched rho |
|:--|:--|--:|--:|--:|--:|
| `n13a0` | pair | 1.419 | 1.432 | 1.009 | 0.486 |
| `n17a1` | pair | 1.335 | 1.344 | 1.007 | 0.454 |
| `n19a0` | pair | 1.221 | 1.228 | 1.006 | 0.468 |
| `n19a1` | pair | 1.436 | 1.445 | 1.007 | 0.448 |
| `n23a0` | pair | 1.762 | 1.776 | 1.008 | 0.492 |
| `n23a1` | pair | 2.720 | 2.742 | 1.008 | 0.450 |
| `n29a1` | pair | 1.698 | 1.702 | 1.002 | 0.489 |
| `n31a0` | pair | 1.779 | 1.789 | 1.005 | 0.442 |
| `n37a0` | triple | 6.044 | 2.554 | 0.423 | 0.482 |
| `n43a1` | triple | 21.90 | 5.294 | 0.242 | 0.587 |
| `n59a0` | triple | 1.925 | 1.348 | 0.700 | 0.951 |
| `n61a1` | triple | 11.98 | 4.562 | 0.381 | 0.602 |

`n59a0`'s lean/matched ratio is near one because curve construction, which
both rhos share, is most of its job (cofactor 5.7·10⁷).

## 5. Predictions

1. **Correctness.** Every report verifies under `evaluator-r25`. All three
   arms recover the same logarithm on every fixture of every stage. No trial
   is unfinished, and the A/A stage passes.
2. **The switch chooses as its rule says.** It picks the pair collector on all
   eight panel cells and the triple table on all four crossover cells, on
   every fixture. The rule reads only `(n, r, points)`, so this is fixed per
   cell.
3. **The switch costs about nothing where it runs pairs, and saves where it
   runs triples.** `switch`/`incumbent` in instructions is in
   [1.000, 1.015] on every panel cell, and in [0.30, 0.65] at `n37a0`,
   [0.15, 0.40] at `n43a1`, [0.55, 0.85] at `n59a0` and [0.20, 0.55] at
   `n61a1`.
4. **Neither arm beats the lean rho anywhere.** In instructions,
   `incumbent`/rho and `switch`/rho are above one at every cell on both final
   stages. The panel cells fall in [1.1, 3.5], and the crossover cells are
   above 1.2.
5. **The decision.** `beats_rho_strict` is false and `rho_parity` is false. The
   status is `retained` whichever way the incumbent gate goes, because
   `rho_gate` fails at every cell. Whether `switch` passes the incumbent gate
   is not predicted. It saves at four cells and is flat at eight, and the
   cell-bootstrapped interval may or may not clear one.

Native time is reported, not predicted.

**How this round re-reads round 0024, fixed now.** If prediction 4 holds, the
record reads: *against a rho in the IC arm's own arithmetic, no IC arm this
campaign has built is below rho at any of the twelve cells.* That includes the
six small cells where round 0024's pair arm was below the matched rho. It is
recorded as an additive correction scoped to this seed, and no earlier record
is rewritten.

## 6. The third direction: a different decomposition

[DECOMPOSITION-SURVEY.md](DECOMPOSITION-SURVEY.md) closes generic higher-arity
tables as to the exponent. That covers the pair and triple arms, and so this
round's switch too. Only constants remain for them. The one open family is
algebraic decomposition with `m ≥ 4` summands. Its next step is an `m = 4`
exponent audit with a fixed decision rule, pre-registered separately under
`research/ic_m4_exponent_audit_20260928/`. It is not a tournament arm and is
not part of this round.

## 7. Falsification and scope

This round is refuted by any of the following:

* any bad, missing or unmatched certificate;
* any unfinished trial;
* any logarithm that differs between arms on a fixture;
* any switch choice other than prediction 2's;
* `switch`/`incumbent` above 1.02 on a panel cell.

Every conclusion is scoped to these twelve cells, one seed, these two
collectors, this rho and this checker. The lean rho is not offered as the
leanest possible rho. It uses the IC arm's own primitives, so that rho is not
given arithmetic the IC arm lacks.

## 8. Outcome (appended after the round; §1–7 unchanged)

The round ended `retained`, with `beats_rho_strict` and `rho_parity` false:
prediction 5. The audit verified all 4,656 receipts, and all five predictions
are confirmed.

- **Rho.** In instructions, both IC arms are above the lean rho at every
  cell on both stages: 1.36–2.56× on the panel, and at the crossover cells
  2.1–13.7× for the pair arm and 1.35–5.0× for the switch. Round 0024's six
  small-cell wins against the matched rho therefore do not hold against the
  lean rho, as §5 fixed in advance.
- **Switch.** It chose as its rule says on every fixture. It passed the
  incumbent gate, which was not predicted: 0.771 [0.601, 0.957] in
  instructions.
- **Native time**, which was not predicted, puts both IC arms below the lean
  rho at the four smallest cells (0.87–0.98).

Full account:
[RESULTS.md](RESULTS.md#round-0025-a-lean-rho-and-the-pairtriple-switch),
[report](../runs/round-0025/REPORT.md).
