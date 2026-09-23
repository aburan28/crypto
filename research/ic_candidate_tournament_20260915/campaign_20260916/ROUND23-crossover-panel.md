# Round 0023 — the crossover cells enter the tournament

Pre-registered before any round-0023 stage ran. Nothing below reads round 0023.

## 1. Why this round exists

The tournament has only ever been scored where the collector wins. Its
eight-cell panel tops out at `r` = 4.2·10⁶, and rounds 0019 and 0020 hold
`beats_rho_strict` there under two seeds. Round 0022 measured what lies past
it, in ad-hoc scripts rather than in the harness: the collector crosses rho
between `n23a1` and `n37a0` and loses at `n43a1`, and **the harness as built
would have measured that loss wrong in two ways at once**:

* every case pins the factor base at `points = 6·degree`, one batch of eight
  orbits. That is optimal on the panel and nowhere past it — the best base is
  16 orbits at `n37a0` and 24 at `n43a1`, and the pinned batch costs the
  candidate 2.6× at `n43a1`;
* every arm carries `max_trials = 4096`, which cuts rho off on some `n37a0`
  fixtures. The harness rejects an unfinished trial as unverified, so one
  such fixture disqualifies the comparison outright.

So "keep running tournaments until rho is decisively beaten" has a precondition
the harness did not meet: it has to be able to score the cells where rho is
not beaten. This round puts them on the panel, fixes both defects, and fields
the one candidate round 0022 measured to move them.

**This round is not expected to beat rho, and says so here rather than after.**
Under objective `rho` a challenger is promoted only if every cell is strictly
below rho, and at the four crossover cells no arm the campaign has measured is.
The round's value is the audited, replayed, per-cell measurement of how far
the collector trails rho past the crossing once it stops being handicapped —
and whether the base rule holds outside the range it was fitted on.

## 2. The round

| | |
|:--|:--|
| seed | `2026092323` |
| cells | `--cells 13a0,17a1,19a0,23a0,23a1,31a0,37a0,43a1` |
| holdouts | `--holdout-cells 19a1,29a1,59a0,61a1` |
| allocation | `--confirmation-cases n23a1=40` (carried from round 0019); 172 confirmation cases |
| profile | pilot, one target, objective `rho`, `--require-native-progress` |
| comparison | `--comparison-kind factor-base-policy` — support varies across arms, fixed within an arm |
| limits | `--timeout 300` (worst measured trial 11.6 s under valgrind), `--max-processes 5000` |
| arms | `incumbent`, `scaled` (`round-0023-single-candidates.json`, built by `round23_candidates.py`) |
| trials | about 3,650 |

**Arms.** `incumbent` is round 0020's promoted `both` plus
`round21-wide-pair-table.patch`, with `max_trials` raised to 65,536; its base
is the case's single batch. `scaled` adds `round23-scaled-base.patch`: the
orbit count is `8·round((r/4,196,903)^0.157)`, never below one batch. That rule
reproduces round 0022's measured optima — 8, 16, 24 orbits at `n23a1`, `n37a0`,
`n43a1` — and puts one batch on every panel cell.

**The raised cap is a protocol change, and it is carried by both arms and by
rho.** Round 0022 measured it at ≤ 14 ppm on either arm, both signs, with
identical logarithms and factor bases at every cell.

**Holdouts, and why these two.** `n59a0` (`r` = 1.0·10¹⁰) and `n61a1`
(1.2·10¹⁰) are the certified cells past `n43a1` where rho finishes at 65,536:
8/8 on a pre-round probe, where `n47a1` finished 5/8 and is excluded — one
unfinished rho trial disqualifies a comparison. They lie outside the range the
base rule was fitted on, so they are the round's out-of-sample test of it.

**The arms' own tests were run before this pre-registration was committed.**
`cargo test --lib cryptanalysis::koblitz_tiny_ic` on the `scaled` tree, which
contains every line of the incumbent's; the result is recorded in §5. The
harness still does not run it — the round-0022 next-proposal carries that as
open — so this round does it by hand.

## 3. Predictions

Point values below are round 0022's 64-fixture ladder. The tournament's
confirmation carries twelve cases at each new cell, so its bands will be wider,
and the predictions are bands.

1. **The policy is inert on the panel.** `scaled`/`incumbent` in
   [0.995, 1.005] in instructions at every one of the eight panel cells, on
   both final stages. Both arms build the same single batch there; anything
   outside that band means the patch changed something besides the base size.
2. **It moves the crossover cells.** `scaled`/`incumbent` at `n37a0` in
   [0.65, 0.95] (ladder 0.792) and at `n43a1` in [0.28, 0.50] (ladder 0.381).
3. **It holds outside its fit.** `scaled`/`incumbent` below 0.70 at both
   `n59a0` and `n61a1`. No point value: the rule puts 24 orbits there, the
   incumbent's single batch was 2.6× off at similar `r`, and nothing more
   specific is known.
4. **The panel win survives the widened build and raised cap.** `scaled`/rho
   below one at all eight panel cells in instructions, and `n23a1` below 0.95.
5. **Past the crossing, rho wins by the ladder's margins.** `scaled`/rho at
   `n37a0` in [1.15, 1.70] (ladder 1.396) and at `n43a1` in [1.8, 2.8]
   (2.207). At the holdouts, above one; the rate `r^[0.13, 0.19]` from `n43a1`
   points to about 2.3–2.6, but the holdouts differ in degree and cofactor
   (`n59a0`'s cofactor is 5.7·10⁷), so only "above one" is registered.
6. **The decision.** `scaled` passes every gate against the incumbent on both
   final stages, fails `rho_gate` at the four crossover cells, and is not
   promoted: status `retained`, `beats_rho_strict` false. **This is the
   expected outcome, not a failure of the round**, and it is written down so
   that it cannot be re-read afterwards as one.

## 4. Falsification and scope

Any bad, missing or unmatched certificate; any unfinished trial in either arm
or in rho; any recovered logarithm that differs between the arms on the same
fixture. Prediction 1 failing is a defect in the patch and voids the
crossover-cell comparisons until it is explained.

Every conclusion is scoped to these twelve cells, one seed, the pilot profile
and this collector. Two of the four crossover cells are in the sample the base
rule was fitted to, so their `scaled`/`incumbent` ratios are in-sample for the
rule, though out-of-sample for fixtures and seed. Only `n59a0` and `n61a1`
test it. Nothing here bears on sect163k1 or on any curve outside the census.

## 5. Arm tests (filled in before the first stage runs)

