# The depth-1 refusal wired into the §17 walk

Run 2026-10-08 on branch `research/jv-cover-end-to-end-20261007` with
`--refuse-depth1` added to `examples/jv_isogeny_walk.rs` and `run_walk2`
(`src/cryptanalysis/jv_isogeny_walk.rs`).  Same seed, trials, cap, jumps
and samples as §17.5's grid (`--v2 --sizes 7,11,13,17,23,31,53 --trials 40
--cap-mult 3 --jumps 3,5,7 --samples 20000 --seed 1`), one run with the
refusal and one without.  **Class: engineering for the walk** (the test
is a congruence on the trace the walk already has; no algorithm changed
on the admitted trials).  Nothing about the route's `S / rho` changes.

## The test

`t = q³ + 1 − #E` is read from the point count (or the public order in
the end-to-end setting); `D = t² − 4q³`; the class's 2-volcano has depth 1
iff `v₂(−D) = 5` (`t = 2t′`, `t′` odd, `q³ ≡ 1 mod 8`, depth 1 iff
`v₂(t′² − q³) = 3`).  Such a class holds no weak curve (PR #1556: 3,690
classes, 4,871,062 curves, no exception, `p ≤ 19`), so the walk refuses it
before enumerating anything.  A refusal is charged its point count.

## Results

| p | trials | refused | start-weak | found | success, as §17 counts it | success on admitted trials | muls per walk, all (refusing) | muls per walk, all (baseline) | baseline / refusing | muls per admitted walk | muls per refusal | §17.5's success |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 7 | 40 | 12 | 1 | 22 | 0.54 | 0.78 | 9.53e5 | 4.28e+06 | 4.5× | 1.34e6 | 6.1e4 | 0.54 |
| 11 | 40 | 17 | 1 | 15 | 0.36 | 0.64 | 1.78e6 | 6.02e+06 | 3.4× | 3.02e6 | 8.4e4 | 0.36 |
| 13 | 40 | 16 | 1 | 17 | 0.41 | 0.70 | 1.46e6 | 1.89e+07 | 12.9× | 2.37e6 | 1.0e5 | 0.41 |
| 17 | 40 | 19 | 0 | 13 | 0.33 | 0.62 | 1.26e6 | 2.94e+07 | 23.3× | 2.28e6 | 1.3e5 | 0.33 |
| 23 | 40 | 14 | 0 | 17 | 0.42 | 0.65 | 1.26e6 | 1.36e+08 | 108.4× | 1.84e6 | 1.8e5 | 0.42 |
| 31 | 40 | 15 | 0 | 15 | 0.38 | 0.60 | 7.60e6 | 3.94e+08 | 51.9× | 1.20e7 | 2.6e5 | 0.38 |
| 53 | 40 | 16 | 0 | 10 | 0.25 | 0.42 | 1.63e7 | 2.57e+08 | 15.8× | 2.68e7 | 4.5e5 | 0.25 |

Reading:

- The refusing run reproduces §17.5's success rates exactly at every
  size, as it must: the refused trials are exactly trials that §17's
  walk ran to exhaustion or cap, and the admitted trials are unchanged.
- Between 30% and 48% of random full-2-torsion starts are refused, at
  2% to 6% of the cost of an admitted walk (the point count), where the
  baseline spent a full exhaustion or cap on them.
- The mean cost per trial, every trial inside, falls by `4.5×` (`p = 7`), `3.4×` (`11`), `12.9×` (`13`), `23×` (`17`) and `109×` (`23`): an exhausted walk enumerates the start's whole reachable closure, and the refusal replaces that with one point count.  At `p = 31` and `53` the factors are `51.9×` and `15.8×`.
- On admitted trials the walk succeeds in 60% to 78% of starts at
  `p ≤ 31` and 42% at `p = 53`.  What remains is the component-reach
  limit of §17.5 (the 2,3,5,7-reachable set of the start is a small part
  of the class), which the refusal does not touch.
- The ascent prefix of PR #1563's policy A has no analogue here: this
  walk enumerates whole 2-isogeny components, so every height in a
  component is visited anyway.

Files: `refuse.json`, `baseline.json` (the `Walk2Report`s, with
`refused`, `success_fraction_admitted`, `mean_muls_admitted`,
`mean_muls_refused` on the refusing run; every trial row carries
`v2_disc` and `refused`).
