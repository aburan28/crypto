#set page(paper: "a4", margin: (x: 1.8cm, y: 1.5cm))
#set text(size: 9.5pt)
#set par(justify: true)
#show heading: set block(above: 1em, below: 0.5em)

= Leaving a blocked twist tower: 3-isogeny walks and the class ceiling

_2026-10-14. Frozen protocol `PROTOCOL.md` (sha256 `86a67503d4b5c4c5…`). Native instrument `examples/twist_walk_blocked_towers.rs`. Both positive controls passed before any data. Scope: 5 catalog curves, 13 towers, depth 4, degree 3 only._

== Catalog result (exact computation)

Twist of each of the 13 blocked towers; depth counts are 3-isogeny steps, Frobenius free. `m` is the MMT magic number of the start; `m'` is the least magic number among unblocked nodes within depth 4.

#table(
  columns: (auto, auto, auto, auto, auto, auto, 1fr, auto),
  align: (left, right, right, right, right, left, left, right),
  stroke: 0.4pt + luma(180),
  inset: 3pt,
  [*Curve*], [*N*], [*l*], [*n*], [*m*], [*nodes d1,d2,d3,d4*], [*verdict*], [*m'*],
  [icv1-f2m1…d2c04152], [18], [1], [18], [13], [0, 0, 0, 0], [no rational 3-isogeny], [—],
  [icv1-f2m1…d2c04152], [18], [3], [6], [5], [0, 0, 0, 0], [no rational 3-isogeny], [—],
  [icv1-f2m2…54539071], [24], [1], [24], [23], [0, 0, 0, 0], [no rational 3-isogeny], [—],
  [icv1-f2m2…2d0c0e8c], [24], [1], [24], [22], [1, 0, 0, 0], [no unblocked twist within 4 steps], [—],
  [icv1-f2m2…2d0c0e8c], [24], [2], [12], [11], [1, 0, 0, 0], [unblocked at depth 1], [12],
  [c2pnb208w1], [208], [1], [208], [14], [0, 0, 0, 0], [no rational 3-isogeny], [—],
  [c2pnb208w1], [208], [2], [104], [7], [0, 0, 0, 0], [no rational 3-isogeny], [—],
  [c2pnb208w1], [208], [13], [16], [14], [0, 0, 0, 0], [no rational 3-isogeny], [—],
  [c2pnb208w1], [208], [26], [8], [7], [0, 0, 0, 0], [no rational 3-isogeny], [—],
  [c2pnb368w1], [368], [1], [368], [13], [1, 3, 0, 0], [no unblocked twist within 4 steps], [—],
  [c2pnb368w1], [368], [2], [184], [7], [1, 3, 0, 0], [unblocked at depth 1], [8],
  [c2pnb368w1], [368], [23], [16], [13], [1, 3, 0, 0], [no unblocked twist within 4 steps], [—],
  [c2pnb368w1], [368], [46], [8], [7], [1, 3, 0, 0], [unblocked at depth 1], [8],
)

Of 13 towers, 3 are unblocked after one step (`N=24, l=2`; `c2pnb368w1, l=2` and `l=46`). 7 have no rational 3-isogeny at all (three curves, including `c2pnb208w1`). 3 have neighbours and stay blocked: the two `l = 1` towers and `c2pnb368w1, l = 23`.

== Random blocked starts, for context (not tested against the catalog)

2000 random starts per configuration, same field, `Tr_(K\/k)(b) = 0` forced; 95% Clopper–Pearson intervals. The last column counts escapes whose least `m'` equals `n` (the generic maximum).

#table(
  columns: (auto, auto, auto, auto, 1fr, 1fr, auto),
  align: (right, right, right, right, left, left, right),
  stroke: 0.4pt + luma(180),
  inset: 3pt,
  [*N*], [*l*], [*n*], [*starts*], [*has a 3-isogeny*], [*escaped within 4 steps*], [*escapes with m' = n*],
  [18], [1], [18], [2000], [74.2% (0.722–0.761)], [0.0% (0.0000–0.0018)], [0/1],
  [18], [3], [6], [2000], [75.4% (0.735–0.773)], [64.6% (0.6246–0.6670)], [1225/1292],
  [24], [1], [24], [2000], [74.3% (0.723–0.762)], [0.0% (0.0000–0.0018)], [0/1],
  [24], [2], [12], [2000], [73.7% (0.717–0.756)], [36.8% (0.3468–0.3896)], [630/736],
  [208], [1], [208], [2000], [75.2% (0.732–0.770)], [0.0% (0.0000–0.0018)], [0/1],
  [208], [2], [104], [2000], [74.3% (0.723–0.762)], [36.0% (0.3394–0.3820)], [721/721],
  [208], [13], [16], [2000], [75.0% (0.731–0.769)], [75.0% (0.7309–0.7693)], [1501/1501],
  [208], [26], [8], [2000], [74.7% (0.727–0.766)], [74.7% (0.7273–0.7659)], [1494/1494],
  [368], [1], [368], [2000], [75.6% (0.737–0.775)], [0.0% (0.0000–0.0018)], [0/1],
  [368], [2], [184], [2000], [77.1% (0.752–0.789)], [39.9% (0.3775–0.4208)], [798/798],
  [368], [23], [16], [2000], [74.8% (0.728–0.767)], [74.8% (0.7284–0.7669)], [1496/1496],
  [368], [46], [8], [2000], [74.4% (0.724–0.763)], [74.4% (0.7243–0.7630)], [1488/1488],
)

#figure(image("escape_by_config.svg", width: 88%), caption: [Escape rate of random blocked starts within 4 steps (blue, 95% interval) and the catalog towers (orange). Source: `summary.csv`.])

== Class ceiling (exact computation, finite)

The twist-side isogeny class of `b` is `{b' : t(b') = t(b)}`. For every field `N = 8…22` and every `l | N` with `n = N/l` even we counted the trace classes holding a blocked curve and those that also hold an unblocked one (`class_ceiling.csv`). For `l = 1` no class is mixed in any of the 8 field sizes (0 mixed classes in total), so no isogeny of any degree leaves an `l = 1` blocked tower in those fields. For `l >= 2`, 60% to 99.9% of blocked curves lie in a mixed class, so an escape exists inside the class.

#figure(image("class_ceiling.svg", width: 88%), caption: [Share of blocked curves whose isogeny class also contains an unblocked curve, by `l`. Source: `class_ceiling.csv`.])

== Claims and scope

- *Exact computation:* the 13 catalog outcomes above; the class counts for `N = 8…22`.
- *Observed on a sample:* 0 of 2000 random starts escape an `l = 1` tower in each of four fields (upper limit 0.0018), as the class counts require; 36–40% escape at `l = 2`, 65% at `N = 18, l = 3`, and 74–75% at `l >= 13`.
- *Not tested:* that the catalog differs from random starts. Five independent graphs cannot support that test at the corrected level (`PROTOCOL.md`). One catalog escape is far from generic: `c2pnb368w1` at `l = 2` reaches `m' = 8` where random escapes all have `m' = n = 184`; the other two escapes have `m' = n`. That contrast is descriptive.
- *Conjecture, not proved here:* `Tr_(K\/F_2)(b)` is constant on each trace class for all even `N`; the exhaustive counts support it for `N <= 22`.
- *Not shown:* that any twist is attackable. An unblocked twist tower removes the Lemma 6 obstruction only; the twist's group order and cofactor, the cover's genus and the descent cost are separate. Degrees other than 3 are not walked, and depth is capped at 4.
