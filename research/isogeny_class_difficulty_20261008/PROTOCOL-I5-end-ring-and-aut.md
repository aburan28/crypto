# Protocol I-5: the endomorphism ring, vertical large-degree edges, and automorphisms

Frozen 2026-10-08, before any instrument is built.  **Control experiment;
boundary.**  `S`, end-to-end cost and speedup are **unset**.

## Derivation (stated before measuring)

Within a class, the only group-theoretic quantity that changes is
`End(E)`: horizontal edges keep it, vertical edges (`ℓ | f_π`) move
between orders of conductor differing by `ℓ`.  Two consequences are
testable:

1. **Automorphisms.**  `Aut(E) ⊋ {±1}` only at `j = 0` (`|Aut| = 6`) and
   `j = 1728` (`|Aut| = 4`) in characteristic `> 3`, i.e. only on the
   surface of a `D = −3` or `D = −4` volcano.  Rho with `Aut` folding
   gains exactly `√(|Aut|/2)`; a factor base folds by `|Aut|/2`.  Nowhere
   else in the class does `End(E)` give an automorphism, so nowhere else
   can `End(E)` change a generic or a folded cost.
2. **Vertical large `ℓ`.**  When `f_π` has a large prime factor `ℓ`, the
   `ℓ`-volcano has a surface and a floor one vertical step apart, and that
   step is the one large-prime-degree edge that no small-degree route
   reproduces (it changes the order).  The companion thread's C7 shows no
   endomorphism-based shortcut exists for it on the floor.  Prediction:
   nothing in the DLP's cost depends on which side of that step the curve
   sits, because `E[r]` and `Aut` are the same on both sides when
   `j ∉ {0, 1728}`.

A reader who expects "deeper in the volcano is weaker" (a recurring
suggestion) is making a claim about `End(E)` that this protocol is
designed to falsify at toy size: the only known DLP relevance of
`End(E)` is through the automorphisms it contains and the endomorphisms
usable for folding, both of which this protocol enumerates.

## Instrument (Rust, follow-on PR)

- CM construction of `D = −3` and `D = −4` classes at `p ≈ 2^{20}` and
  `2^{24}` (prime-field CM for class number 1 is a root of `x³ + ax + b`
  with `a = 0` or `b = 0`), plus a random class with `f_π` divisible by a
  prime `ℓ ∈ [61, 2^{12}]`, found by trial over random `p` and traces.
- `isogeny_walk` for the horizontal part; the vertical large-`ℓ` step from
  an `ℓ`-torsion point over `F_{p^k}` with `k` the eigenvalue order (the
  companion's C4 window), certified by the existing edge verifier.
- `ecbench` rho arms with and without `Aut` folding on the special-`j`
  nodes; the I-4 IC arm with endomorphism folding on the same nodes.

## Frozen inputs

| item | value |
|:--|:--|
| CM classes | `D = −3` and `D = −4` at `p ≈ 2^{20}`, `2^{24}`; registered in the follow-on PR |
| vertical class | one class per size with `ℓ | f_π`, `61 < ℓ < 2^{12}` |
| arms | surface node, floor node, a horizontal neighbour of each, the special-`j` node where present |
| targets per arm | 24 at `2^{20}`, 12 at `2^{24}` |
| seed | 20261015 |

## Predictions (pass/fail)

- **Q1 (Aut factor).**  On the `j = 0` node, rho with folding costs
  `S(E)/√3` within the A/A interval; on `j = 1728`, `S(E)/√2`; on every
  other node the folded and unfolded arms agree.
- **Q2 (vertical step).**  Across the vertical large-`ℓ` edge, rho `S`,
  IC yield `γ` and `S_IC` agree within the A/A interval on both sides.
- **Q3 (certificate).**  The vertical edge's kernel polynomial is verified
  by the existing division-polynomial check and the image points have
  order `r`; zero failures.
- **Q4 (folding).**  Factor-base folding by `Aut` on the special-`j` node
  reduces the factor base by exactly `|Aut|/2` and leaves `γ` per
  *class-of-points* unchanged.

## Decision rule (registered)

- Q1–Q4 pass: **boundary (control confirmed)**.  `End(E)` matters only
  through `Aut`, the special `j` are the whole story, and vertical large
  `ℓ` is a reach question with no difficulty payoff.  The methodology
  records `|Aut|` as the one `End`-derived screen.
- Q2 fails: trace the cells under the I-4 control sequence; an effect that
  survives is a **reproducible unexplained anomaly** attached to the
  conductor, and gets its own protocol.  It is not classed higher here.
- Q1 or Q4 off by a constant other than the predicted one: an accounting
  error in the folding; fix and rerun before anything is read.

## Stop condition and inadmissible moves

Bounded: four CM classes, two vertical classes, the arm set above, one
run.

Inadmissible: reading a special-`j` gain as anything but `Aut`; a
vertical edge without a verified kernel; comparing surface and floor
with different `r` (the group order is the same, but the subgroup chosen
must be too); any statement about P-256, whose `D_π` is fundamental and
whose class has no vertical edge at all.
