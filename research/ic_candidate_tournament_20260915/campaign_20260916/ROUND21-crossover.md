# Round 0021 — the crossover, measured

> **CORRECTION, round 0022.** §4's `n37a0` reading of **1.533 [1.238, 1.899]
> is upper-biased and superseded**; see `ROUND22-budget-and-curvature.md`.
> `round21_crossover.py` verified that the IC arm completed and never checked
> rho's status before charging it — `instructions()` reads callgrind's
> `Collected:` line, not the worker's JSON — so **4 of the 64 rho runs at
> `n37a0` were cut off at `max_trials = 4096`** and their truncated counts went
> into the denominator. An under-charged denominator inflates IC/rho, so the
> correction runs **downward**. `n23a1`'s 0.831 is unaffected: 64 of 64
> completed there. The claim in §4 that "both arms complete at both cells" is
> withdrawn — it was asserted, not checked. Everything below is left as
> published; §§1–3 (the widening, its cost, the red test) are unaffected.

Rounds 0019 and 0020 established a strict win over rho on the eight-cell panel
and replicated it, and both classified it `engineering` for one reason: the
`Θ(r^{2/3})` against `Θ(r^{1/2})` boundary was **derived**, its only measured
support a same-degree contrast between two cells. `PROBE-degree-ceiling.md`
then found the ninth cell that would settle it — and that the promoted
implementation could not run there.

This round removes that obstacle and takes the measurement.

## 1. The ceiling was eight bytes

`koblitz_tiny_ic` declared `MAX_DEGREE = 31`. Nothing about its arithmetic
needed that: the field has always been one `u64`. The ceiling lived in the pair
table, which packed each stored sum's coordinates into `u32` and said so —
*"twelve bytes, since every coordinate fits 31 bits."*

`round21-wide-pair-table.patch` widens those two fields to `u64` and lifts both
degree bounds that mirrored them:

| bound | was | now |
|:--|--:|--:|
| `koblitz_tiny_ic::MAX_DEGREE` | 31 | 61 |
| worker dispatch guard | `(5..=31)` | `(5..=61)` |
| `oracle.py` `Curve.__init__` | `5..31` | `5..61` |

Those are three independent bounds. Lifting the module's alone left the worker
still refusing degree 37, which cost a build to discover — they have to move
together, and `round21_build.sh` now checks all three before it will hand back
a binary.

The oracle's bound was never a width limit — its arithmetic is Python integers
over the fixture's own irreducible polynomial. It exists so the checker refuses
exactly what the collector refuses, so it follows. The amendment is additive:
it accepts strictly more and reads every earlier fixture identically, which
`test_certificate.py::DegreeBoundTests` asserts by *which* failure each degree
produces rather than by whether one does — `Curve.__init__` goes on to check
the group order, cofactor, generator and eigenvalue, so a fixture with only its
degree swapped fails on the group order no matter what the bound says.

## 2. What the widening cost

Identical answers, measured rather than assumed. Every one of the eight panel
cells, three confirmation fixtures each, widened build against the promoted
worker:

**24 of 24 fixtures returned the same logarithms and the same factor-base
hash.** The instruction cost of carrying eight more bytes an entry:

| cell | n13a0 | n17a1 | n19a0 | n19a1 | n23a0 | n23a1 | n29a1 | n31a0 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| widened / promoted | 1.0011 | 1.0008 | 1.0015 | 1.0006 | 1.0003 | 1.0002 | 1.0004 | 1.0004 |

**0.07% panel-wide**, worst cell 0.15%. The table is a few thousand entries at
panel sizes, so doubling an entry's width costs almost nothing.

## 3. A test that had been red since round 0018

Building the widened tree surfaced a failure at `n13a0` — the recovered scalar
coming back exactly 4× expected, and `n13a0`'s cofactor is 4. It reproduces
**identically on the unmodified round-0020 winner**, so it is not this round's.

`complete_solve_verifies_in_general_arithmetic` asserted

```rust
assert_eq!(mulmod(h, (w.a + mulmod(w.b, d, r)) % r, r), logsum, ...)
```

— the row reading `h·(a + b·d)`, which is the **cofactor** convention. Round
0018's `column` patch moved the solver to the **representative** convention,
where a column is the orbit representative itself and the row reads
`a + b·d` with the factor of `h` gone. `oracle.py` was amended for that at the
time. This test was not, and has been failing at every cell whose cofactor is
not one ever since — through rounds 0019 and 0020, both of which promoted the
arm that broke it.

**It does not invalidate those rounds, and the reason is specific.** The
assertion that failed is an internal convention cross-check between the tiny
path and the general one. The assertion in the *same test* that the recovered
logarithm equals the planted scalar — `assert_eq!(found, BigUint::from(d))` —
passes, and passed throughout. Every published result was verified by
`oracle.py` against the planted secret, 2,646 receipts a round, and that
checker had been correctly amended.

What it does expose is a real gap in the process: **the campaign verified
through `oracle.py` and the tournament harness and never ran `cargo test` on a
candidate arm's own source.** A convention change that the checker was taught
about, and the code's own tests were not, could sit red for three rounds.
Fixed here, in the patch, to assert the representative convention.

## 4. The crossover

`n37a0` — degree 37, `r = 230,603,167`, **55× the panel's largest cell** —
against `n23a1` at `r = 4,196,903`. Both arms complete at both cells; every IC
report is checked by `oracle.py`, which is what makes these admissible where
`PROBE-degree-ceiling.md`'s numbers were not (that probe measured the general
fallback path, whose output the checker rejects as not certified).

Sixty-four independent fixtures a cell, pooled over two seed streams:

| cell | r | IC/rho | 95% band | per-case sd(log) |
|:--|--:|--:|:--|--:|
| `n23a1` | 4,196,903 | **0.831** | [0.772, 0.894] | 0.300 |
| `n37a0` | 230,603,167 | **1.533** | [1.238, 1.899] | 0.872 |

**Both bands exclude one, in opposite directions. The collector crosses rho
between these two subgroup orders.** Measured, on cells where both arms
complete and every answer is certified — not extrapolated.

The `n23a1` value is a cross-check as well as an anchor: 0.831 [0.772, 0.894]
here, against 0.7498 and 0.7948 from rounds 0019 and 0020 at forty fixtures
each. The widened build reproduces the published panel.

Against the derived boundary: the ratio grows as **r^0.153** across those two
cells, where `round19_model.py` derives **r^{1/6} = r^0.167** for the balanced
optimum. Two cells fix one rate and no curvature, so that is a rate between two
points and not a fit, and the two cells differ in degree as well as in `r`. The
agreement is worth exactly what a two-point rate is worth — but it is the first
time the boundary has had any measured support at all beyond the panel.

### The sample size was arrived at the hard way

Three fixtures read `n23a1` at 0.882, then 1.097 on a second draw. Sixteen read
`n37a0` at 1.782 [1.394, 2.278] on one stream and **1.116 [0.578, 2.154]** on
another — a band containing one, which is no answer. `n37a0`'s per-case log
spread is 0.872, nearly triple `n23a1`'s 0.300.

That is round 0019's finding — the spread grows with the cell because rho's
collision search is a growing share of its cost — holding two decades of `r`
further out than round 0019 could reach. It is also the exact error rounds
0018b through 0020 exist to have caught, and this round nearly published it a
second time. Sixty-four fixtures bring both standard errors under 0.11 and both
bands clear of one.

## 5. What is not claimed

**`n41a0` yields no comparison.** At `r = 5.5·10¹¹` IC completes with a large
enough base (104 orbits, 3.8s) and rho does not complete at all — returning in
0.06s, far too fast for the ~72,600 steps it would need. That is the frozen
trial budget cutting rho off, not rho being out-run. An arm that does not
complete has no cost, and a budget exhaustion is never evidence about an
algorithm. Reported here only so the absence is on the record.

**The factor base was not tuned toward the answer.** The base is swept at
`n37a0` and every size reported: 8 orbits 1.533, 32 orbits 1.815, 40 orbits
1.884, 56 orbits 2.105 on the sixteen-fixture pass. Bigger is monotonically
worse, which is what round 0019's sweep found at all eight panel cells,
reproduced two decades of `r` higher. The cheapest configuration is the one the
sampler builds by default, and that is the one the headline uses.

**A prediction of mine was wrong and is withdrawn.** Before writing any code I
computed that `n37a0` would need 32 orbits and that 8 could not solve it inside
`max_trials`, and I let that shape the design. Eight orbits solve it and are
the *best* configuration. The error: `max_trials` bounds relation-producing
trials, while the quantity in the cost model is the pair-table scan count —
two different things. The same mistake produced the claim that `n41a0` "needs
~39,922 orbits"; the empirical fact there is only that rho does not complete.

**This does not overturn rounds 0019 and 0020.** Their strict win stands
exactly as scoped: eight cells, subgroup orders 2·10³ to 4·10⁶, two seeds. This
round adds the thing that scope was always waiting for — a measured cell beyond
it, where the same collector loses. By the AGENTS.md §3 test the campaign's
gains remain **engineering**, and there is now a measurement rather than a
derivation behind the reason why.
