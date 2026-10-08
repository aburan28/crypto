# Round 0019 — the boundary the method is running into, and the gate that cannot see it

Pre-registered before any round-0019 stage ran. Every number below is read
from **frozen rounds 0017 and 0018b** — 24 independent fixtures a cell, under
two independent seeds, for the two executables this round runs — or derived
from the solver's own published cost model. Nothing here reads round 0019.

Baseline (`incumbent`): `runs/round-0018b/source`, worker sha256
`e9f263b839d439424f0d…`, byte-identical to `runs/round-0017/source_candidates/orbits/worker`,
the executable rounds 0017 and 0018b both measured.

---

## 1. What round 0018b actually established, and what it did not

Round 0018b found `beats_rho_strict = false` where round 0017, three days
earlier, had found it true **for the same executable on the same eight-cell
panel**. One cell moved: `n23a1`, where the winner/rho instruction ratio read
0.7765 under seed 2026091717 and 1.0231 under seed 2026091818.

Round 0018b recorded that as a replication failure and stopped there. It is
worth one more step, because the two readings are not in conflict. Pooling
both rounds' 24 independent fixtures at `n23a1` gives a winner/rho ratio of
**0.8913, 95% band [0.7802, 1.0182]** — a single estimate that contains both
readings comfortably. The per-case log spread of that ratio at `n23a1` is
**0.329**. Over twelve fixtures that is a standard error of 0.095, and a true
ratio of 0.89 lands above 1.0 about **9% of the time**. Round 0018b spent that
9%; round 0017 did not.

So the campaign's central claim was never as strong as round 0017 stated nor
as dead as round 0018b implied. **It was measured with an instrument too coarse
to resolve it, at exactly one cell.** Everything below follows from asking why
that one cell, and what the answer implies.

## 2. Why `n23a1`: the variance is rho's, and rho's share of its own cost grows

Pollard rho's cost is a collision time — a random variable whose relative
spread is roughly scale-free. The index-calculus arm's cost is nearly
deterministic: it builds a fixed factor base, a fixed pair table, and collects
a fixed number of relations. So essentially all of the ratio's variance is
rho's, and how much of it reaches the ratio depends on how much of rho's total
cost is the collision search rather than the shared setup both arms pay.

| cell | r | rho solve gm (Ir) | sd(log solve) | solve share of rho | sd(log winner/rho) |
|:--|--:|--:|--:|--:|--:|
| `n13a0` | 2,003 | 594,915 | 0.060 | 0.710 | 0.037 |
| `n17a1` | 65,587 | 1,144,374 | 0.162 | 0.802 | 0.140 |
| `n19a0` | 130,873 | 1,484,762 | 0.216 | 0.792 | 0.171 |
| `n19a1` | 262,543 | 1,621,012 | 0.203 | 0.810 | 0.166 |
| `n23a0` | 2,095,853 | 2,512,920 | 0.203 | 0.852 | 0.214 |
| `n23a1` | 4,196,903 | 3,047,402 | 0.313 | 0.883 | 0.329 |
| `n29a1` | 42,457 | 1,391,238 | 0.122 | 0.736 | 0.088 |
| `n31a0` | 1,439,393 | 3,028,614 | 0.256 | 0.863 | 0.225 |

The spread of the per-cell ratio rises by a factor of nine across the panel,
ordered by subgroup order, and `n23a1` has the largest subgroup order on it.
A flat twelve fixtures a cell therefore buys a failure probability of about
`1e-8` at one end of the panel and 9% at the other. **The gate's resolution
degrades exactly as the cells get cryptographically interesting**, which is the
worst possible place for it to degrade.

## 3. The boundary: what the method's own cost model says it costs

`koblitz_tiny_ic.rs` publishes its cost model and minimises it at run time:

```
W(t) = t·size + (K + TABLE_PROBES) / (λ · cov(t))
λ = size(size+1) / 2r        cov(t) = 1 − ((K−t)/K)²
```

`size` is the number of signed base points and `K` the column count after the
signed-Frobenius orbits are folded, so `size = 2nK` for a base of `K` orbits on
a degree-`n` curve. The first term builds the pair table; the second scans for
relations, and it carries `r` in the numerator through `λ`.

**First, what `K` actually is.** The contract asks for `6·degree` points and
every round from 0002 on has been read as running with that base. It does not.
`build_subgroup_orbit_factor_base` samples orbits in batches of eight and stops
at the first rebuild that reaches the target, so any request at or below one
batch returns that whole batch. Measured (`round19_base_sweep.py`): at degree
23, `points=1`, `points=138` — the contract's value — and `points=368` all
return **eight orbits**. The panel's factor base is seven or eight orbits at
every cell, `size ≈ 16·degree` signed points, and **the contract's `points`
parameter has been inert since round 0002.** It is a batch boundary, not a
chosen size. The base can be enlarged in whole batches; it cannot be shrunk
below one.

This was found by measuring. The first version of this sweep varied
`6·degree` down to `2·degree` and up to `12·degree` and returned instruction
counts identical to six figures at every value — a flat line that the
pre-registration had predicted would fall by 2.6× at `n13a0`. The prediction
was wrong because the model was parameterised on a quantity the code ignores.
Both the model and this section are the corrected ones; the original is
recorded in this file's history and claimed nothing that survived.

**The boundary.** Treating `K` as free, the scan term is
`2r(K+3) / (4n²K²·cov(t))`. At fixed `t` and large `K`, `cov(t) → 2t/K` and the
scan term tends to `r/(4n²t)`: enlarging the base buys nothing, because
coverage thins exactly as fast as the pair density grows, while the table cost
`2ntK` grows linearly. Progress needs `t` to grow with `K`; at `t = αK` the
table costs `2nαK²` and the scan `~ cr/K`, so

```
K* ~ r^{1/3}/n,   size* ~ r^{1/3},   W* ~ r^{2/3}
```

against rho's `Θ((r/2n)^{1/2})`. **At its own optimal factor-base size this
method is `r^{2/3}` where rho is `r^{1/2}`: the ratio grows as `r^{1/6}` and
must eventually exceed one.** This is a boundary in the AGENTS.md §1 sense —
derived, not measured, and not movable by tuning constants. It is a statement
about this pair-table collector, not about index calculus in general.

Against the continuous optimum, on the base the panel really builds:

| cell | r | orbits built | size | r^{1/3} | size/optimal |
|:--|--:|--:|--:|--:|--:|
| `n13a0` | 2,003 | 7 | 182 | 12.6 | 14.44 |
| `n29a1` | 42,457 | 8 | 464 | 34.9 | 13.30 |
| `n17a1` | 65,587 | 8 | 272 | 40.3 | 6.74 |
| `n19a0` | 130,873 | 8 | 304 | 50.8 | 5.99 |
| `n19a1` | 262,543 | 8 | 304 | 64.0 | 4.75 |
| `n31a0` | 1,439,393 | 8 | 496 | 112.9 | 4.39 |
| `n23a0` | 2,095,853 | 8 | 368 | 128.0 | 2.88 |
| `n23a1` | 4,196,903 | 8 | 368 | 161.3 | 2.28 |

That last column falls monotonically, and it is the crossover mechanism drawn
out: `size` grows linearly in the degree while `r` grows exponentially in it,
so a base fourteen times over-provisioned at `n13a0` is barely twice the
optimum at `n23a1` and would be under-provisioned at the next degree up.
**`n23a1` is not an unlucky cell. It is the cell where the panel runs out of
the over-provisioning the other seven cells are winning on.**

**Measured, against the model.** Three fixtures a cell, every run
oracle-verified to return the incumbent's own logarithm, with rho as an
invariance control (rho never reads the factor base; its Ir moved by at most
3.05e-4, which is the cost of parsing a longer `points` field):

```
measured Ir over the orbit count, relative to the base the panel builds
cell     built      K=7      K=8     K=14     K=16     K=20     K=24     K=28     K=32     K=36     K=39     K=40
n13a0        7   1.0000        -   1.5303        -   2.0456        -   3.2116        -   4.0124        -        -
n17a1        8        -   1.0000        -   1.5190        -   2.0922        -   2.8670        -        -   3.5620
n19a0        8        -   1.0000        -   1.3995        -   1.9047        -   2.4492        -        -   3.3658
n19a1        8        -   1.0000        -   1.3955        -   1.7611(K=23)  2.2228        -   2.9743        -
n23a0        8        -   1.0000        -   1.6278        -   1.8477        -   3.0500        -        -   2.6596
n23a1        8        -   1.0000        -   1.2262        -   1.8076        -   2.2390        -        -   2.7965
n29a1        8        -   1.0000        -   1.6147        -   2.2510        -   2.9238        -        -   3.6533
n31a0        8        -   1.0000        -   1.4630        -   1.9865        -   2.7131        -        -   3.0829

measured / model, at the orbit counts both cover
cell       K=16     K=24     K=32     K=40
n17a1    0.8788   0.8427   0.8841   0.8892
n19a0    0.8497   0.8183   0.8120   0.9077
n19a1    0.9574        -        -        -
n23a0    1.2633   1.2468   1.7933   1.3784
n23a1    0.9984   1.1903   1.3527   1.5527
n29a1    0.8257   0.7722   0.7544   0.7553
n31a0    1.0531   1.0814   1.1804   1.1148
```

**The cost rises monotonically with the base at every one of the eight cells,
and the base the panel builds is the cheapest one reachable.** The model's
conclusion is confirmed by measurement. Its accuracy is not: it agrees within
25% at five cells and is out by up to 80% at `n23a0`, so it is reported as what
it is — a derivation that gets the direction and the rough magnitude right,
not a fit.

So there is no factor-base lever, at `n23a1` or anywhere on this panel: the
base is already at its floor, and its floor is its optimum. This is why round
0019 carries no new algorithmic arm. The measurement that would have justified
one says there is nothing there to take.

## 4. The phase ledger this rests on

Geometric means over the 24 frozen fixtures a cell, instructions:

| cell | r | fb+tables | collection | individual log | shared fixed | rho solve |
|:--|--:|--:|--:|--:|--:|--:|
| `n13a0` | 2,003 | 244,197 | 39,500 | 21,946 | 222,543 | 594,915 |
| `n17a1` | 65,587 | 333,086 | 85,753 | 36,827 | 265,637 | 1,144,374 |
| `n19a0` | 130,873 | 425,710 | 133,362 | 43,720 | 378,549 | 1,484,762 |
| `n19a1` | 262,543 | 431,441 | 228,912 | 57,498 | 368,160 | 1,621,012 |
| `n23a0` | 2,095,853 | 755,875 | 780,210 | 117,159 | 426,310 | 2,512,920 |
| `n23a1` | 4,196,903 | 983,731 | 1,366,827 | 134,674 | 402,291 | 3,047,402 |
| `n29a1` | 42,457 | 879,527 | 75,300 | 44,491 | 482,629 | 1,391,238 |
| `n31a0` | 1,439,393 | 901,375 | 556,216 | 87,611 | 477,464 | 3,028,614 |

The two same-degree pairs are the clean `r`-only contrast, because `size` and
the field arithmetic are identical within a pair:

| pair | r ratio | collection ratio | rho solve ratio | winner/rho ratio |
|:--|--:|--:|--:|--:|
| `n19a0` → `n19a1` | 2.006 | 1.716 | 1.092 | 1.045 |
| `n23a0` → `n23a1` | 2.002 | 1.752 | 1.213 | 1.166 |

Doubling the subgroup order at a fixed field costs the candidate 72–75% more
collection and rho only 9–21% more, so the ratio rises 4.5% per doubling at
degree 19 and 16.6% at degree 23. From 0.8913 at `n23a1`, **one further
doubling of `r` at degree 23, at the degree-23 rate, reads 1.040: above rho.**

Two honesty notes on that extrapolation. It rests on two points per degree, so
it is a rate between two cells and not a fit; AGENTS.md §6 marks it as an
extrapolation and this is that mark. And it is *slower* than the model's
`r^{1/2}` prediction for the same contrast — collection grows as `r^{0.80}` at
degree 19 and `r^{0.81}` at degree 23, not `r^{1.0}`, so the model of §3 is an
upper bound on how fast the candidate loses ground at these sizes rather than a
description of it. The direction is what both agree on.

## 5. The protocol amendment: a per-cell confirmation allocation

`tournament.py prepare` gains `--confirmation-cases cell=count,…`. It is
**additive and raises only**: a count below the profile floor is refused, an
unnamed cell keeps the flat profile, and every round before 0019 reproduces
byte-for-byte with the flag absent. The estimator and both gates are untouched
— the cross-cell mean is unweighted, so a cell measured 40 times carries
exactly the weight of a cell measured 12 times, and the per-cell gate stays a
maximum rather than an average. Two tests in `test_tournament.py` hold those
two properties.

The allocation is computed from frozen rounds 0017 and 0018b alone, by
`campaign_20260916/round19_allocate.py`, as the smallest multiple of four whose
one-sided failure probability is at most **1%** in *both* metrics:

| cell | r | Ir `both`/rho | per-case sd | native `both`/rho | per-case sd | P(fail) at 12 | cases | P(fail) at that count |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| `n13a0` | 2,003 | 0.6916 | 0.037 | 0.8612 | 0.053 | 0.0000 | **12** | 0.0000 |
| `n17a1` | 65,587 | 0.5689 | 0.140 | 0.8625 | 0.058 | 0.0000 | **12** | 0.0000 |
| `n19a0` | 130,873 | 0.5737 | 0.171 | 0.8426 | 0.057 | 0.0000 | **12** | 0.0000 |
| `n19a1` | 262,543 | 0.5911 | 0.166 | 0.8631 | 0.062 | 0.0000 | **12** | 0.0000 |
| `n23a0` | 2,095,853 | 0.7339 | 0.214 | 0.8984 | 0.080 | 0.0000 | **12** | 0.0000 |
| `n23a1` | 4,196,903 | 0.8809 | 0.329 | 0.9577 | 0.108 | 0.0910 | **40** | 0.0074 |
| `n29a1` | 42,457 | 0.7780 | 0.088 | 0.8862 | 0.061 | 0.0000 | **12** | 0.0000 |
| `n31a0` | 1,439,393 | 0.5804 | 0.225 | 0.8308 | 0.098 | 0.0000 | **12** | 0.0000 |

**124 confirmation cases against the flat panel's 96. No cell is measured
less than round 0018b measured it.** The whole increase goes to `n23a1`.

Why this cannot buy a pass, stated plainly: more fixtures move a cell's
estimate toward its true value in whichever direction that lies. If the
candidate is really above rho at `n23a1`, forty fixtures make that verdict
*more* certain than twelve, not less. The allocation removes a coin flip; it
does not choose which way the coin was weighted.

One thing the allocation cannot fix, and which round 0019 states rather than
hides: **confirmation and replay share fixtures, and instruction counts are
deterministic, so the two stages are one instruction test performed twice, not
two independent ones.** Round 0018b's per-cell Ir agreed to four decimals
across the two stages for exactly this reason. Replay tests native
reproducibility; it adds nothing to the instruction evidence.

## 6. The round

| | |
|:--|:--|
| seed | `2026092119` (fresh; rounds 0017 and 0018b used 2026091717 and 2026091818) |
| panel | `--cells 13a0,17a1,19a0,23a0,23a1,31a0 --holdout-cells 19a1,29a1` |
| allocation | `--confirmation-cases n23a1=40` |
| profile | pilot, one target, objective `rho`, `--require-native-progress` |
| arms | `incumbent`, `both` |
| trials | 2,646 — 36 aa, 54 smoke, 162 development, 162 selection, 1,116 confirmation, 1,116 replay |

`both` is round 0018's two levers together, rebuilt from the committed patches
on the frozen baseline: `round18-block.patch` (the scan block's clamp ceiling,
64 → 16) and `round18-column.patch` (a column is the orbit representative
itself rather than `[h]` of it, with the matching factor of `h` gone from the
row and the descent, declaring `column_convention: representative`). Round
0018b measured it at 0.9652 [0.9448, 0.9828] against the incumbent in
instructions and 0.9613 [0.9288, 0.9895] natively, passing the no-regression
gate on both stages and failing only the strict rho gate at `n23a1`.

No new arm. Round 0019 is a resolution round, and §3 is the measurement that
says a new arm would have nothing to take.

## 7. Predictions

Registered before the round. Each is falsifiable from round 0019's own frozen
outputs.

1. **`both` beats the incumbent again.** Instruction ratio in [0.95, 0.98]
   with the 95% upper limit below 1, native upper limit below 1, on both final
   stages. Round 0018b: 0.9652 and 0.9613.
2. **The `n23a1` instruction ratio lands inside the pooled band.**
   `both`/rho at `n23a1` in [0.78, 1.00], i.e. below both round 0018b's 1.0112
   and above round 0017's implied value. Point prediction 0.881.
3. **`beats_rho_strict` holds.** Every cell below 1 in both metrics on both
   final stages. Pre-registered probability ≥ 0.97 if §5's margins and spreads
   are right; a failure at any cell other than `n23a1` falsifies the spread
   model outright.
4. **The per-cell spread reproduces.** Each cell's per-case sd(log winner/rho)
   in instructions within ±40% of the §2 column. This is the check that §5's
   allocation was drawn from a stable quantity and not from two lucky rounds.
5. **The ordering by subgroup order survives.** sd(log winner/rho) is
   monotonically increasing in `r` across `n13a0`, `n17a1`, `n19a0`, `n19a1`,
   `n31a0`, `n23a1` — allowing `n29a1` and `n23a0` out, since 0017/0018b
   already put `n23a0` marginally above `n19a1` and `n31a0`.
6. **The same-degree contrast holds.** winner/rho at `n23a1` divided by
   winner/rho at `n23a0` is above 1.05 (frozen value 1.166), and winner/rho at
   `n19a1` over `n19a0` is above 1.00 (frozen value 1.045 — the weaker bound is
   because the degree-19 rate is already close to 1 and forty fixtures do not
   go to that pair). This is §4's crossover mechanism, re-measured on fresh
   fixtures.
7. **Confirmation and replay agree exactly in instructions.** Every cell's Ir
   ratio equal to four decimal places across the two stages. If they differ,
   the instruction metric is not deterministic and §5's whole variance account
   is wrong.

## 8. Falsification, and what the round is not allowed to claim

Any bad, missing or unmatched certificate; any arm returning a different
logarithm or factor-base hash from the incumbent on the same fixture; a `both`
report failing to declare `column_convention: representative`; the A/A control
outside [0.95, 1.05] at any cell; a prepare whose contract differs from §6 in
seed, panel, allocation, profile, objective or limits.

**What a pass would and would not mean.** A `beats_rho_strict` on this panel
would say the candidate is cheaper than matched rho, in instructions and in
cold native wall, on eight Koblitz cells with subgroup orders from 2·10³ to
4·10⁶, under one solver, one compiler and one profiler. It would **not** be an
exponent result, and §3 is the reason: the method's own cost model puts it at
`r^{2/3}` against rho's `r^{1/2}`, its factor base is already at its optimum at
the largest cell on the panel, and the same-degree contrast in §4 measures the
ratio rising 9–18% per doubling of `r`. By the AGENTS.md §3 test the campaign's
gains from round 0006 to here are **engineering** — `S` fell, the ratio to the
boundary at the largest cell did not — and a win here must be labelled that way
on the scoreboard whichever way the gate goes.

Round 0019 also may not claim that round 0017's or round 0018b's records were
wrong. Both stand. They were two draws from one distribution this round is
built to resolve, and the correction is to the campaign's reading of them, not
to the records.
