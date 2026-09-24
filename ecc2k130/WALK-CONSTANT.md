# Iterations per solve: the σ walk against the table walk

Question: the campaign has two iteration functions: the shipping walk
`R ← R + σʲ(R)` and the table walk `R ← R + ε·σᵏ(T_h)` of
[ITERATION-FUNCTION.md](ITERATION-FUNCTION.md).  The table walk runs
9–15% faster per update in every paired measurement, but the choice between
them was made on rate alone, assuming both need the same number of
iterations up to a small branch constant.  How many iterations does each
really need, and which is cheaper per solve?  The choice has to be made
before the first distinguished point is collected, because the two walks
cannot share a corpus.

Answer, after round 2 (§11): **`maxIters` is raised, and the table walk's
cycle rule now refuses the fruitless cycles that disqualified it.  Priced end
to end it is projected at 0.81–0.86× the σ walk's cost per solve.  The
campaign still runs σ until one measurement is in: the new kernel's paired
rate on the card.**

- **Round 2.** The rule now refuses a step that undoes any of the last four,
  or that closes a τ-relation with the last three.  No fruitless cycle with
  four or fewer determined tags survives it, and every one the old rule let
  through is gone.  What it leaves is counted at `2 × 10⁻⁵` of trails trapped
  at `dpWeight = 32`: mostly six-step relations such as
  `τ⁴ + τ³ + τ² − τ − 2 = 0`, one of which was seen.  The table walk's `c` is
  unchanged.  Walks that meet with different histories part 0.97% of the time,
  which costs ×1.0048 in iterations.  The projection uses the old kernel's
  paired rates; the switch pays unless the new rule costs the kernel 14–19% of
  its rate.  `aws/campaign.json` still names σ, and the switch is the campaign
  owner's.
- **`maxIters`** is `2^32` in `aws/campaign.json` and in the witness
  generator's default, up from `2^30`.  The bucket copy, which is what
  workers read, is the owner's to republish.

Round 1's answer was **keep the σ walk, and raise `maxIters`**, on these
findings (the fruitless-cycle ones are about the rule as it was then):

- **Iterations.** The table walk needs `c = 1.00` times a random mapping's
  iterations at every degree measured.  With ECC2K-130's branch distribution
  the σ walk needs `c = 1.08–1.10` at `n = 23–59`, falling with the degree
  towards its first-order value `1.070`.  At `n = 131` it lies between
  `1.070` and `1.082` (an extrapolation), so Bailey et al.'s `2^60.9` holds
  for it.  In iterations alone, the table walk needs 6–7% fewer.
- **Fruitless cycles.** The table walk as built lets fruitless cycles through
  its cycle rule.  Frobenius on these curves satisfies `σ² + σ + 2 = 0`, so
  four steps in one branch, adding `σᵏ⁺²T`, `σᵏ⁺¹T`, `σᵏT` and `σᵏT`, return
  a walk to where it started.  No two of those tags cancel, so the rule, which
  looks only for pairs, never fires.  Neither does it catch six-step pairwise
  cycles.  The rate of both is measured on the device's own walk and on an
  emulation, and tracks a count of the patterns: 1145 τ-relation returns
  seen against 1312 predicted.  The count runs about 14% high.
- **What the cycles cost.** At the live bucket's `dpWeight = 32` those
  cycles trap **about half** of the table walk's trails (50–54%, H = 8).
  Each trapped trail runs to the `maxIters = 2^30` guard and is discarded,
  so a completed trail costs **7.1–8.6×** its steps.  Priced end to end, the
  table walk as built costs **4.85–6.26× the σ walk per solve**.  With every
  fruitless cycle caught and escaped (a check that does not exist yet), it
  would cost **0.81–0.86×**; that figure is a projection.
- **The σ walk's own loss.** The same guard costs the σ walk **15.6% of its
  steps** today.  At `dpWeight = 32` a trail averages `2^28.41` steps, so
  `2^30` is only three trail lengths.  About 5% of honest trails reach the
  guard first, and they are the longest.  `maxIters = 2^32` brings that loss
  to 0.01%.  No distinguished point has been collected
  (`docs/ecc2k130-status/history.json` is empty), so the change costs
  nothing now.

Everything below follows `AGENTS.md`: §1 gives the boundary and unit, §2 the
method, §3 why the walks differ, §4 the single table, §5 the fruitless
cycles, §6 the cost to solve, §7 the falsification target, §8 the
classification, and §9 what this changes elsewhere.  §11 is round 2: the
cycle rule extended, with its own target, results, cost and classification.

## 1. Boundary and unit

**Unit.** `c` is iterations over a random mapping's on the same class set.
Rho on the classes of `⟨σ, −1⟩` (`N = (ℓ − 1)/2n` of them) visits
`1 + Q(N) ≈ √(πN/2)` distinct classes before its first collision if the
walk is a random mapping; `c` is the measured mean over that.  At ECC2K-130,
`√(πN/2) = 2^60.809`.

**Floor.** `c = 1`, derived.  A walk that keeps the class structure and
whose collisions are not trivial cannot beat a random mapping on average
without making its images more concentrated than a random mapping's, and
neither walk tries to.  Two first-order models sit above the floor:

- a walk that is **injective on each branch** pays `1/√(1 − Σp²)` (Teske;
  Brent–Pollard), where `p` is the branch distribution.  The σ walk is one:
  `R ↦ (1 + λʲ)R` is a bijection;
- the table walk adds `T_h` in a frame fixed by the class, so two classes on
  the same branch collide unless their frames agree.  Only `1/2n` of
  same-branch collisions are suppressed, and the penalty is
  `1/√(1 − Σp²/2n) ≈ 1.0002` at `n = 131`.

**Reference.** The σ walk, measured the same way on the same curves.
[ITERATION-FUNCTION.md](ITERATION-FUNCTION.md) §4.4 modelled the table walk
by the r-adding constant `1 + 1/(2H)`, which is the injective model; §3 shows
it does not apply to this walk.

**Cost to solve**, the quantity §6 decides on: iterations times one over
rate, with every loss priced.  That is `c × (1/rate) × (steps per step on a
completed trail)`, and the last factor is where the fruitless cycles and the
guard enter.

## 2. Method

A **trial** starts `W` walks from random starts `[r]P`, steps them in turn,
and stops when one lands on a class any of them has visited (an exact set
of canonical class keys, so there are no distinguished points and no
tail).  The count is the number of distinct classes visited.  Many walks
rather than one, because a single walk's rho length is set by the cycle of
the component it lands in, which one mapping does not average away.  Each
trial draws a fresh mapping where the walk has one to draw (the table
`T_h`, the branch-hash salt), so trials average over mappings as well as
starts; the σ walk with the device's own branch choice has exactly one
mapping per curve.

Two harnesses:

- **Emulation**, `examples/ecc2k130_walk_constant.rs`, degrees 19–59.  It
  builds the table walk as the device does (the point itself, the step tag
  `(h, k, ε)`, and the device's tag-history rule), except for the frame: `k`
  and `ε` come from the canonical representative of the class, where the
  device uses a phase and a pivot bit.  Both frames satisfy
  `k(σR) = k(R) + 1` and `ε(−R) = 1 − ε(R)`, which is all the class structure
  and the rule use.
- **Device**, `src/walkconstant.cpp` (`make walk-constant`), at `n = 23` and
  `41`.  These are ecc2k130's own reference walks exactly as
  `Solver::rewalk` calls them: `Ref::step` for σ and `TableWalk::step` for
  the table, with branches from the device's type-II normal-basis weight.

**Branch distributions.** Both harnesses take three:

- `native`: `(HW(x)/2) mod H`, the device's choice;
- `uniform`: a salted hash of the class;
- `ecc2k130`: a salted hash mapped onto the exact probabilities
  `(HW/2) mod H` has at `n = 131`.

`native` is not ECC2K-130's distribution at the degrees that run: `HW/2`
has a standard deviation of `√n/4`, so `Σp²` is 0.23 at `n = 23` against
0.127 at `n = 131`.  The `ecc2k130` rows are the ones that carry over.

**Fruitless returns.** Every table step compares the new point with the
walk's last 16 points.  A return whose tags sum to zero formally is counted
by kind (§5) and length, and the walk restarts; it is not scored as a
collision.  Returns that are not formal are genuine collisions and are
scored as such.

**Controls.**

- An adding walk with 256 uniform branches must give `c = 1`.
- `W = 4, 16, 64` at `n = 37` must agree.
- The device and the emulation must agree at `n = 23`.
- The coset index of `⟨1 + λʲ, λ, −1⟩` in `Z_ℓ^*`, which would confine the σ
  walk to a coset and cost it a factor no branch statistic shows, is **1**
  at every degree run and at ECC2K-130
  ([coset_index.txt](benchmarks/walk-constant/coset_index.txt)).

## 3. Why the walks differ

The σ walk multiplies the scalar by `1 + λʲ`.  On one branch that is a
bijection of the class set, so two classes on the same branch never
collide, and the walk pays the full injective penalty
`1/√(1 − Σp²) = 1.070`.  Because the multipliers commute it also pays a
higher-order penalty: `(1 + λᵃ)(1 + λᵇ)` is the same step in either order,
which forbids the two-step collisions a random mapping would have.  The
emulation measures that excess at `c/model = 1.011–1.028` with hashed
branches, shrinking with the degree (§4).

The table walk adds `ε σᵏ(T_h)` with `k, ε` read from the class.  Two
classes on one branch land in the same class whenever their frames differ,
so on the class set it is not injective per branch.  It behaves like a
random mapping up to `Σp²/2n`.

## 4. The single table

All rows are in unit `c`; the model column is §1's model for that walk.
Sources are the frozen files in
[benchmarks/walk-constant](benchmarks/walk-constant/), printed by
`summarize.py`.  The σ rows at `n = 23, 37, 41` come from `matrix-v1.jsonl`,
whose σ code is unchanged; everything else comes from `matrix-v2.jsonl` and
the device files.  v1's table and adding rows used only a previous-class
cycle rule and scored fruitless 4-cycles as collisions.  At `n = 41`, with
102k steps per trial, that ended about 14% of trials early and gave `c =
0.909 ± 0.006`, which is the prediction for that artefact (`1 − 0.64 λμ`
with `λ = q²`).  Those rows are superseded, not deleted.

| n | walk | H | branches | harness | W | trials | c | ± | model | c / model | source |
|---:|---|---:|---|---|---:|---:|---:|---:|---:|---:|---|
| 19 | table | 8 | ecc2k130 | emulation | 8 | 199,911 | 1.0028 | 0.0012 | 1.0017 | 1.0011 | matrix-v2.jsonl |
| 19 | table | 8 | native | emulation | 8 | 199,891 | 1.0064 | 0.0012 | 1.0034 | 1.0030 | matrix-v2.jsonl |
| 19 | sigma | 8 | ecc2k130 | emulation | 8 | 200,000 | 1.0063 | 0.0012 | 1.0700 | 0.9405 | matrix-v2.jsonl |
| 19 | sigma | 8 | native | emulation | 8 | 200,000 | 1.2130 | 0.0014 | 1.1602 | 1.0455 | matrix-v2.jsonl |
| 23 | adding | 256 | uniform | emulation | 8 | 199,974 | 1.0002 | 0.0012 | 1.0000 | 1.0001 | matrix-v2.jsonl |
| 23 | table | 8 | ecc2k130 | device | 8 | 19,995 | 1.0076 | 0.0037 | 1.0014 | 1.0062 | device-v2.jsonl |
| 23 | table | 8 | native | device | 8 | 59,989 | 1.0080 | 0.0021 | 1.0026 | 1.0054 | device-v2.jsonl |
| 23 | table | 16 | native | device | 8 | 19,991 | 1.0041 | 0.0037 | 1.0026 | 1.0016 | device-v1.jsonl |
| 23 | table | 8 | ecc2k130 | emulation | 8 | 199,974 | 1.0037 | 0.0012 | 1.0014 | 1.0023 | matrix-v2.jsonl |
| 23 | table | 8 | native | emulation | 8 | 199,975 | 1.0050 | 0.0012 | 1.0026 | 1.0024 | matrix-v2.jsonl |
| 23 | table | 8 | uniform | emulation | 8 | 199,978 | 1.0026 | 0.0012 | 1.0014 | 1.0012 | matrix-v2.jsonl |
| 23 | table | 16 | ecc2k130 | emulation | 8 | 199,974 | 1.0021 | 0.0012 | 1.0011 | 1.0010 | matrix-v2.jsonl |
| 23 | table | 16 | native | emulation | 8 | 199,978 | 1.0046 | 0.0012 | 1.0026 | 1.0020 | matrix-v2.jsonl |
| 23 | sigma | 8 | ecc2k130 | device | 8 | 20,000 | 1.0948 | 0.0040 | 1.0700 | 1.0232 | device-v2.jsonl |
| 23 | sigma | 8 | native | device | 8 | 20,000 | 1.1595 | 0.0041 | 1.1415 | 1.0157 | device-v1.jsonl |
| 23 | sigma | 8 | native | device | 8 | 20,000 | 1.1705 | 0.0041 | 1.1415 | 1.0255 | device-v2.jsonl |
| 23 | sigma | 8 | uniform | device | 8 | 20,000 | 1.0951 | 0.0040 | 1.0690 | 1.0244 | device-v2.jsonl |
| 23 | sigma | 8 | ecc2k130 | emulation | 8 | 200,000 | 1.0995 | 0.0013 | 1.0700 | 1.0276 | matrix-v1.jsonl |
| 23 | sigma | 8 | native | emulation | 8 | 200,000 | 1.1537 | 0.0012 | 1.1422 | 1.0100 | matrix-v1.jsonl |
| 23 | sigma | 8 | uniform | emulation | 8 | 200,000 | 1.0971 | 0.0013 | 1.0690 | 1.0262 | matrix-v1.jsonl |
| 31 | table | 8 | ecc2k130 | emulation | 8 | 199,971 | 1.0035 | 0.0012 | 1.0010 | 1.0025 | matrix-v2.jsonl |
| 31 | table | 8 | native | emulation | 8 | 199,975 | 1.0029 | 0.0012 | 1.0016 | 1.0013 | matrix-v2.jsonl |
| 31 | sigma | 8 | ecc2k130 | emulation | 8 | 200,000 | 1.0964 | 0.0013 | 1.0700 | 1.0247 | matrix-v2.jsonl |
| 31 | sigma | 8 | native | emulation | 8 | 200,000 | 1.1560 | 0.0013 | 1.1204 | 1.0318 | matrix-v2.jsonl |
| 37 | adding | 256 | uniform | emulation | 16 | 40,000 | 1.0030 | 0.0026 | 1.0000 | 1.0030 | matrix-v2.jsonl |
| 37 | table | 8 | ecc2k130 | emulation | 4 | 39,999 | 1.0052 | 0.0026 | 1.0009 | 1.0043 | matrix-v2.jsonl |
| 37 | table | 8 | ecc2k130 | emulation | 16 | 40,000 | 1.0036 | 0.0026 | 1.0009 | 1.0027 | matrix-v2.jsonl |
| 37 | table | 8 | ecc2k130 | emulation | 64 | 39,999 | 1.0041 | 0.0026 | 1.0009 | 1.0033 | matrix-v2.jsonl |
| 37 | table | 8 | native | emulation | 16 | 39,998 | 1.0015 | 0.0026 | 1.0013 | 1.0002 | matrix-v2.jsonl |
| 37 | sigma | 8 | ecc2k130 | emulation | 4 | 40,000 | 1.0939 | 0.0028 | 1.0700 | 1.0223 | matrix-v1.jsonl |
| 37 | sigma | 8 | ecc2k130 | emulation | 16 | 40,000 | 1.0949 | 0.0028 | 1.0700 | 1.0233 | matrix-v1.jsonl |
| 37 | sigma | 8 | ecc2k130 | emulation | 64 | 40,000 | 1.0870 | 0.0028 | 1.0700 | 1.0159 | matrix-v1.jsonl |
| 37 | sigma | 8 | native | emulation | 16 | 40,000 | 1.2170 | 0.0033 | 1.1075 | 1.0988 | matrix-v1.jsonl |
| 41 | adding | 256 | uniform | emulation | 16 | 8,000 | 1.0089 | 0.0059 | 1.0000 | 1.0089 | matrix-v2.jsonl |
| 41 | table | 8 | ecc2k130 | device | 16 | 300 | 1.0051 | 0.0299 | 1.0008 | 1.0043 | device-v2.jsonl |
| 41 | table | 8 | ecc2k130 | emulation | 16 | 8,000 | 1.0125 | 0.0059 | 1.0008 | 1.0118 | matrix-v2.jsonl |
| 41 | table | 8 | native | emulation | 16 | 8,000 | 1.0084 | 0.0059 | 1.0011 | 1.0073 | matrix-v2.jsonl |
| 41 | sigma | 8 | ecc2k130 | device | 16 | 300 | 1.0997 | 0.0336 | 1.0700 | 1.0278 | device-v2.jsonl |
| 41 | sigma | 8 | ecc2k130 | emulation | 16 | 8,000 | 1.0910 | 0.0064 | 1.0700 | 1.0196 | matrix-v1.jsonl |
| 59 | table | 8 | ecc2k130 | emulation | 16 | 20,000 | 0.9949 | 0.0037 | 1.0005 | 0.9943 | matrix-v2.jsonl |
| 59 | table | 8 | native | emulation | 16 | 20,000 | 1.0050 | 0.0037 | 1.0006 | 1.0043 | matrix-v2.jsonl |
| 59 | sigma | 8 | ecc2k130 | emulation | 16 | 20,000 | 1.0817 | 0.0040 | 1.0700 | 1.0109 | matrix-v2.jsonl |
| 59 | sigma | 8 | native | emulation | 16 | 20,000 | 1.1221 | 0.0042 | 1.0850 | 1.0342 | matrix-v2.jsonl |

**Reading it.**

- The table walk is at the floor: `c = 0.995–1.013` over every emulation
  row (pooled `1.0035 ± 0.0006` with ECC2K-130's branches at H = 8).  It shows
  no dependence on `H`, on the branch distribution (`Σp²` from 0.10 to
  0.26), on `W`, or on the degree, and the device and the emulation agree.
- The σ walk with ECC2K-130's distribution is at `c = 1.0995, 1.0964,
  1.0919, 1.0910, 1.0817` at `n = 23, 31, 37, 41, 59`.  That is the
  first-order `1.070` times a higher-order excess (the commuting multipliers
  of §3) which shrinks with the degree, from 2.8% to 1.1%.  The first-order
  value is its limit, so at `n = 131` the constant lies in `[1.070, 1.082]`.
  That is an extrapolation, and §6 carries the whole bracket.
- `n = 19` is excluded for σ.  Two of its eight multipliers fall in one
  class (`1 + λ¹⁰ = λ⁻⁹(1 + λ⁹)` when `n = 19`; every `n ≥ 21` has eight
  distinct ones), and with 3,444 classes its short walks return to their
  own classes early.  A pure-arithmetic model of the σ walk on classes
  ([sigma_classes.py](benchmarks/walk-constant/sigma_classes.py),
  [output](benchmarks/walk-constant/sigma_classes.txt)) reproduces the
  `1.006`, and gives `1.087 ± 0.011` at `n = 23`: a third implementation,
  agreeing with the other two.
- One mapping at a time, both walks' constants stray from the average over
  mappings far more than their standard errors allow: σ 1.108 with a
  standard deviation of 0.034 over 6 mappings, each ± 0.003; table 0.995
  with a standard deviation of 0.041 over 3 mappings, each ± 0.003 at `n =
  37`, and σ 1.075 with a standard deviation of 0.043 over 4 mappings, each
  ± 0.005; table 1.021 with a standard deviation of 0.040 over 4 mappings,
  each ± 0.005 at `n = 59`.
- Most of that spread belongs to the statistic, not to either walk.
  Uniformly random functions on 45,562 points, under the same count, stray
  by a 6% standard deviation at `W = 8`, 2.3% at `W = 16` and 1.1% at
  `W = 64` ([mapping_spread.py](benchmarks/walk-constant/mapping_spread.py),
  [output](benchmarks/walk-constant/mapping_spread.txt)).  `W = 16` there
  has the same ratio of walk length to `√N` as the walks' `n = 37` rows.
  The walks' 3.4–4.3% is of the same order and somewhat above; each spread
  rests on only 4–6 mappings.  A few walks, each a sizeable fraction of `√N`
  long, sample only a few neighbourhoods of one mapping's trees, and more,
  shorter walks average them.
- A campaign is at the far end of that: about `2^33` trails of `2^28`
  steps against `√N ≈ 2^61`.  So the mapping average, the main table's
  rows, is the constant a campaign pays, and every figure below uses it.
- With the device's own weight branches, σ is at `1.122–1.217` at
  `n = 23–59`, with `c/model` 1.010–1.099.  The device's branch function
  is fixed, so each of those rows is a single mapping.  Their scatter,
  including `n = 37`'s 1.217, is within the one-mapping spread at their
  `W`, so they establish no effect of weight-based branches beyond the
  model column's `Σp²`.
- §6 uses the hashed σ bracket.  That is the conservative choice for
  switching, since any extra cost of the native branches would only make σ
  dearer.

| n | walk | mappings (seed: c) | mean | spread (sd) | per-mapping ± |
|---:|---|---|---:|---:|---:|
| 37 | sigma ecc2k130 H = 8 | 1: 1.0810, 2: 1.0752, 3: 1.1057, 4: 1.0960, 5: 1.1680, 6: 1.1239 | 1.1083 | 0.0341 | 0.0029 |
| 37 | table ecc2k130 H = 8 | 1: 0.9810, 2: 1.0412, 3: 0.9630 | 0.9951 | 0.0409 | 0.0026 |
| 59 | sigma ecc2k130 H = 8 | 1: 1.0226, 2: 1.0573, 3: 1.1161, 4: 1.1046 | 1.0752 | 0.0433 | 0.0055 |
| 59 | table ecc2k130 H = 8 | 1: 1.0667, 2: 0.9886, 3: 1.0427, 4: 0.9850 | 1.0208 | 0.0405 | 0.0054 |

## 5. Fruitless cycles

The table walk is additive, so a run of steps whose addends sum to zero
returns a walk to a point it has left, with no collision.  The device's rule
(`eccTagFruitless` in `include/tablewalk.h`) refuses two such runs: a step
that undoes the last one, and a pairwise 4-cycle.  Enumerating the patterns
it passes
([fruitless_patterns.py](benchmarks/walk-constant/fruitless_patterns.py),
[output](benchmarks/walk-constant/fruitless_patterns.txt)) finds, at the
leading order `(2n)⁻³`:

- **24 four-step τ-relation cycles.** Frobenius satisfies `τ² + τ + 2 = 0`
  on these curves, hence also `τ³ + τ − 2 = 0`.  So four steps in one
  branch whose addends are `σᵏ⁺²T, σᵏ⁺¹T, σᵏT, σᵏT` (all one sign), or
  `σᵏ⁺³T, σᵏ⁺¹T, −σᵏT, −σᵏT`, in any of 12 orders, sum to `O`.  No two of
  the tags cancel, so the rule never fires.  All four steps share one
  branch, so a walk enters one with probability `24 Σp⁴/(2n)³` per step.
- **4 six-step pairwise cycles**, such as `a, b, c, −a, −b, −c`, entered
  with probability `4 (Σp²/2n)³` per step.
- nothing at 2, 3 or 5 steps.

A trapped walk repeats its cycle forever: the tag history at the start of
each lap is the same.  It never becomes distinguished.  The σ walk has no
fruitless cycles, since `∏(1 + λʲ) = ±λᵃ` has no formal solutions.

Measured against that count:

| n | walk | H | branches | harness | steps | pairwise seen | expected | τ-relation seen | expected |
|---:|---|---:|---|---|---:|---:|---:|---:|---:|
| 19 | table | 8 | ecc2k130 | emulation | 13,285,870 | 0 | 2.0 | 7 | 12.2 |
| 19 | table | 8 | native | emulation | 13,338,624 | 9 | 16.5 | 74 | 140.4 |
| 23 | adding | 256 | uniform | emulation | 52,046,602 | 0 | 0.0 | 0 | 0.0 |
| 23 | table | 8 | ecc2k130 | device | 5,243,610 | 0 | 0.4 | 3 | 2.7 |
| 23 | table | 8 | native | device | 15,737,837 | 5 | 8.3 | 65 | 70.0 |
| 23 | table | 8 | ecc2k130 | emulation | 52,232,792 | 0 | 4.4 | 36 | 27.0 |
| 23 | table | 8 | native | emulation | 52,301,181 | 18 | 27.6 | 227 | 231.6 |
| 23 | table | 8 | uniform | emulation | 52,173,746 | 5 | 4.2 | 18 | 25.1 |
| 23 | table | 16 | ecc2k130 | emulation | 52,149,589 | 0 | 2.1 | 16 | 17.3 |
| 23 | table | 16 | native | emulation | 52,281,030 | 19 | 27.6 | 224 | 231.5 |
| 31 | table | 8 | ecc2k130 | emulation | 36,858,714 | 0 | 1.3 | 8 | 7.8 |
| 31 | table | 8 | native | emulation | 36,836,568 | 6 | 5.1 | 42 | 43.2 |
| 37 | adding | 256 | uniform | emulation | 88,149,303 | 0 | 0.0 | 0 | 0.0 |
| 37 | table | 8 | ecc2k130 | emulation | 88,821,893 | 1 | 1.8 | 6 | 11.0 |
| 37 | table | 8 | ecc2k130 | emulation | 88,204,046 | 2 | 1.8 | 8 | 11.0 |
| 37 | table | 8 | ecc2k130 | emulation | 86,328,362 | 2 | 1.7 | 10 | 10.7 |
| 37 | table | 8 | native | emulation | 88,011,299 | 12 | 5.5 | 46 | 46.4 |
| 41 | adding | 256 | uniform | emulation | 828,193,268 | 0 | 0.0 | 0 | 0.0 |
| 41 | table | 8 | ecc2k130 | device | 30,937,515 | 0 | 0.5 | 4 | 2.8 |
| 41 | table | 8 | ecc2k130 | emulation | 831,148,159 | 16 | 12.2 | 77 | 75.9 |
| 41 | table | 8 | native | emulation | 827,739,841 | 28 | 32.9 | 257 | 275.2 |
| 59 | table | 8 | ecc2k130 | emulation | 229,988,097 | 2 | 1.1 | 6 | 7.1 |
| 59 | table | 8 | native | emulation | 232,323,976 | 1 | 1.9 | 11 | 15.2 |

**Reading it.**

- The τ-relation rate follows `24Σp⁴/(2n)³` on the device walk and in the
  emulation, across a factor of 344 in predicted rate: 1145 seen against
  1312 predicted, 0.87 ± 0.03 overall, and 425 seen against 503 predicted,
  0.85 ± 0.04 on the rows at `n ≥ 37`, where walks run thousands of steps.
  The device walk confirms it directly: 65 seen against 70 predicted at `n =
  23`.
- The pairwise term follows `4(Σp²/2n)³`: 126 seen against 166 predicted,
  0.76 ± 0.07 overall, 64 seen against 67 predicted, 0.95 ± 0.12 at `n ≥
  37`.  At `n = 19–23` short walks cut six-step windows off.
- The measured total runs at 0.86 of the count.  The count is leading
  order, and the rule's advancing of `h` breaks a few entries it does not
  model.  §6 prices both the count and the count scaled to the
  measurements; no conclusion depends on which.

## 6. Cost to solve

**Trap loss.** A walk becomes distinguished with probability `θ` per step
and enters a fruitless cycle at rate `r`.  A trapped walk runs to `maxIters`,
emits nothing (`PACKED.md`: overdue walks restart and emit no distinguished
point), and is restarted.  Only the points of trails that end distinguished
are in the corpus, so the loss is total steps over steps on completed
trails ([trap_cost.py](benchmarks/walk-constant/trap_cost.py),
[output](benchmarks/walk-constant/trap_cost.txt)).  This is a model, built
on the rates §5 measures; its outputs are marked as model rows.  At
`n = 131`, H = 8, ECC2K-130's distribution:

| `dpWeight` (trail) | walk | trails trapped | loss at `maxIters = 2^30` | loss at the best guard |
|---|---|---:|---:|---:|
| **32** (`2^28.41`, the live bucket) | table as built | 53.7% | ×8.634 | ×5.328 (`2^28.5`) |
| 32 | table, rule also refusing τ-relations | 13.8% | ×1.773 | ×1.773 (`2^30`) |
| 32 | **σ** | 0 | **×1.185** | ×1.000 (`2^32`) |
| 34 (`2^25.27`, the benchmarks) | table as built | 11.6% | ×4.951 | ×1.655 (`2^27`) |
| 34 | table, rule also refusing τ-relations | 1.8% | ×1.493 | ×1.126 (`2^27.75`) |
| 34 | σ | 0 | ×1.000 | ×1.000 |

These rows use the count.  Scaled to the measured rates (0.86 of the count,
§5), the as-built losses at `dpWeight = 32` are ×7.1 at `2^30` and ×4.8 at
the best guard, and the ranges below span both.  The σ row at
`dpWeight = 32` is the campaign's own configuration today.
`maxIters = 2^30` is `2^1.59` trail lengths, so a fraction
`e^{−3.01} = 4.9%` of honest trails, the longest ones, are cut before their
distinguished point.  That is 15.6% of all steps.

**Paired rates.** These are table over σ, each pair from one session on one
card: `ITERATION-FUNCTION.md` §6 and `ONE-BLOCK-GEOMETRY.md` §1.

| session | σ B/s | table B/s | σ / table |
|---|---:|---:|---:|
| RTX PRO 6000, automatic workers, 6 alternating reps | 14.41 | 16.56 | 0.870 |
| RTX PRO 6000, audited 385k workers, 3 reps | 14.98 | 16.35 | 0.916 |
| RTX PRO 6000, one-block note's session, 256×2 | 14.00 | 15.84 | 0.884 |

The table kernel's 20.08 B/s (`ONE-BLOCK-GEOMETRY.md`) has no σ partner: the
σ kernel was never retuned at 512×1, so it is not a walk comparison and is
not used here.

**Cost to solve, table over σ**, `(c_t/c_σ) × (σ rate / table rate) ×
(table loss / σ loss)`.  Here `c_t/c_σ = 0.928–0.938` is the pooled table
constant over the σ bracket at `n = 131` (§4, `ecc2k130` rows, H = 8), and
each range spans the bracket and the three paired sessions:

| configuration (`dpWeight = 32`) | table / σ | kind |
|---|---:|---|
| as built, both at `maxIters = 2^30` | 4.85–6.26 | measured rates × model loss |
| as built, each at its best guard | 3.85–4.58 | measured rates × model loss |
| rule also refusing τ-relations, best guards | 1.36–1.52 | projection (rule not built) |
| every fruitless cycle caught and escaped, best guards | 0.81–0.86 | projection (check not built) |

**Verdict.**

- The table walk is not a campaign candidate as built.  It costs five to
  six times the σ walk per solve, which is a larger gap than any rate or
  iteration difference in this tree.
- Refusing the τ-relations in the rule is not enough, because the six-step
  cycles remain.
- A walk that detects any cycle and escapes it deterministically would put
  the table walk at 0.81–0.86 of σ's cost.  Two candidate designs:
  - The cheapest widens the rule's history from three tags to five.  At
    H = 8 a tag is 12 bits, so five fit the 64-bit word the device already
    carries.  The rule could then refuse the 24 τ-relation patterns and the
    4 six-step ones outright.  It would still decide from the last few tags
    only, so merged trails would still re-synchronise, and what remains is
    8-step cycles at order `(2n)⁻⁴`.
  - The general one is a periodic comparison with a stored point.  Its
    escape would be a fixed function of the cycle's least point, so that
    two walks in one cycle leave it the same way.

  Either is worth building before the table walk is reconsidered, and this
  harness, with the paired rate, is the test it must pass.
- Independently of the walk, `maxIters = 2^30` at `dpWeight = 32` throws
  away 15.6% of the σ campaign's steps.  `2^32` recovers them.

**Expected work of the σ walk.** On completed trails it is
`2^60.809 × [1.070, 1.082] = 2^60.91–60.92` iterations (extrapolated
bracket, hashed branches).  With the current guard it is ×1.185 more,
`2^61.15–61.17`.  Bailey et al.'s planning budget `2^60.9`, quoted across
this tree, is the first-order `2^60.809 × 1.069`, and it stands.

## 7. Falsification target

Declared before the v2 matrix ran (after the first `n = 23` and `n = 37`
rows): recommend the table walk only if

- `c_t/c_σ < 1` by more than three standard errors at every degree run
  (every degree with eight distinct σ multipliers, `n ≥ 21`; §4 gives the
  structural reason `n = 19` has seven);
- the device rows agree with the emulation within two standard errors at
  `n = 23`;
- the adding control is `1.00 ± 0.01`;
- `c` is independent of `W` at `n = 37`;
- and the cost to solve, with measured paired rates, is below σ's.

Inadmissible: unpaired rates, a changed `dpWeight`, dropping the cycle rule,
or leaving a phase unpriced.  Result:

- The first four conditions hold.
- The last one fails, on a cost the target did not name when it was
  declared: fruitless cycles longer than the rule's window.  That is
  `AGENTS.md` §5's lesson again.  The rate benchmarks run 32k steps per lane,
  under a thousandth of a DP-34 trail, and count a trapped lane's updates
  like any other; the `GF(2^41)` planted logs run 72-step trails.  Nothing
  measured so far could have seen it.

## 8. Classification

| change | class | why |
|---|---|---|
| σ walk constant `1.069` (assumed) against `1.08–1.10` measured at `n = 23–59`, falling towards `1.070` | **accounting** | the planning figure `2^60.9` is confirmed within 1% (`2^60.91–60.92`, extrapolated); no walk changed |
| table walk `c = 1.00` against σ's `[1.070, 1.082]` | **engineering** | 6–7% fewer iterations, at the floor, not across it; a random mapping already does this |
| table walk as built, at campaign settings: rate +9–15%, cost per solve ×4.85–6.26 | **relabelling** | the headline B/s counts the updates of trapped lanes |
| `maxIters = 2^30` at `dpWeight = 32` discards 15.6% of σ steps | **accounting** | a cost the campaign pays that no page priced; raising the guard is the engineering fix |
| v1 table rows → v2 (fruitless 4-cycles no longer scored as collisions) | **accounting** | a correction to this note's own emulation |

None is an advance: no row moves the ratio to the floor below one.

## 9. What this changes elsewhere

- `aws/campaign.json` keeps `"walk": "sigma"`.  *(Round 2 made the
  `maxIters` change below in the repository, §11.)*  Round 1 left the
  `maxIters` change to the campaign owner, for two reasons: the
  bucket copy is what workers read, and `maxIters` is part of the campaign
  contract (`campaignContract` in `aws/protocol.py`), so changing it starts a
  new corpus.  That costs nothing now and costs the corpus later, so the
  change has the same deadline as the walk choice.  Recommended before
  collection: `maxIters ≥ 2^32` at `dpWeight = 32`.  The witness generator
  `src/witness.cpp` has to move with it.  It refuses a trail longer than its
  own `--max-iters`, which defaulted to `2^30`, and a cairn job's
  `max_steps_per_walker` can lower that.  Unless both are raised, the
  distinguished points the longer guard recovers (trails of `2^30` to `2^32`
  steps, about 5% at `dpWeight = 32`) would be refused there.
- [ITERATION-FUNCTION.md](ITERATION-FUNCTION.md) said that choosing the
  table walk later "is a configuration change and not a code change".  That
  did not hold: the walk needed a cycle check first, which §11 builds.  Its §6.2 ratio
  `0.96 ± 0.07` on `GF(2^41)` is consistent with the `0.928–0.938` measured
  here; the `1 + 1/(2H)` it was compared with is the injective model, which
  §3 shows does not apply.
- `README.md` quotes `2^60.9` for ECC2K-130 and now points here.  Its
  "measured for this walk" figure read `2^60.94`, from a draft of §6; round 2
  corrects it to §6's `2^60.91–60.92` (accounting).

## 10. Reproduce

```
cargo build --release --example ecc2k130_walk_constant
target/release/examples/ecc2k130_walk_constant --n 37 --walk table --dist ecc2k130 --branches 8 --walks 16 --trials 40000 --seed 37
make -C ecc2k130 walk-constant          # the device rows at n = 23, H = 8 (about ten minutes on four cores)
cd ecc2k130/benchmarks/walk-constant
python3 summarize.py                    # this note's tables, from the frozen files
python3 fruitless_patterns.py           # the 24 + 4 patterns of §5
python3 trap_cost.py                    # §6's loss table
python3 coset_index.py                  # index 1 at every degree
python3 sigma_classes.py 19 130873 41811 8,2 15000   # σ on classes: why n = 19 is degenerate
python3 mapping_spread.py 45562 16 6000 5             # one-mapping spread of random mappings (W = 8, 16, 64 in the output)
```

Round 2 (§11), from the same directory:

```
./round2.sh                                           # merge, matrix-v3, device-v3 and the residual probe (about two hours on four cores)
python3 fruitless_patterns.py --residual 8            # what rule v2 passes, by signature (fruitless_patterns_v2.txt)
python3 fruitless_patterns.py 8 4 4 --rule v2         # brute force: nothing at four determined tags or fewer (L = 8 runs for hours)
make -C ../.. test-cycle-rule                         # the rule alone: every cycle it must break, and covariance
python3 summarize.py                                  # the round-2 tables, after round 1's
```

`matrix-v3` and `device-v3` repeat `matrix-v2`'s and `device-v2`'s table
rows with the same seeds, under rule v2.

## 11. Round 2: the cycle rule extended

Declared before any code changed (`AGENTS.md` §4), for the change §6's
verdict asks for.

**The change.** The history word already holds the last four step tags in
16-bit slots; the rule reads three.  It will read all four, and refuse a step
`t` when either

- `t` undoes any of the last four steps: `t = −tᵢ`, `i = 1..4`.  This
  replaces both current checks.  Every pairwise cycle of at most 8 steps has
  a pair at most 4 steps apart along the walk, so none can close.  The first
  pairwise survivor has 10 steps (`a, b, c, d, e, −a, …, −e`), at order
  `(2n)⁻⁵`.  The check also refuses some harmless steps (`a, b, −a`), at
  about `4Σp²/2n` per step, by advancing `h` as the rule already does.
- `t, t₁, t₂, t₃` share `h` and form one of the 24 τ-relation patterns of §5.
  That means one repeated `(k, ε)`, with the other two at `(k+1, ε), (k+2, ε)`
  or at `(k+1, −ε), (k+3, −ε)`, offsets mod `n`.

Both tests read only differences of `k` and the parity of `ε`, so the rule
stays a function on classes: `σ` shifts every `k` and negation flips every
`ε`.  The history layout is unchanged, so the device's `& 0xFFFF` readers are
untouched.  `aws/protocol.py` names the rule inside the walk's definition, so
the rule change is a new walk identity.  No table-walk point has been
collected, so nothing forks.

The sketch in §6's verdict, a five-tag history in 12-bit slots, was wrong
about what is needed.  Four full tags are already stored, and cycles longer
than six steps are caught by reading all four.  This is an accounting
correction to §6's sketch, and §6's projection is unaffected by it.

**Target.**  The rule is a success only if all of the following hold.

1. `fruitless_patterns.py` under the new rule finds no pattern with at most
   three determined tags at `L ≤ 8`, and its count at four determined tags
   gives under 0.1% of DP-32 trails trapped at `n = 131`.  That is a count,
   not a measurement: the residual rate at the degrees that run is below
   what `10⁹` steps can see.
2. No τ-relation 4-step return and no pairwise return of 8 steps or fewer
   in either harness, where the old rule showed 1,145 and 126.
3. The table walk's `c` stays within two standard errors of the old rule's
   rows at `n = 23, 37, 41`.
4. Two walks that meet at one point, carrying different histories, part
   within 16 steps no more often than would cost 1% of collisions.  This is
   measured, since a parted merge is a lost collision under distinguished
   points.
5. The device and the reference pick the same tag on every one of 4,096
   points, with histories that trip every new check.  `test-production` and
   the planted small-curve logs pass with `WALK_TABLE=1`, and the CUDA
   sources compile.
6. The switch stays the campaign owner's.  This round prices the table walk
   with the new rule end to end.  The GPU-paired rate of the new kernel is
   the one input it cannot measure here.

Inadmissible: changing `H`, `dpWeight` or the history layout's meaning
without saying so; scoring a formal return as a collision; dropping the
16-point return detector the harnesses use as the independent check.

### 11.1 Results

Sources are the frozen files `matrix-v3.jsonl`, `device-v3.jsonl`,
`merge.jsonl`, `residual-check.jsonl` and `fruitless_patterns_v2.txt` in
[benchmarks/walk-constant](benchmarks/walk-constant/), printed by
`summarize.py`.  [round2.sh](benchmarks/walk-constant/round2.sh) holds the
command lines.  `matrix-v3` and `device-v3` repeat round 1's table rows with
the same seeds.  Where trials are long (`n ≥ 37`) the two rules' walks part
early and a v1/v2 pair is two independent measurements.  At `n = 23` many
trials never reach a step where the rules differ, so the pairs are
correlated: the emulation's differences there are a tenth of the
independent standard error.  The standard errors below treat every pair as
independent, which overstates the uncertainty of the correlated ones.

**1. The count: met.**

- Nothing with four or fewer determined tags survives the rule, at any
  length.  The brute force finds nothing through `L = 7`.  The count by
  signature covers every signature with at most four determined tags, the
  last at `L = 8`, and finds zero in each.  A branch's steps come in even
  numbers, since every power of `τ` is `1 mod τ̄` and `τ̄` has norm 2, so
  there are no other signatures.  A pair needs ten steps, since both of its
  tags must be at least five steps apart around the cycle.
- The first survivors have five determined tags:
  - **1,920 six-step cycles on one branch**, the orders of the nine
    irreducible six-term relations in `Z[τ]`.  An example is
    `τ⁴ + τ³ + τ² − τ − 2 = (τ² + τ + 2)(τ² − 1) = 0`.  None contains a
    pair, and in none do four consecutive steps form a τ-relation.
  - one **10-step pairwise cycle**.
- At six determined tags there are 17,856 interleavings of two τ-relations
  on two branches, and 258 more patterns.  Seven or more are not counted:
  each further determined tag costs a factor of `2n = 262`.
- The brute force over five determined tags finds the same 1,920, and spans
  up to 18 add none.
- At `n = 131`, H = 8, ECC2K-130's distribution, the survivors are entered at
  `5.7 × 10⁻¹⁴` per step, 99.6% of it through the six-step cycles.  That traps
  `2.0 × 10⁻⁵` of DP-32 trails, fifty times under the target's 0.1%.

**2. Returns: met as declared, and the declaration missed a class.**

- **What v2 had to remove.** There were no τ-relation 4-step returns and no
  pairwise returns in 2.36 × 10⁹ table steps across both harnesses.  Rule v1's
  count predicts 1150 on the same rows.
- **What it leaves.** One return was of a kind item 2 does not name: a
  six-step τ-relation return at `n = 23`, native, H = 16.  The residual count
  of item 1 predicts 0.81 there and 2.13 over all rows.  Item 2 listed
  what the rule must remove and not what it leaves.  Item 1's count found the
  class, and the probe below tests it.
- **The probe.** At `n = 23` with two uniform branches, where the six-step
  cycles are entered about a thousand times more often than at H = 8, 52.6M steps show 21
  six-step and 3 eight-step τ-relation returns.  The count gives 15.3 and 1.5.
- **Why the probe is only a bound.**
  - At H = 2 the rule can refuse both branches, and the step is then taken
    anyway.  The probe's 15,646 two-step and 58 four-step pairwise returns,
    and its 25 four-step τ-relation returns, are such steps.
  - At H = 8 this cannot happen.  The rule refuses at most five branches:
    four by negation, one per history slot, and one by the τ test, which
    needs the last three tags on one branch.
  - The count does not model the forced steps, so the probe's excess of 1.4×
    bounds the count rather than checks it.  At twice the count the trap loss
    would be ×1.0006, and item 6's projection would move in its fourth digit.

**3. `c`: not met as declared; no shift pooled.**

| n | branches | H | harness | W | c (v2) | ± | c (v1) | ± | difference / SE |
|---:|---|---:|---|---:|---:|---:|---:|---:|---:|
| 23 | ecc2k130 | 8 | device | 8 | 1.0121 | 0.0037 | 1.0076 | 0.0037 | +0.9 |
| 23 | ecc2k130 | 8 | emulation | 8 | 1.0038 | 0.0012 | 1.0037 | 0.0012 | +0.0 |
| 23 | ecc2k130 | 16 | emulation | 8 | 1.0019 | 0.0012 | 1.0021 | 0.0012 | -0.1 |
| 23 | native | 8 | device | 8 | 1.0016 | 0.0021 | 1.0080 | 0.0021 | -2.1 |
| 23 | native | 8 | emulation | 8 | 1.0044 | 0.0012 | 1.0050 | 0.0012 | -0.4 |
| 23 | native | 16 | emulation | 8 | 1.0055 | 0.0012 | 1.0046 | 0.0012 | +0.6 |
| 23 | uniform | 8 | emulation | 8 | 1.0027 | 0.0012 | 1.0026 | 0.0012 | +0.1 |
| 37 | ecc2k130 | 8 | emulation | 16 | 1.0032 | 0.0026 | 1.0036 | 0.0026 | -0.1 |
| 37 | native | 8 | emulation | 16 | 1.0032 | 0.0026 | 1.0015 | 0.0026 | +0.5 |
| 41 | ecc2k130 | 8 | device | 16 | 1.0638 | 0.0302 | 1.0051 | 0.0299 | +1.4 |
| 41 | ecc2k130 | 8 | emulation | 16 | 0.9906 | 0.0058 | 1.0125 | 0.0059 | -2.6 |
| 41 | native | 8 | emulation | 16 | 1.0060 | 0.0059 | 1.0084 | 0.0059 | -0.3 |
| 59 | ecc2k130 | 8 | emulation | 16 | 1.0031 | 0.0037 | 0.9949 | 0.0037 | +1.6 |

- **Pooled.** Over the 13 matched rows, `c(v2) − c(v1) = -0.0002 ± 0.0007`.
- **Per row.** Two rows miss the declared two standard errors, both on the
  low side:
  - the device's `n = 23` native row (−2.1: 1.0016 against v1's 1.0080);
  - the emulation's `n = 41` ECC2K-130 row (−2.6: 0.9906 against 1.0125).

  In both, the v1 row is the outlier against its model (1.0080 against
  1.0026, 1.0125 against 1.0008).  The v2 rows sit within two standard errors
  of theirs.
- **Chance.** If the 13 comparisons were independent, they would miss
  two standard errors 2 or more times with probability 0.12.
- **Device against emulation under v2.**
  - At `n = 23`: +2.1 standard errors for ECC2K-130's distribution (1.0121
    against 1.0038) and −1.1 for native.
  - At `n = 41`: +2.4 for ECC2K-130's distribution (1.0638 ± 0.0302 over 300
    trials, against 0.9906).
  - Under v1 the same comparisons were +1.0, +1.2 and −0.2.
  - Round 1 declared device/emulation agreement at `n = 23` as a control
    (§7), so a replicate was declared before it ran.  It repeats the device's
    `n = 23` ECC2K-130 row with 40,000 trials on a new seed (231).  If it
    also sits more than two standard errors above the emulation, the two
    harnesses disagree under v2, and the pooled `c` below is not used until
    that is explained.
- **The lesson.** A per-row test was the wrong one to declare for a change
  expected to do nothing across a dozen rows.  The pooled difference is the
  right test.
- **The pooled constant.** The table walk under v2 is at the floor, as under
  v1: `c = 1.0039 ± 0.0010` pooled over the ECC2K-130 rows at H = 8.

**4. Merge parting: met, narrowly.**

- Two walks that meet at one point with different histories part when the
  rule refuses a step for one and not the other.
- Under v2 the first step after the meeting has four history slots still
  holding different tags, the next three, then two, then one.  Each slot
  refuses with probability `q = Σp²/2n`, for either walk.  That predicts
  `2(4 + 3 + 2 + 1)q = 20q`.  Under v1, only the step that undoes the last
  one (at first order), so `2q`.
- Measured over nine v2 rows, including the device's own walk at `n = 23`
  and `41`: 0.93–1.00 of `20q`, parted in the proportions 4 : 3 : 2 : 1 over
  the first four steps.  For v1, 0.98–1.10 of `2q`.

| n | branches | rule | harness | merges | parted | rate | ± | 20q or 2q |
|---:|---|---|---|---:|---:|---:|---:|---:|
| 23 | ecc2k130 | v1 | emulation | 399,987 | 2,227 | 5.57e-03 | 1.2e-04 | 5.50e-03 |
| 23 | ecc2k130 | v2 | device | 40,000 | 2,209 | 5.52e-02 | 1.1e-03 | 5.50e-02 |
| 23 | ecc2k130 | v2 | emulation | 399,989 | 21,438 | 5.36e-02 | 3.6e-04 | 5.50e-02 |
| 23 | native | v1 | emulation | 399,995 | 4,054 | 1.01e-02 | 1.6e-04 | 1.02e-02 |
| 23 | native | v2 | device | 40,000 | 3,799 | 9.50e-02 | 1.5e-03 | 1.02e-01 |
| 23 | native | v2 | emulation | 399,996 | 38,111 | 9.53e-02 | 4.6e-04 | 1.02e-01 |
| 41 | ecc2k130 | v1 | emulation | 400,000 | 1,355 | 3.39e-03 | 9.2e-05 | 3.09e-03 |
| 41 | ecc2k130 | v2 | device | 4,000 | 122 | 3.05e-02 | 2.7e-03 | 3.08e-02 |
| 41 | ecc2k130 | v2 | emulation | 400,000 | 12,313 | 3.08e-02 | 2.7e-04 | 3.09e-02 |
| 41 | native | v1 | emulation | 400,000 | 1,755 | 4.39e-03 | 1.0e-04 | 4.30e-03 |
| 41 | native | v2 | emulation | 400,000 | 16,599 | 4.15e-02 | 3.2e-04 | 4.30e-02 |
| 59 | ecc2k130 | v1 | emulation | 400,000 | 889 | 2.22e-03 | 7.4e-05 | 2.15e-03 |
| 59 | ecc2k130 | v2 | emulation | 400,000 | 8,528 | 2.13e-02 | 2.3e-04 | 2.15e-02 |
| 59 | native | v1 | emulation | 400,000 | 1,002 | 2.50e-03 | 7.9e-05 | 2.55e-03 |
| 59 | native | v2 | emulation | 400,000 | 9,966 | 2.49e-02 | 2.5e-04 | 2.55e-02 |

- At `n = 131`, H = 8: 0.97% of merges part under v2, against 0.10% under
  v1.  A parted merge is a lost collision, so the expected number of
  collisions a solve needs grows by `1/(1 − δ)`, and its iterations by
  `1/√(1 − δ) = ×1.0048`.
- The 1% limit holds by a margin of 3%.  The rule costs this because
  refusing more patterns means refusing on more of each walk's history.

**5. Agreement and tests: met.**

- `make test-cycle-rule` runs the rule alone, from random histories, at
  `m = 23, 41, 131`, and every check passes:
  - all 19,200 τ-relation 4-cycles (both relations, all 12 orders each,
    phases across the wrap) are refused within the first lap;
  - every pairwise cycle of 2, 4, 6 and 8 steps, in all its matchings, is
    refused within the first lap;
  - the 10-step one is refused in 2–7 of 200 laps, only by accidental
    matches with the history;
  - an empty history never fires;
  - the τ test fires 0 times on 200,000 random quadruples;
  - conjugating every tag by `σᶜ` or negating it changes the rule's answer
    0 times in 200,000 per degree, about 67,000 of which fire.
- **Device against reference.**
  - The host probe (`make test-table-walk-host`), in both pivot layouts,
    finds 0 mismatches in phase, pivot, sign, tag and addend words on
    4,096 points.  The new histories trip the rule on 2,560 of them.
  - The CUDA probe, `testtablewalkcuda.cu`, carries the same new
    histories.  It compiles here with `nvcc` 12.9, but no GPU runs it.
  - The client's covariance self-test adds a τ-relation history.
- **Other tests.**
  - `test-production` passes in both builds.  The table build needed a new
    collision fixture, since the rule change moves trajectories;
    `test-production --find-fixture` found one that resolves to the same
    discrete log.
  - `certify-local`, the planted small-curve logs (8 of 8) and the table
    client's `--test` pass.
  - The CUDA sources compile, including the 20b-knob preset with
    `PACKED_CLMAD=0`; its `PACKED_CLMAD` needs CUDA 13.3.

**6. Cost to solve: a projection.**

- The inputs, at `dpWeight = 32` with both walks at `maxIters = 2^32`:
  - the pooled v2 constant over the σ bracket `[1.070, 1.082]`;
  - the merge factor ×1.0048;
  - the table walk's trap loss ×1.00031 from item 1's count, against σ's
    guard loss ×1.00007;
  - the v1 kernel's three paired rate sessions (§6).
- Together they give:

| configuration (`dpWeight = 32`, `maxIters = 2^32`) | table / σ | kind |
|---|---:|---|
| rule v2 (this round) | 0.81–0.86 | projection: the v1 kernel's paired rates |
| every fruitless cycle caught and escaped (§6) | 0.81–0.86 | projection (§6) |
| rule v1, both at `maxIters = 2^30` (§6) | 4.85–6.26 | measured rates × model loss (§6) |

- **What the projection does not know.** The v2 kernel's own rate on the
  card is unmeasured.  Its rule reads four history slots instead of three
  and adds the τ test, a few integer comparisons per step.  It fires about
  four times as often as v1's, `4q` against `q` per step (0.43% against 0.11%
  at `n = 59`, about 0.2% at `n = 131`), and each firing repeats the check on
  the next branch.
  None of that touches the point arithmetic, but none of it is measured.  If the v2
  kernel loses less than 14–19% of the v1 kernel's rate, the table walk costs
  less per solve than σ.
- **Before it is chosen,** that paired measurement is required: v2 table
  kernel against the σ kernel, alternating, in one session on one card, as
  §6's rows are.

### 11.2 Classification

| change | class | why |
|---|---|---|
| rule v1 → v2 | **engineering** | `c` is at the floor before and after (pooled difference `-0.0002 ± 0.0007`); the cost per solve against σ falls from 4.85–6.26 to 0.81–0.86 (projection) by removing a loss, not by crossing the floor |
| `maxIters` `2^30` → `2^32`, witness default with it | **engineering** | σ's guard loss ×1.185 → ×1.00007; expected work on the campaign's configuration `2^61.15–61.17` → `2^60.91–60.92` (§6 bracket) |
| §6's five-tag sketch → four 16-bit tags | **accounting** | a correction to this note's own sketch; nothing measured changes |
| the residual six-step relations | **accounting** | §6's 0.81–0.86 assumed every cycle caught; what is built leaves `2.0 × 10⁻⁵` of trails trapped, a ×1.0002 difference, so the projection is unchanged to two digits |
| `README.md`'s `2^60.94` → `2^60.91–60.92` | **accounting** | a figure from a draft of §6 left in the README |

None is an advance: no row moves the ratio to the floor below one.

### 11.3 What changes elsewhere

- `aws/campaign.json`:
  - `maxIters` is `2^32`, and `src/witness.cpp`'s default `--max-iters` is
    `2^32` with it.  The bucket copy is still what workers read, and
    republishing it is the campaign owner's.
  - `maxIters` is in the campaign contract, so the change is free only
    before the first distinguished point.  `docs/ecc2k130-status/history.json`
    is still empty.
  - Every counter the guard compares, and every record field, is 64-bit.
    The witness's per-branch step counts are 32-bit.  A trail of at most
    `2^32` steps plus one guard period overflows one only by taking nearly
    every step on one branch, and the record's own sum check catches a
    wrapped count.
  - `"walk"` stays `"sigma"`.
- `aws/protocol.py` names the rule inside the table walk's definition, so
  the table walk under v2 is a new walk identity.  No table-walk point
  exists, so nothing forks.
- Switching needs three things, in order:
  1. the paired rate above;
  2. the owner's decision;
  3. both before the first distinguished point, since the walk and
     `maxIters` are both in the contract.
- [ITERATION-FUNCTION.md](ITERATION-FUNCTION.md) and the [README](README.md)
  now point here, not at round 1's verdict.
