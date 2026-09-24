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

Answer: **keep the σ walk, and raise `maxIters`.**

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
  emulation, and matches a count of the patterns: 1145 τ-relation
  returns seen against 1297 predicted.
- **What the cycles cost.** At the live bucket's `dpWeight = 32` those
  cycles trap **54%** of the table walk's trails (H = 8).  Each trapped trail
  runs to the `maxIters = 2^30` guard and is discarded, so a completed trail
  costs **8.6×** its steps.  Priced end to end, the table walk as built costs
  **5.88–6.26× the σ walk per solve**.  With every fruitless cycle caught
  and escaped (a check that does not exist yet), it would cost
  **0.81–0.86×**; that figure is a projection.
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
classification, and §9 what this changes elsewhere.

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
- The coset index of `⟨1 + λʲ, λ, −1⟩` in `Z_ℓ^*`, which would confine the
  σ walk to a coset and cost it a factor no branch statistic shows, is **1**
  at every degree run and at ECC2K-130 ([coset_index.txt](benchmarks/walk-constant/coset_index.txt)).

## 3. Why the walks differ

The σ walk multiplies the scalar by `1 + λʲ`.  On one branch that is a
bijection of the class set, so two classes on the same branch never
collide, and the walk pays the full injective penalty
`1/√(1 − Σp²) = 1.070`.  Because the multipliers commute it also pays a
higher-order penalty: `(1 + λᵃ)(1 + λᵇ)` is the same step in either order,
which forbids the two-step collisions a random mapping would have.  The
emulation measures that excess at `c/model = 1.016–1.028` with hashed
branches (§4).

The table walk adds `ε σᵏ(T_h)` with `k, ε` read from the class.  Two
classes on one branch land in the same class whenever their frames differ,
so on the class set it is not injective per branch.  It behaves like a
random mapping up to `Σp²/2n`.

## 4. The single table

All rows are in unit `c`; the model column is §1's model for that walk.
Sources are the frozen files in [benchmarks/walk-constant](benchmarks/walk-constant/),
printed by `summarize.py`.  σ rows come from `matrix-v1.jsonl`, whose σ code
is unchanged.  v1's table and adding rows used only a previous-class cycle
rule and scored fruitless 4-cycles as collisions.  At `n = 41`, with 102k
steps per trial, that ended about 14% of trials early and gave
`c = 0.909 ± 0.006`, which is the prediction for that artefact
(`1 − 0.64 λμ` with `λ = q²`).  Those rows are superseded, not deleted.

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

- The table walk is at the floor: `c = 1.002–1.008`, with no dependence on
  `H`, on the branch distribution (`Σp²` from 0.10 to 0.23) or on `W`.
  The device and the emulation agree.
- The σ walk with ECC2K-130's distribution is at `c = 1.0995, 1.0964, 1.0919, 1.0910, 1.0817` at
  `n = 23, 31, 37, 41, 59`.  That is the first-order `1.070` times a
  higher-order excess (the commuting multipliers of §3) which shrinks with
  the degree, from 2.8% to 1.1%.  The first-order value is its limit, so at
  `n = 131` the constant lies in `[1.070, 1.082]`.  That is an
  extrapolation, and §6 carries the whole bracket.
- `n = 19` is excluded for σ.  Two of its eight multipliers fall in one
  class (`1 + λ¹⁰ = λ⁻⁹(1 + λ⁹)` when `n = 19`; every `n ≥ 21` has eight
  distinct ones), and with 3,444 classes its short walks return to their
  own classes early.  A pure-arithmetic model of the σ walk on classes
  ([sigma_classes.py](benchmarks/walk-constant/sigma_classes.py),
  [output](benchmarks/walk-constant/sigma_classes.txt)) reproduces the
  `1.006`, and gives `1.087 ± 0.011` at `n = 23`: a third implementation,
  agreeing with the other two.
- One mapping at a time, both walks' constants stray from the average over
  mappings far more than their standard errors allow: σ 1.108 with a standard deviation of 0.034 over 6 mappings, each ± 0.003; table 0.995 with a standard deviation of 0.041 over 3 mappings, each ± 0.003.  A random
  mapping on 3.1M classes would stray by about 0.1%.  The suspected cause is
  accidental short relations mod `ℓ` among a mapping's steps, since
  `(2Hn)³ ≈ ℓ` at `n = 37`; at ECC2K-130's `ℓ ≈ 2^129` no short relation
  can be accidental.  At `n = 59` (`ℓ ≈ 2^33`) the spread is (running).
- With the device's own weight branches, σ is at `1.122–1.217` at
  `n = 23–59`, with `c/model` 1.010–1.099.  Each of those rows is a single
  mapping, and the scatter, including `n = 37`'s 1.217, is within the
  one-mapping spread just described.  So these rows establish no effect of
  weight-based branches beyond the model column's `Σp²`.
- §6 uses the hashed σ bracket.  That is the conservative choice for
  switching, since any extra cost of the native branches would only make σ
  dearer.

| n | walk | mappings (seed: c) | mean | spread (sd) | per-mapping ± |
|---:|---|---|---:|---:|---:|
| 37 | sigma ecc2k130 H = 8 | 1: 1.0810, 2: 1.0752, 3: 1.1057, 4: 1.0960, 5: 1.1680, 6: 1.1239 | 1.1083 | 0.0341 | 0.0029 |
| 37 | table ecc2k130 H = 8 | 1: 0.9810, 2: 1.0412, 3: 0.9630 | 0.9951 | 0.0409 | 0.0026 |

## 5. Fruitless cycles

The table walk is additive, so a run of steps whose addends sum to zero
returns a walk to a point it has left, with no collision.  The device's rule
(`eccTagFruitless` in `include/tablewalk.h`) refuses two such runs: a step
that undoes the last one, and a pairwise 4-cycle.  Enumerating the patterns
it passes ([fruitless_patterns.py](benchmarks/walk-constant/fruitless_patterns.py),
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
  emulation, across a factor of 343 in predicted rate: 1145 seen against 1297 predicted, 0.88 ± 0.03 overall.
- Where it falls short, as at `n = 19`, the walks are short.  A `W = 8`
  trial there ends after about nine steps per walk, which cuts off windows
  at both ends.  The ratio is 425 seen against 488 predicted, 0.87 ± 0.04 on the rows at `n ≥ 37`, where walks
  run thousands of steps.
- The pairwise count is small and noisy: 126 seen against 164 predicted, 0.77 ± 0.07 overall, 64 seen against 65 predicted, 0.99 ± 0.12 at
  `n ≥ 37`.
- The pairwise term is 14% of the `n = 131` rate, so the pricing below
  rests mainly on the τ-relation law, which the device walk confirms
  directly: 65 seen against 70 predicted at `n = 23`.

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

The σ row at `dpWeight = 32` is the campaign's own configuration today.
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
| as built, both at `maxIters = 2^30` | 5.88–6.26 | measured rates × model loss |
| as built, each at its best guard | 4.30–4.58 | measured rates × model loss |
| rule also refusing τ-relations, best guards | 1.43–1.52 | projection (rule not built) |
| every fruitless cycle caught and escaped, best guards | 0.81–0.86 | projection (check not built) |

**Verdict.**

- The table walk is not a campaign candidate as built.  It costs about six
  times the σ walk per solve, which is a larger gap than any rate or
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

- `c_t/c_σ < 1` by more than three standard errors at every degree run;
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
| table walk as built, at campaign settings: rate +9–15%, cost per solve ×5.88–6.26 | **relabelling** | the headline B/s counts the updates of trapped lanes |
| `maxIters = 2^30` at `dpWeight = 32` discards 15.6% of σ steps | **accounting** | a cost the campaign pays that no page priced; raising the guard is the engineering fix |
| v1 table rows → v2 (fruitless 4-cycles no longer scored as collisions) | **accounting** | a correction to this note's own emulation |

None is an advance: no row moves the ratio to the floor below one.

## 9. What this changes elsewhere

- `aws/campaign.json` is unchanged: `"walk": "sigma"` stays.  The
  `maxIters` change is the campaign owner's to make, for two reasons: the
  bucket copy is what workers read, and `maxIters` is part of the campaign
  contract (`campaignContract` in `aws/protocol.py`), so changing it starts a
  new corpus.  That costs nothing now and costs the corpus later, so the
  change has the same deadline as the walk choice.  Recommended before
  collection: `maxIters ≥ 2^32` at `dpWeight = 32`.
- [ITERATION-FUNCTION.md](ITERATION-FUNCTION.md) said that choosing the
  table walk later "is a configuration change and not a code change".  That
  no longer holds: the walk needs a cycle check first.  Its §6.2 ratio
  `0.96 ± 0.07` on `GF(2^41)` is consistent with the `0.928–0.938` measured
  here; the `1 + 1/(2H)` it was compared with is the injective model, which
  §3 shows does not apply.
- `README.md` quotes `2^60.9` for ECC2K-130 and now points here.

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
```

Every row names its seed.  The full matrices took about an hour on four
cores: `matrix-v2` in the emulation, and `device-v2` on the reference walks.
