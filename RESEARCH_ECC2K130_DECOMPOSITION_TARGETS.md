# Targets: six experiments on ECC2K-130 point decomposition

**Boundaries:** `scripts/ecc2k130_decomposition_targets.py`
**Frozen artefact:** `experiments/ecc2k130_decomposition_targets.json`
**Background:** [`RESEARCH_ECC2K130_DECOMPOSITION.md`](RESEARCH_ECC2K130_DECOMPOSITION.md)
**Pipeline under test:** `src/cryptanalysis/koblitz_index_calculus.rs`,
`src/cryptanalysis/semaev_decomp.rs`, the `ic` tool

Six scoped experiments, each with a boundary derived *before* anything is run, a
primary metric in one unit, and a falsifier specific enough that a run either
meets it or does not.  Written to be picked up one at a time.

The background note answered "can an ECC2K-130 point be decomposed?" with
*exists, admissible, not findable*, and closed the fix-some-summands oracle
family at `m·2^131`.  These are the six places that answer could still be wrong,
ordered by how much would change if it were.

**One of them already paid.**  E6 was written as a check on a line in the
background note and instead corrected it: the GGMP collapse needs `π(F) = F`,
not a subspace, so it *is* available at `n = 131` from a union of Frobenius
orbits, and the bottom of the family is `2^124.99` rather than `2^132.58`.  That
is `2^64.18×` rho instead of `2^71.77×` — an **accounting** correction by §3 of
`AGENTS.md`, recorded in §6.1 of the background note, and the reason the rest of
this list is worth running rather than assuming.

## The shared boundary

Everything is scored against Pollard rho with the `⟨−1⟩ × ⟨π⟩` speed-up on the
same subgroup:

```text
rho reference at n = 131:   2^60.8090        S = 0.077430
```

and the family's derived line, from the background note's §5.1:

```text
    2^l relations · 2^n/C(|F|,m) targets · C(|F|,m−1) oracle  =  m · 2^n
                                    ... and m · 2^n / n for a Frobenius-stable base.
```

Every experiment below reports **`Λ = total operations / 2^n`**, which is the
unit that makes the product law a constant (`Λ = m`, or `m/n`) and turns "did
anything move" into reading one column.  `S = ops/√r` is reported alongside so
the rows drop into `docs/index-calculus-scoreboard.html` unchanged.

---

## E1 — The scale model: a ladder of curves shaped like ECC2K-130

**Question.**  Everything at `n = 131` in the background note is derived.  Does
the product law actually hold on curves with the same structure?

**The ladder.**  Prime `n` where `2` is a primitive root, so `ord_n(2) = n−1`,
the 2-cyclotomic cosets are `{1, n−1}`, and the only invariant subspace
dimensions are `0, 1, n−1, n` — the exact obstruction that defines the ECC2K-130
instance.  Below 70 that is

```text
    n = 11, 13, 19, 29, 37, 53, 59, 61, 67          (then 83, 101, 107, 131)
```

on `K_0 : y² + xy = x³ + 1`, with a subspace factor base of dimension
`l = ⌈(n + log₂ m!)/m⌉` and `m = 3`.

**Feasibility, and the honesty line it forces.**  The predicted total is `3·2^n`,
so end-to-end runs stop at `n ≈ 37` (`2^38.6` group operations, hours) and the
rest must be **composed** from separately measured phase rates.  Report them in
different columns and never in the same one:

| class | rungs | what is measured |
|---|---|---|
| end-to-end | 11, 13, 19, 29, 37 | one planted logarithm recovered per run |
| composed | 53, 59, 61, 67 | measured yield, measured oracle rate, measured matrix, multiplied |

The composed rungs exist to extend the fit, and the end-to-end rungs exist to
validate the composition: at `n = 29` and `n = 37` the composed projection and
the real run must agree to within `1.3×` or the composition is not admissible
evidence.

**A distortion to control for.**  The odd part of `#E(F_2^n)` is not prime at
every rung — `n = 11` is `23²`, `n = 67` has a 26-bit largest prime factor.
Pose the discrete logarithm in the largest prime-order subgroup, report
`log₂` of it, and use `Λ` (which does not see the subgroup) as the primary
metric; `S` is the secondary and must be read with the subgroup column beside
it.  The faithful rungs, where the largest prime is within 4 bits of `2^n`, are
`n = 13, 19, 59`.

**Primary metric.**  `Λ = total operations / 2^n`, predicted flat at `m = 3`.

**Falsifier.**  A least-squares slope of `log₂(total)` against `n` below `0.95`
over four or more rungs, or any rung with `Λ < 0.5·m`.  Either says the product
law is not what governs these curves and the `n = 131` derivation does not
follow.

**Correctness gate.**  Every end-to-end rung recovers its planted logarithm;
every relation is re-verified in the group before entering the matrix; zero
trivial relations counted.  A run failing any of these is a bug report.

**What it costs.**  A few CPU-days for the end-to-end half.  This is the
experiment that most directly tests the background note, and the one to run
first.

---

## E2 — Is the total really flat in the factor-base dimension?

**Question.**  The product law's whole content is that `l` is not a lever.  The
note measures that as arithmetic; measure it as a run.

**Design.**  One curve from E1's ladder (`n = 59` for the faithful subgroup, or
`n = 29` if the budget is tight), `m = 3`, sweep `dim V = 8, 12, 16, 20, 24, 28`
against a saturating dimension of `20.53`.  Measure targets tried, oracle
seconds per target, relations, and linear algebra separately, and compose.

**Predicted**, from the frozen artefact: `Λ` constant up to the saturating
dimension, then rising as `2^{ml−n}`.

**Falsifier.**  Any dimension whose measured total is below half the flat line.
That would mean the cancellation between "fewer targets" and "dearer targets" is
not exact, and the exponent at 131 is not `n`.

**What it costs.**  Hours, on top of E1's harness.  Cheap, and it is the
experiment that would catch a modelling error in E1 before it propagates.

---

## E3 — What a free detector buys, and why deciding is localising

**Question.**  The background note's free-oracle floor grants a hypothetical
oracle *everything* — it decides, and hands over the witness — and comes out at
`2^56.40` for `m = 4`, **below rho**.  That is what says the counting argument
does not protect this curve.  It is tempting to hope that an oracle answering
only **yes or no** is a much weaker object, and that the gap is where the real
protection lives.  This section used to say exactly that.  It was wrong twice,
and the second correction is the interesting one.

**It is not a weaker object.**  §3.2 of the note measures the reduction: let `R`
decompose, walk the candidates `P ∈ F`, draw a base point `Q` of the same
cofactor class, and ask the detector about `R − P + Q`.  The class match keeps
the query inside `⟨G⟩`, so even a subgroup-only detector suffices.  If `P` is a
summand then `R − P + Q` is literally an `m`-sum of base points and the answer
is yes by construction; if it is not, the point is generic and the answer is yes
with probability `λ`.  Two independent `Q` per candidate, and the summands fall
out of `2|F|` **whole-base** queries — no sub-base query anywhere.  The queries
are free by hypothesis; the group operations that build them cost `k·2^{2l}`
against the `m·2^{2l}` the linear algebra already pays, so they are absorbed.

So the column to compare is not "decides" against "localises" but
**target-agnostic** against **target-restricted**: whether the detector is an
algorithm on the target's coordinates, or answers only for some distinguished
family of targets and refuses `R − P + Q`.

| `m` | free target-agnostic detector | witness route | free target-restricted detector | gap | beats rho? |
|---:|---:|---|---:|---:|---|
| 2 | `2^89.45` | table | `2^89.45` | `2^0.00` | neither |
| 3 | `2^68.88` | table | `2^68.88` | `2^0.00` | neither |
| 4 | `2^56.76` | swap | `2^68.29` | `2^11.53` | agnostic only |
| 5 | `2^48.76` | swap | `2^59.75` | `2^10.99` | **both** |
| 6 | `2^43.15` | swap | `2^62.00` | `2^18.85` | agnostic only |
| 7 | `2^39.01` | swap | `2^56.93` | `2^17.92` | **both** |
| 8 | `2^35.86` | swap | `2^59.27` | `2^23.41` | **both** |

An agnostic detector may still build the restricted one's table, so it takes
whichever witness route is cheaper: the meet-in-the-middle table at `m = 2, 3`,
where the two columns coincide exactly, and the swap from `m = 4` up, where it
costs `0.2`–`0.35` bits over the free-oracle floor of §5.3.

**The agnostic column is under rho from `m = 4`, monotonically, with no table.**
The restricted column saws — it dips under at `m ∈ {5, 7, 8}` and back above at
`6`, because its `⌈m/2⌉ / ⌊m/2⌋` meet-in-the-middle splits more evenly at odd
`m` — and every restricted cell that reaches rho carries a `2^45`–`2^57`-entry
table, which at `m = 5` is most of its total.  But the gap between the columns
measures an **artificial restriction**, not the curve: no proposed detector has
that shape.  A resultant-vanishing test, a trace condition, a partial Gröbner
refutation at fixed degree, the Nagao/Riemann–Roch coefficient search are all
algorithms on the target's coordinates and all swap for free.

*(Two superseded numbers, kept so the correction is legible.  The first revision
charged the witness a naive `C(|F|, m−1)` search and reported the second column
bottoming out at `2^73.29`, "above rho at every `m`" — wrong by up to sixteen
bits.  The second revision fixed the price but kept the premise, and reported a
`2^11`–`2^24` gap as the curve's protective margin.  The margin is not there.)*

**Design.**  Two measurements, both at toy sizes on E1's harness:

1. For each oracle the repository has — pairs-and-solve, matrix-F4, SAT —
   measure the cost of a query on `R − P + Q` against a query on `R`. Predicted:
   equal, because all three take the target as input like any other. Anything
   that is *not* equal is the finding.
2. Take any *proposed* cheap detector and measure the same ratio, plus the two
   swap failure rates against the collision law of §3.2 (`miss ≈ k·m/|F|`,
   `false positive ≈ (m/|F| + λ)^k`). A detector that swaps for free is scored
   against §5.3's free-oracle floor directly, cost and memory together.

The Nagao/Riemann–Roch encoding landed on main by #316
(`research/nagao_relations/`) is the natural second case.  It is a different
*relation representation* rather than a subspace-membership oracle — it searches
for a function in `L(4O)` and reads the decomposition off its zeroes — so the
question lands on it unchanged: does its coefficient search cost the same on
`R − P + Q` as on `R`?  That note is explicit that it is "implementation and toy
correctness work, not an ECDLP speedup result", so nothing here is scored
against it yet; E3 is the frame it would be scored in.

**Falsifier / what would count.**  A proposed detector whose cost on `R − P + Q`
exceeds its cost on `R` by more than a constant — that, and only that, puts the
restricted column back in play.  Failing that, the thing to exhibit is the
detector itself: a full-base query cheaper than `C(|F|, m−1)`, scored against
`2^56.40` with its memory priced against BSGS at the same memory.

**Why this one is worth designing carefully.**  "A fast decomposition oracle"
covers three different objects, and two revisions of this section mispriced the
distinctions between them.  What survives is the shortest of the three claims
and the only one that holds: **for any detector that takes its target as input,
deciding and localising are the same problem**, so there is one free-oracle floor
and it is the one in §5.3 of the note.

---

## E4 — Do large primes move the product law, or relabel it?

**Question.**  Large primes are the one classical index-calculus lever the
background note does not price.  Allow the last summand to land anywhere in a
larger subspace `V' ⊃ V`, store the partial relation under that large prime, and
pair partials off by birthday.

**The model, and its guard.**  With `dim V = l`, `dim V' = l'`:

```text
    partial relations per target  =  C(2^l, m−1) · 2^{l'} / (m · 2^n)
    partials needed               =  2^{(l + l' + 1)/2}            (to yield 2^l full ones)
    collection                    =  m · 2^{n + (l − l' + 1)/2}
    memory                        =  2^{(l + l' + 1)/2}
```

The guard is **partial yield ≤ 1 per target**.  Without it `l' → n`, the oracle
stops filtering anything, and what is left is a birthday search on differences of
targets — a generic algorithm, which has to be scored as one and loses to rho.
Guarded, the derived optima are `2^70.50` at `m = 2` (store `2^66.50`), `2^72.31`
at `m = 3`, `2^73.73` at `m = 4`; memory-capped they run

| store | large-prime variant | BSGS, same memory |
|---:|---:|---:|
| `2^30` | `2^107.58` | `2^99.00` |
| `2^40` | `2^97.58` | `2^89.00` |
| `2^50` | `2^87.58` | `2^79.00` |
| `2^60` | `2^77.58` | `2^69.00` |

**Prediction.**  A genuine and large improvement on `m·2^131` — sixty bits at the
unguarded optimum — that nonetheless loses to plain baby-step giant-step at every
memory budget by a near-constant `2^8.6`, and to rho by more.  Class:
**relabelling**, by §3 of `AGENTS.md`.

**Falsifier.**  A guarded cell below the BSGS line at the same memory.

**Caveat stated in advance.**  The model above is crude on purpose: it does not
price collision bookkeeping, and it does not check that the full relations
produced by pairing are **independent** over the factor base.  Independence is
the most likely place for it to be wrong, and the experiment must measure the
rank of the relation matrix, not assume it.

---

## E5 — The yield distribution, not just its mean

**Question.**  Every extrapolation in the background note rests on
`λ = C(|F|,m)/#E` and a Poisson tail.  The twelve measured cells came in at
`0.95` to `1.38` times the prediction, eleven of them above `1.0`.  That is the
signature of **under-dispersion** — the decompositions spread more evenly across
targets than independence would give — and if it is real it is a small
security-negative correction that nobody has written down.

**Design.**  On the toy rungs where the full curve is enumerable (`n ≤ 19`),
compute the **exact** decomposition count for every target rather than a yes/no,
and test the index of dispersion `Var/mean` against its null value of `1`.  Its
standard error is `√(2/T)`, so resolving a 20% effect at 3σ needs **450 targets
per cell** — well inside budget.

**Primary metric.**  `Var/mean` per `(n, m, l)` cell.

**Falsifier.**  `|Var/mean − 1| > 0.2` at 3σ on three or more cells, with a
consistent sign.  That would mean the Poisson tail is the wrong model and the
`m·l ≥ n + log₂ m!` saturation threshold needs a correction term.

**Bound on what it can matter.**  A constant, not an exponent: the mean is what
the product law consumes, and E5 can only move the constant in front.  It is on
this list because it is cheap and because the background note's central table
quietly depends on it.

---

## E6 — Frobenius-stable bases without a subspace anywhere

**Status: this one already fired, and corrected the background note.**

**Question.**  The classification "the only invariant dimensions at `n = 131`
are `0, 1, 130, 131`" is a classification of invariant **subspaces**.  GGMP's
collapse asks only for `π(F) = F`.  Does the saving survive without the
subspace?

**Answer, derived.**  Yes.  `131` is prime, so every `x ∉ F_2` has an orbit of
size exactly `131`, and a union of `k` orbits is `π`-stable at any size.  The
collapse is therefore available at ECC2K-130, the product law becomes
`m·2^n/n`, and the bottom of the family moves from `2^132.58` to `2^124.99` —
`2^64.18×` rho instead of `2^71.77×`.  What it costs is the subspace root-find:
an orbit union has no low-degree membership polynomial, so the base must be
**materialised** (`2^44.53` points at `m = 3`, about 300 TB; `2^23.42`, about
90 MB, at `m = 6`).  **`131` being prime forbids having both the collapse and
the implicit base, not the collapse itself.**

**The concrete shape to build it out of.**  "A union of orbits" is not a
construction.  The natural one is a **Hamming-weight shell in a normal basis**:
in the basis `{β, β², β⁴, …}` the Frobenius is a cyclic shift of the 131
coordinates, so the weight-`w` shell is shift-invariant by definition, and
because 131 is prime every shell with `0 < w < 131` is exactly `C(131,w)/131`
orbits of size 131 — no fixed points to special-case.  The sizes land on what
§6 of the background note asks for almost exactly:

| `m` | orbit base wants | weight shell | shell size |
|---:|---:|---:|---:|
| 3 | `2^44.53` | `w = 9` | `2^44.43` |
| 4 | `2^33.90` | `w = 6` | `2^32.54` |
| 5 | `2^27.58` | `w = 5` | `2^28.15` |
| 6 | `2^23.42` | `w = 4` | `2^23.48` |

Shells jump three to five bits apart, so finer sizes need a union of adjacent
shells or a shift-invariant subset of one; the membership test is a popcount
either way.  **The representation is already in the repository and the factor
base is not**: `ecc2k130/FROBENIUS-NETWORK.md` applies Frobenius as a fixed
bit-permutation network over the permuted type-II ONB, and
`src/bin/ic/params.rs` carries ECC2K-130 only as an abstract normal-basis
profile with "point coordinates not imported".  So E6 needs the base built, not
merely selected.

*(A correction to something adjacent: `research/nagao_relations/README.md` on
main states that this repository "already uses normal-basis Hamming-weight
sets" as factor bases.  It does not — what exists is the ONB as a **field
representation** for the rho client.  The observation underneath that sentence
is right and is the same one as §6.1 of the background note; only the claim that
it is implemented is not.)*

**What is still to run.**  The derivation assumes `|F|/n` relations suffice,
which needs them to be **independent over the orbit unknowns**.  That is the
part most likely to be wrong, and it is measurable on E1's ladder:

- build an orbit-union base of `k` orbits at `n = 19, 29, 37`;
- collect relations, rewrite them onto orbit representatives with their `λ^{k_i}`
  weights, and measure the **rank** of the resulting matrix against `|F|/n`;
- report the realised saving as `rank-limited relations needed / (|F|/n)`.

**Falsifier.**  A realised saving differing from `n` by more than 20%, which
would mean orbit relations carry a dependency the unknown-count argument misses
— exactly the failure mode GGMP's own "multiple invariant factor bases" remark
worries about (`RESEARCH_KOBLITZ_INDEX_CALCULUS.md`, open problems: "a single
one need not yield `n` *independent* relations").

---

## Ordering, and what each one would change

| | experiment | cost | if it confirms | if it falsifies |
|---|---|---|---|---|
| 1 | **E1** ladder | CPU-days | the `n = 131` derivation is evidence, not arithmetic | the whole note's exponent is wrong |
| 2 | **E6** rank | hours on E1 | `m·2^n/n` stands as the family's floor | the seven-bit correction reverses |
| 3 | **E3** target-agnosticism | days | proposed detectors are scored against one floor, cost *and* memory, with no credit for answering yes/no only | a cheap detector exists — it reaches rho from `m = 4` whether or not it hands over the witness |
| 4 | **E2** flatness | hours | `l` is confirmed not to be a lever | the cancellation is inexact and the exponent moves |
| 5 | **E4** large primes | days | one more classical lever priced and closed | the first sub-rho decomposition cell in this repository |
| 6 | **E5** dispersion | hours | the Poisson model is safe to extrapolate | a correction term to the saturation threshold |

E3 is the one whose design moved most while being written, and in the end it
moved *against* the curve: its headline went from "a free detector never reaches
rho", to "it reaches rho at three of seven summand counts with a table the size
of the algorithm", to "the distinction it was built on does not exist for any
detector anyone would propose". Two of those three were mine and both were
accounting, not measurement — which is the argument for pre-registering a
boundary and then checking it against a run rather than against intuition.  E1 is the one to run first,
because every other row on this page is scored against a law that only E1 can
confirm actually governs these curves.

## What is deliberately not on this list

- **Anything that changes the curve or the field.**  Isogeny class search (L5),
  extension fields, and quasi-subfield polynomials are each measured and closed
  in their own notes; re-opening them needs a new mechanism, not a new run.
- **A better Gröbner engine.**  The deciding degree grows as `l/2`
  (`RESEARCH_SEMAEV_DECOMPOSITION.md`), and F5's criteria remove zero reductions
  without lowering the degree of regularity, so a better implementation moves the
  constant and not the slope.
- **Factor bases that are neither subspaces nor Frobenius-stable.**  A base
  needs either a cheap algebraic membership test or materialisation; E6 covers
  the second, and the first is exactly the quasi-subfield census that came back
  empty at 131.  Note this excludes less than it first appears: normal-basis
  weight shells, which look like a third shape, are shift-invariant and so are
  orbit unions, already inside E6.  A base outside both would be a new idea, and
  this list is for experiments, not for ideas.
