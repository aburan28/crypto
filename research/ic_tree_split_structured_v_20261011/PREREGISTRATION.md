# Is the `m = 4` tree-split set cut out at low degree on a *structured* factor base?

Registered before any registered cell runs. §5 lists everything that was run before
registration and what was seen. The decision rule is registered twice: in prose (§4) and as
a Lean function, [`lean/TreeSplit/Rule.lean`](lean/TreeSplit/Rule.lean), whose hash is in
§7 and in the cairn registration record (§8).

## 1. Why

- **The predecessor.** [`ic_tree_split_cutout_20260930`](../ic_tree_split_cutout_20260930/RESULTS.md)
  asked whether the tree-split set `T2 = {x(P + Q) : P, Q ∈ F}` is cut out at lower degree
  than a random set of its size. Verdict **CLOSED**: gap 0 in all 48 cells. Every cell used
  a **random** GF(2)-subspace `V` for the factor base.
- **What it left open.** Its own scope line: "one random subspace per seed". The
  first-fall-degree literature, and this repository's descent-expansion predictor
  (`src/cryptanalysis/descent_expansion.rs`, prediction P2), say low degree is plausible
  exactly where the descent has structure: a low-expansion basis or a Frobenius-stable
  subspace, not a random one. A random `V` is the case most likely to look random.
- **The question.** Does a structured `V` give a `T2` with a low-degree description that a
  random `V` does not?
- **The prior is closure.** `T2` is the image of `2ℓ` bits under a map of high algebraic
  degree in the `n` field bits. Structure in `V` is linear, and the group law is not.

## 2. Objects, arms and controls

The curve is `K_a : y² + xy = x³ + a·x² + 1` over `GF(2^n)`, `n` odd. The field polynomial
is the first irreducible trinomial, else the first pentanomial in lexicographic order. That
is the predecessor's `find_irreducible`, reimplemented. `F` is the set of points with
`x ∈ V \ {0}`, and `T2` as above, with `P + Q ≠ O`.

| arm | `V` | role |
|:--|:--|:--|
| `poly` | `span{1, z, …, z^{ℓ−1}}` in the polynomial basis, with `2ℓ − 1 ≤ n`, so a product of two elements of `V` never reduces mod the field polynomial | **primary** |
| `frob` | `ker g(σ)`, where `σ` is squaring, `g` divides `t^n − 1` with `deg g = ℓ`, and `g ≠ t^d − 1` for `d \| n` (so not a subfield). `V` is Frobenius-stable, so `F` and `T2` are unions of Frobenius orbits. `g` is the `(seed mod #candidates)`-th admissible product of irreducible factors, in factor-subset bitmask order | **primary** |
| `rand` | a uniformly random `ℓ`-dimensional subspace: the predecessor's object at a new `n` | replication |
| `subfield` | `GF(2^ℓ) = ker(σ^ℓ − 1)`, `ℓ \| n`. On a Koblitz curve `F ⊂ E(GF(2^ℓ))` (because `n/ℓ` is odd), so `T2 ⊆ GF(2^ℓ)`, a linear subspace | positive control (validity gate) |

The **control `R`** for each cell is a uniformly random subset of `GF(2)^n \ {0}` of size
`|T2|`.

Every random choice comes from a `splitmix64` stream keyed by `(seed, n, a, ℓ, tag)`, with
tag 1 for `V`, 2 for `R`, 3 for the curve arm's samples and 4 for `R`'s. So a cell is a
pure function of its arguments.

## 3. Metric

This is the predecessor's metric, unchanged. For a set `S ⊂ GF(2)^N` and a degree `D`, `h_S(D)` is
the GF(2) rank of the evaluation matrix of all multilinear monomials of degree `≤ D`
(`M(D)` columns) at the points of `S`.

- `deficit(D) = min(|S|, M(D)) − h_S(D)`.
- The **sampled cut-out degree `D_s`** is the smallest `D ≤ 6` at which none of `K = 2000`
  points drawn uniformly from `GF(2)^N \ S` satisfies every degree-`≤ D` polynomial
  vanishing on `S`. A point satisfies them all exactly when its monomial vector lies in the
  row space.
- **Censoring.** A degree with `M(D) > 90,000` is skipped. An arm that stops before reaching
  `D_s`, or reaches `D = 6` with false zeros left, is censored. A cell that exits non-zero
  or exceeds its 1,800 s limit is censored too. Censored cells are excluded from the verdict
  and listed. Censoring is never evidence either way.
- **What `D_s` resolves.** It is a screen, not an exact cut-out degree. With 2,000 samples it
  detects a degree-`D` closure larger than `S` by at least 0.15% of the space with
  probability 95%, and can miss a smaller excess.

## 4. Cells and decision rule

[`run.sh`](run.sh) runs every cell once, in the order below. The harness is
[`examples/tree_split_structured.rs`](../../examples/tree_split_structured.rs).

| arm | `n` | `a` | seeds | `ℓ` | cells |
|:--|:--|:--|:--|:--|--:|
| `subfield` | 21 | 0, 1 | 20261011 | 7 | 2 |
| `poly` | 17, 19, 21 | 0, 1 | 20261011 | 4–8 | 30 |
| `frob` | 21 | 0, 1 | 20261011, 20261012 | 4–8 | 20 |
| `frob`, supplementary | 17 | 0, 1 | 20261011 | 8 | 2 |
| `rand` | 21 | 0, 1 | 20261011 | 4–8 | 10 |

- **Why these `n` for `frob`.** At `n = 21`, `t^21 − 1` factors with degrees
  `1, 2, 3, 3, 6, 6`. So a Frobenius-stable non-subfield `V` exists at every `ℓ` from 4 to 8,
  with 2, 2, 5, 2 and 3 candidates. At `n = 17` (degrees `1, 8, 8`) one exists only at
  `ℓ = 8`, and at `n = 19` (`1, 18`) none does. The `n = 17` cells give no slope, so they
  are reported but do not enter the verdict.
- **Series and gap.** A series is one `(arm, n, a, seed)`. The arms have six `poly`, four
  `frob` and two `rand` series. For each uncensored cell, the gap is
  `D_s(curve) − D_s(R)`.

**Verdict, applied separately to `poly`, `frob` (n = 21) and `rand`:**

- **ALIVE** if the gap at the largest uncensored `ℓ` is `≤ −2` in every series, **and** at
  least half of the series (rounded up) have an OLS slope of `D_s(curve)` on `ℓ` of `≤ 0.15`.
  That is, the lower median of the slopes is `≤ 0.15`.
- **CLOSED** if every uncensored gap is `≥ −1` in every series. A curve set that sits above
  its control counts as closed: the lever needs a lower degree.
- **INCONCLUSIVE** otherwise.

This is the predecessor's rule. One clause is restated so that it is computable in
integers: its "median slope" becomes the lower median, `#{slope ≤ 0.15} ≥ ⌈k/2⌉`. With an
even number of series this is weaker than the midpoint median, so it can only make ALIVE
easier. That direction is disclosed, not hidden.

**Validity gate.** The run is valid only if every `subfield` cell has `T2 ⊆ V` and a
degree-1 deficit of `≥ n − ℓ` (`= 14`) in the curve arm. The gate uses the deficit, not
`D_s`. A `T2` that fills most of a small subfield leaves non-members of `T2` inside the
subfield that satisfy every degree-`≤ D` polynomial vanishing on `T2`. That is a correct
reading of `D_s`, not a failure of the screen (§5 shows it once). The harness self-test
(§5) must also pass. A failed gate makes the run **invalid**, whatever the verdicts say.

**What each verdict would mean.**

- **`poly` or `frob` ALIVE.** That family's `T2` has a low-degree description at
  `n ≤ 21`, `ℓ ≤ 8`. It motivates a pre-registered solving-degree audit of `Z ∩ A_R` on that
  family, against the `c* = 0.25` gate of the `m = 4` audits. It is **not** an exponent
  claim.
- **CLOSED.** Structure in `V` of that kind does not give `T2` a low-degree shortcut at
  these sizes. The tree-split reopening condition stays unmet for that family. It says
  nothing about `Z`, about other structured `V` (normal-basis windows, other `n`), or about
  larger `n` and `ℓ`.
- **`rand` not CLOSED.** This disagrees with the predecessor at `n = 21`. It is reported as
  a replication failure, to be explained before any structured-arm verdict is used.

## 5. What ran before registration

- **Harness self-test** (`tree_split_structured selftest`). It checks:
  - the field inverse;
  - commutativity, associativity and closure of the group law on 30 random points, for
    `(a, n) = (0, 13), (1, 17), (0, 7), (1, 21)`;
  - `#K_0(GF(2^7)) = 116`, with `116·P = O` on five points;
  - the factorisations of `t^21 − 1` (degrees `1, 2, 3, 3, 6, 6`) and of `t^17 − 1`
    (`1, 8, 8`);
  - `GF(2^7) = ker(σ^7 − 1)` in `GF(2^21)`, and that every `frob` `V` at `n = 21` has
    dimension `ℓ` and is closed under squaring;
  - the predecessor's three linear-algebra checks: a dimension-8 affine subspace of
    `GF(2)^20` gives `D_s = 1`; the quadric `x₀x₁ + x₂ = 0` in `GF(2)^14` gives `D_s = 2`
    with deficit 3 at `D = 2`; 1,000 random points of `GF(2)^20` give `D_s = 3`.

  All passed. The self-test printed the `frob` candidate counts quoted in §4.
- **Timing (not a registered cell).** `rand`, `n = 21`, `a = 1`, `ℓ = 8`, seed 1, both arms:
  15 s wall on 4 vCPUs. The output was discarded unseen.
- **Lean replay development, `n = 13` and `15` (not registered sizes), seed 7.** This was
  used to test `Replay.lean` against the harness, and was seen in full:
  - `poly`, `ℓ = 4, 5`, `a = 0, 1`: gap 0 in all eight cells. `|T2|` = 49, 194, 62, 280
    (`n = 13`) and 25, 64, 100, 515 (`n = 15`).
  - `subfield`, `n = 15`, `a = 0`, `ℓ = 5`: `T2 ⊆ V`, degree-1 deficit 10 (`= n − ℓ`),
    `|T2| = 22`. But `D_s(curve) = 3` against `D_s(R) = 2`. This is the effect described
    in the gate paragraph of §4, and it is why the gate reads the deficit.
  - The Lean replay reproduced all eight `poly` cells exactly. Changing one recorded `D_s`
    made it fail.
- **Seen but not a result.** None of the registered `(n, a, seed)` cells has been run, for
  any arm.

## 6. Implementation and scope

- **Code.** `examples/tree_split_structured.rs` is native Rust (AGENTS.md "no Python"). It
  builds `V`, `F`, `T2` and `R`, runs the screen with its own GF(2) elimination, and has
  three more subcommands:
  - `readout` prints the table and the verdicts;
  - `lean` writes `lean/TreeSplit/Data.lean` from the records;
  - `selftest`.
- **Not timed.** The metric is a degree, so the AGENTS.md §10 isolation rule for timed runs
  does not apply. The per-degree `secs` fields are logged and are not measurements.
- **Scope.** Two Koblitz curves, `n ∈ {17, 19, 21}`, `ℓ ≤ 8`, the polynomial-basis window,
  Frobenius-stable kernels at `n = 21`, and `T2` only (not `Z`). Anything beyond that is
  untested.

## 7. Lean checks (registered)

Lean has no clock: none of this is about time. It checks that the reported result is the
one the registered rule gives on the recorded data, and that part of the recorded data is
reproducible by a second implementation.

| file | what it proves | how |
|:--|:--|:--|
| `Rule.lean` | the decision rule and validity gate of §4 as total functions, with sanity examples on boundary tables | kernel `decide` |
| `Data.lean` | the recorded cells, generated from `runs/registered/cells.jsonl` by the harness | generated |
| `Verdict.lean` | `report cells = ⟨gate, poly, frob, rand⟩` for the observed values, so `RESULTS.md` cannot state a verdict the rule does not give | kernel `decide` |
| `Replay.lean` | a separate implementation recomputes `|T2|`, `D_s(curve)` and `D_s(R)` for the 12 `poly` cells with `ℓ ≤ 5` and matches the records exactly. Its elimination uses the highest-bit pivot; the harness uses the lowest | `native_decide`: trusts the compiler (axiom `Lean.ofReduceBool`) |

**If `Replay.lean` disagrees with the records, the run is invalid** until the cause is found.

Registered hashes:

- `lean/TreeSplit/Rule.lean`: sha256 `714601c110d88e23202206b741c92f9d198554410b1b087841f3b2e36cd0760d`
- `lean/TreeSplit/Replay.lean` and the harness: their hashes at the registration commit, as
  recorded by git.

## 8. cairn registration

The rule is posted to a cairn log (`cairn/`) before the run, as a `lean` objective. Its
preamble is `Rule.lean`. Its statement is the conjunction of the rule's boundary-case
behaviours, and its statement text carries this file's registration commit. A cairn
objective's id is the hash of its content, so the log fixes the rule's text at a time
before the run. After the run, a second `lean` objective states the observed `report`
over `Rule.lean` and `Data.lean`, closed by `decide`, and is verified in cairn's Lean jail.
Both ids are in `RESULTS.md`.
