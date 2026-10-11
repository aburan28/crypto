# Tree-split cut-out screen, round 2: subfield-free Frobenius kernels, normal-basis windows, `Z`, larger `n`

Registered before any registered cell runs. §5 lists what ran before registration. The
decision rule is registered twice: in prose (§4) and as the Lean function
[`lean/TreeSplit2/Rule.lean`](lean/TreeSplit2/Rule.lean), whose hash is in §7 and in the
cairn registration record.

## 1. Why

[Round 1](../ic_tree_split_structured_v_20261011/RESULTS.md) found the tree-split set `T2`
CLOSED on three families: the polynomial-basis window (`poly`), Frobenius-stable kernels
(`frob`), and a random-`V` replication. Its scope and post-hoc notes left four follow-ups,
and this round runs all four:

1. **`frob2`.** Round 1's `frob` excluded only the subfields themselves. Six of its cells
   had `V` inside `GF(2^7)`, and those were exactly its six gap −1 cells. `frob2` excludes
   every `g` dividing `t^d − 1` for a proper divisor `d` of `n`, so no `V` lies in a proper
   subfield.
2. **`nbw`.** A normal-basis window: `V = span{θ, θ², …, θ^{2^{ℓ−1}}}`, for a normal element
   `θ`. These are the low-expansion coordinates of the descent-expansion predictor
   (`src/cryptanalysis/descent_expansion.rs`). The window is not Frobenius-stable, because
   `ℓ < n`.
3. **`Z`.** The object `Z = {(e + f, e·f) : e, f ∈ T2} ⊂ GF(2)^{2n}` on structured `V`.
   For `Z`, the final tree-split equation is affine. The 2026-09-30 study ran `Z` on random
   `V` only.
4. **Larger `n`.** `poly` at `n = 23, 25`.

The prior is closure, from round 1 and from the 2026-09-30 study.

## 2. Objects, arms, controls

- **Curve, field polynomial, `F`, `T2`, controls and screen.** All exactly as in round 1
  (its PREREGISTRATION §2 and §3). Round 1's code is copied unchanged into
  [`examples/tree_split_structured_v2.rs`](../../examples/tree_split_structured_v2.rs), and
  round 1's harness stays frozen.

| arm | `V` |
|:--|:--|
| `frob2` | `ker g(σ)`, `g \| t^n − 1`, `deg g = ℓ`, and `t^d − 1 mod g ≠ 0` for every `d \| n` with `d < n`. `g` is the `(seed mod #candidates)`-th admissible factor-subset product, in bitmask order. At `n = 21` the candidate counts are 0, 2, 4, 2 and 3 for `ℓ = 4…8`, as asserted by the self-test. So `ℓ` starts at 5. At `n = 35` the deg-7 and deg-8 candidates are cubic·`Φ₅` and `(t + 1)`·cubic·`Φ₅` |
| `nbw` | `span{θ^{2^i} : i < ℓ}`, where `θ` is the first value of the `TAG_V` stream whose `n` conjugates have GF(2) rank `n` |
| `poly`, `subfield` | as in round 1 |

- **Controls.** For `T2`, `R` is a random subset of `GF(2)^n \ {0}` of size `|T2|`. For
  `Z`, `R` is a random subset of `GF(2)^{2n} \ {0}` of size `|Z|`. `ZR` is `Sym²` of a
  random set of `|T2|` nonzero field elements (stream tag 5; its samples use tag 6).
- **Gap.** `D_s(curve) − min` over the cell's controls.

## 3. Metric

Round 1's metric, unchanged:

- the sampled cut-out degree `D_s`, with `K = 2,000` samples, `D ≤ 6`, and at most
  90,000 monomials;
- censoring as in round 1, except that the per-cell limit here is 3,600 s.

## 4. Cells and decision rule

[`run.sh`](run.sh) runs every cell once, in order. That is 102 cells.

| group | object | arm | `n` | `a` | seeds | `ℓ` | cells | role |
|:--|:--|:--|:--|:--|:--|:--|--:|:--|
| G0 | T2 | `subfield` | 21 | 0, 1 | 20261011 | 7 | 2 | validity gate |
| G1 | T2 | `frob2` | 21 | 0, 1 | 20261011, 20261012 | 5–8 | 16 | **primary** |
| — | T2 | `frob2` | 35 | 0, 1 | 20261011 | 7, 8 | 4 | supplementary (no slope) |
| G2 | T2 | `nbw` | 17, 19, 21, 23 | 0, 1 | 20261011 | 4–8 | 40 | **primary** |
| G3 | T2 | `poly` | 23, 25 | 0, 1 | 20261011 | 4–8 | 20 | **primary** |
| G4 | Z | `poly` | 13, 15 | 0, 1 | 20261011 | 4, 5 | 8 | secondary |
| G5 | Z | `nbw` | 13, 15 | 0, 1 | 20261011 | 4, 5 | 8 | secondary |
| G6 | Z | `frob2` | 15 | 0, 1 | 20261011 | 4, 5 | 4 | secondary |

**Verdict per group (G1 to G6).** This is round 1's rule.

- **ALIVE** if the gap at the largest uncensored `ℓ` is `≤ −2` in every series, and at
  least half of the series (rounded up) have an OLS slope of `D_s(curve)` on `ℓ` of `≤ 0.15`.
  The `Z` groups G4–G6 drop the slope clause, as the 2026-09-30 study did for `Z`, because
  two `ℓ` values do not make a slope.
- **CLOSED** if every uncensored gap is `≥ −1`.
- **INCONCLUSIVE** otherwise.
- A series is `(group, n, a, seed)`, and censored cells are excluded.

**Validity gate.**

- G0 must give `T2 ⊆ V` with a degree-1 deficit of `≥ n − ℓ` in both cells.
- The harness self-test must pass.
- If either fails, the run is invalid.

**What each verdict would mean.**

- **ALIVE in G1, G2 or G3.** That family's `T2` has a low-degree description at these sizes.
  It motivates a solving-degree audit of `Z ∩ A_R` on that family. It is not an exponent
  claim.
- **ALIVE in G4–G6** is a `Z`-only signal. On its own it neither opens nor closes the route.
- **CLOSED.** That structure gives no low-degree shortcut at these sizes.

## 5. What ran before registration

- **Self-test** (`tree_split_structured_v2 selftest`). It checks everything round 1's
  self-test checks, plus:
  - the `frob2` candidate counts at `n = 21`;
  - that every `frob2` `V` is Frobenius-stable and lies outside `GF(2^3)` and `GF(2^7)`;
  - that every `nbw` `θ` is normal at `n = 17, 19, 21, 23`, with `dim V = ℓ`;
  - that `Sym²` is injective on unordered pairs.

  It passed.
- **Timing, seed 1 (not a registered seed). Outputs discarded unseen; exit status and
  wall time only:**

| cell | wall |
|:--|--:|
| `Z poly n = 15 a = 1 ℓ = 5` | 408 s |
| `Z frob2 n = 15 a = 1 ℓ = 5` | 25 s |
| `T2 frob2 n = 35 a = 0 ℓ = 8` | 30 s |
| `T2 poly n = 25 a = 0 ℓ = 8` | 38 s |

  `poly` `V` does not depend on the seed. So the first timing cell's curve set is the
  registered cell's curve set; its output was not read.
- **Nothing else was seen.**

## 6. Scope

- **Covered.** Koblitz curves, `n ≤ 25` (35 supplementary), `ℓ ≤ 8` for `T2`, and
  `ℓ ≤ 5` for `Z` at `n ≤ 15`.
- **The metric** is a sampled screen.
- **Not timed.** The `secs` fields are logged, not measurements.

## 7. Lean checks and cairn (registered)

- **`Rule.lean`** is the rule of §4. Its sha256 is
  `ec513363df351f6559a6d5df31bf107a993fd425c61c46ef986b0f8918bfb1e3`. Its boundary examples
  are checked by kernel `decide`.
- **`Data.lean`** is generated from `runs/registered/cells.jsonl` by the harness.
- **`Verdict.lean`** proves `report cells = …` for the observed verdicts by kernel `decide`.
  A CI workflow regenerates `Data.lean` and the readout from the records and runs
  `lake build`.
- **No independent replay this round.** Round 1's `Replay.lean` already cross-checked the
  shared screen and curve code on 12 cells.
- **cairn.** The rule is posted to `cairn/log.jsonl` as a `lean` objective before the run.
  After the run, the observed report is posted over `Rule.lean` and `Data.lean`.
