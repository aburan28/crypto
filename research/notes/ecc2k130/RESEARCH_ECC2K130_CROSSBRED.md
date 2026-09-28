# Crossbred on the Koblitz decomposition systems

**Experiments:** `examples/crossbred_bench.rs` (oracle ladder),
`examples/koblitz_crossbred_panel.rs` (whole pipeline)
**Frozen evidence:** `research/crossbred_20260920/oracle_ladder.txt`,
`research/crossbred_20260920/e2e_panel.txt`, `.../e2e_panel.json`
**Under test:** `DecompositionStrategy::Crossbred`, new in this round
**Background:** [`RESEARCH_ECC2K130_ROUTES.md`](RESEARCH_ECC2K130_ROUTES.md)
(Route 1, which declared this round's metric and falsifier),
[`RESEARCH_ECC2K130_DECOMPOSITION.md`](RESEARCH_ECC2K130_DECOMPOSITION.md)
(the `m·2^131` oracle bound this does not touch),
[`RESEARCH_ECC2K130_RR_SOLVER_PANEL.md`](RESEARCH_ECC2K130_RR_SOLVER_PANEL.md)
(the last solver that won its phase and lost the pipeline)

**The question.** Route 1 ranked Crossbred first because the survey found no
publication combining it with binary ECDLP decomposition systems, and because
`src/cryptanalysis/crossbred.rs` was already built and correctness-gated
against matrix-F4 but had never reached relation collection. Two things were
asked for: write up what the existing bench finds, and make Crossbred a
`DecompositionStrategy` so the scaling ladder runs end to end.

**Bottom line. Crossbred is the fastest decomposition oracle in the tree, its
declared falsifier does not fire, and it still loses the pipeline to a double
loop.** Three results, and conflating them is how this route gets over-sold:

- **The oracle (§2).** `xb/F4 = 0.023` in bit operations at `m = 3`, and
  **`0.696` in wall clock** on the largest rung both engines finish — the
  number Route 1 named. The falsifier was `xb/F4 ≥ 1`; it is not met, so the
  route does not close.
- **The trend inside that win (§3).** The wall-clock advantage *shrinks with
  size*: `0.021 → 0.053` at `m = 2` over `n = 9 → 13`, and `0.040 → 0.696` at
  `m = 3` over `n = 7 → 9`. Two rungs are not a fit, but the direction is the
  wrong one and the `m = 3` margin is nearly gone by `n = 9`.
- **The pipeline (§4).** End to end, exhaustive enumeration beats every
  algebraic oracle by `3.4×` to `1800×`. Crossbred is the best *algebraic*
  arm at `m = 2` (`2.1×` to `7.6×` over Gröbner) and the **worst** arm at
  `m = 3`, where it runs `2.7×` slower than Gröbner and `4.4×` slower than
  SAT. The oracle win does not survive contact with relation collection.

## 1. Boundary and unit

Two units, because the round has two phases and `AGENTS.md` §2 wants the
boundary constant in whichever is reported.

*Oracle phase.* **Bit operations**, one 64-bit word XOR = 64 bit operations,
which is the bench's one stated conversion. The boundary is exhaustive search
over the Boolean system: `2^v` points each evaluating `T` monomials, so
`2^v · T`. That is what an oracle must beat to be doing algebra rather than
search. Route 1's own metric is the *wall-clock* `xb/F4` ratio, reported
alongside.

*Pipeline phase.* **Relations per second** over `relation_collection_ns`, with
`DecompositionStrategy::Enumerate` as the reference arm — the null object, the
same role brute-force pair enumeration plays in
[`RESEARCH_ECC2K130_RR_SOLVER_PANEL.md`](RESEARCH_ECC2K130_RR_SOLVER_PANEL.md).

*What neither unit can move.* Relation collection with any enumeration oracle
costs `Θ(2^n)` independent of `ℓ`, measured flat within `1.5×` across
`ℓ = 6…10`. A faster oracle buys reach, not an exponent. Nothing below is a
statement about `2^131`.

## 2. The oracle ladder

`cargo run --release --example crossbred_bench`. Twelve targets per row, same
curve, factor base, arity and targets across engines.

| n | m | v | brute (bit ops) | F4 (bit ops) | crossbred (bit ops) | xb/brute | xb/F4 | xb wall | F4 wall | **xb/F4 wall** | agree |
|--:|--:|--:|----------------:|-------------:|--------------------:|---------:|------:|--------:|--------:|---------------:|:-----:|
| 9 | 2 | 12 | 902,144 | 27,529,653 | 108,117 | 0.120 | 0.004 | 0.1 ms | 3.8 ms | **0.021** | yes |
| 13 | 2 | 24 | 18,570,980,010 | 10,393,078,458 | 24,350,101 | 0.001 | 0.002 | 10.7 ms | 203.5 ms | **0.053** | yes |
| 7 | 3 | 16 | 20,000,395 | 25,174,603 | 320,436 | 0.016 | 0.013 | 0.2 ms | 6.0 ms | **0.040** | yes |
| 9 | 3 | 27 | — | 1,294,445,989 | 30,338,821 | — | 0.023 | 79.3 ms | 114.1 ms | **0.696** | yes |

`agree` is the correctness gate: crossbred's root set equals exhaustive
search's exactly where exhaustive search fits, and past `v = 24` weakens to
two one-sided checks (every crossbred root verifies against the original
equations; every root F4 found is one crossbred found). A row that does not
say `yes` is not a result.

**Route 1's falsifier does not fire.** It was `xb/F4 ≥ 1` across the ladder
with no `(D, k)` doing better. The worst row is `0.696`. The route stays open
on its own terms.

The `(D, k)` sweep is in the frozen output. The shape: at `D = 2` no crossbred
space exists below `k = 5`, and `kernel_dim` saturates at `9` by `k = 6`; at
`D = 4` the kernel grows monotonically with `k` (`12 → 403` over `k = 2…7`)
while preprocessing word ops *fall* (`925,106 → 356,150`), because a larger
enumerated set leaves fewer Macaulay columns to eliminate. No row produced a
filter, at any `(D, k)` — the specialised-degree-0 polynomials that make
Crossbred's sweep cheap on random MQ systems do not appear here at all. That
is the structural reason the margin decays in §3.

## 3. The advantage shrinks with size

| m | n | xb/F4 wall |
|--:|--:|-----------:|
| 2 | 9 | 0.021 |
| 2 | 13 | 0.053 |
| 3 | 7 | 0.040 |
| 3 | 9 | **0.696** |

Two points per `m` is not a fit and this note does not report an exponent.
What it reports is the sign: the ratio rises on both arms, and at `m = 3` it
rises `17×` for one rung of `n`. Extrapolating that shape — which would be an
extrapolation, marked as one — puts the crossing at `n = 10` or `11` at
`m = 3`, which is thirty years of nothing before `n = 131`.

The mechanism is §2's last observation. Crossbred's search phase is cheap
when filters carry most of the sweep; here there are none, so every one of the
`2^k` points reaches a linear solve, and the Macaulay extraction that produced
the crossbred space is paid *per target* because the system's constant terms
move with `x_R`. At `n = 9, m = 3` extraction plus search is 79.3 ms against
F4's 114.1 ms, and the gap closes from the crossbred side.

## 4. The pipeline, where it loses

`cargo run --release --example koblitz_crossbred_panel`. Whole index-calculus
run, `strategy` the only variable: same curve, factor base, `m`, seed, trial
cap and linear algebra, `allow_direct_relation: false` so the generic
`aG + bQ = O` shortcut cannot solve an instance without calling the oracle.

| n | m | strategy | relations | trials | relations/s | **vs. enumerate** | solved |
|--:|--:|:---------|----------:|-------:|------------:|------------------:|:------:|
| 7 | 2 | crossbred | 3 | 3 | 12,126.9 | 0.29 | yes |
| 7 | 2 | groebner | 3 | 3 | 5,821.6 | 0.14 | yes |
| 7 | 2 | sat | 3 | 3 | 10,278.0 | 0.25 | yes |
| 7 | 2 | **enumerate** | 3 | 3 | **41,488.6** | 1.00 | yes |
| 9 | 2 | crossbred | 3 | 3 | 4,672.2 | 0.09 | yes |
| 9 | 2 | groebner | 3 | 3 | 613.9 | 0.01 | yes |
| 9 | 2 | sat | 3 | 3 | 1,654.1 | 0.03 | yes |
| 9 | 2 | **enumerate** | 3 | 3 | **51,103.0** | 1.00 | yes |
| 9 | 3 | crossbred | 5 | 5 | 11.8 | 0.0006 | yes |
| 9 | 3 | groebner | 5 | 5 | 32.3 | 0.0015 | yes |
| 9 | 3 | sat | 5 | 5 | 51.9 | 0.0024 | yes |
| 9 | 3 | **enumerate** | 5 | 5 | **21,330.5** | 1.00 | yes |

Every row recovers the planted scalar. The full panel including the second
`n = 9, m = 2` instance is in the frozen output.

Two readings, both of them the point:

- **Among algebraic oracles, Crossbred wins at `m = 2` and loses at `m = 3`.**
  `7.6×` over Gröbner at `n = 9, m = 2`; `0.37×` at `n = 9, m = 3`. The
  crossover between §2's oracle win and this pipeline loss is the per-target
  extraction: the oracle table charges it once per target and so does the
  pipeline, but the pipeline also runs the arms that do not pay it.
- **No algebraic oracle beats the double loop at any rung measured.** The
  worst case is `1800×` against enumeration at `n = 9, m = 3`. This is the
  same shape as the Riemann–Roch panel's result on the real curve, reached by
  a different solver family, and it is the reason neither is an attack.

### 4.1 One instance refuses every arm

`K_0 / F_2^7` at `m = 2` collects **zero** relations in 20,000 trials under
all four strategies. That is a property of the instance, not of any oracle:
four independent engines agreeing on "no relation" is what an admissibility
failure looks like, and `RESEARCH_ECC2K130_DECOMPOSITION.md` §3 documents the
cofactor class that makes some `m` decompose nothing on small `K_1`. It is
kept in the panel as a matched control rather than dropped, because a row
where every arm fails is evidence about the instance and dropping it would be
choosing favourable seeds.

## 5. Classification

By `AGENTS.md` §3, for each result:

| result | what moved | class |
|:--|:--|:--|
| `xb/F4 = 0.023` bit ops, `0.696` wall, at `m = 3` | oracle-phase constant fell; ratio to the `Θ(2^n)` collection floor flat | **engineering** |
| Crossbred fastest algebraic oracle at `m = 2` end to end | pipeline constant fell against Gröbner and SAT, not against the reference arm | **engineering** |
| Reading §2's oracle win as a pipeline win | oracle cost fell `43×` against F4 while the pipeline got *slower* than Gröbner at `m = 3` | **relabelling** — named here so it is not claimed later |
| Anything about `2^131` | nothing | not measured; §1 |

No row is an **advance**. The floor is `Θ(2^n)` relation collection and every
arm here sits against it identically.

## 6. What would close the route

The falsifier did not fire, so Route 1 stays open, but its `m = 3` margin is
`0.696` and shrinking. It closes if either:

- the `xb/F4` wall-clock ratio at `m = 3` exceeds `1` at `n = 11` or `n = 15`
  — one more rung decides it, and the bench already takes `n` as an argument;
  or
- a `(D, k)` choice producing **filters** on these systems cannot be found.
  Zero filters across the entire sweep is the structural finding of §2, and
  filters are the whole reason Crossbred's sweep is cheap. Without them the
  method is a Macaulay preprocessing step in front of `2^k` linear solves,
  which is what the decaying margin is measuring.

The GPU argument in Route 1 — the sweep is `2^k` independent points with no
shared state — survives all of this, and is also worth nothing until the
margin stops shrinking: porting a `0.696` that is heading for `1.0` to a GPU
buys a constant on a phase that is not the bottleneck.
