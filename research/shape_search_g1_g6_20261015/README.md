# Shape search at `m ≥ 3`: G1, G2, G3, G5 — results

**Stage diagnostic.** These are refutation degrees of descended Semaev
systems on toy binary fields (`n ≤ 13`).

- `S`, end-to-end cost and speedup are **unset**. Nothing here moves a ratio
  to rho.
- Lane `shape-search`, 2026-10-15.
- Protocol: [`PROTOCOL.md`](PROTOCOL.md), committed in `3bdf9b06e` before
  any measured cell.
- Instrument: `examples/factor_base_shape_ladder.rs`. The measuring code is
  that commit's. The `score` subcommand was added in `048403766`, with no
  change to measurement.
- Goals: frozen in commit `110117d8c` of
  `research/notes/index-calculus/RESEARCH_FACTOR_BASE_SHAPE_SEARCH.md` §7.

Labels: **measured** values come from the JSONL records in `runs/`.
**Derived** values are computed from them by `factor_base_shape_ladder
score`, recorded in [`summary.md`](summary.md), or worked by hand where the
text says so.

## Verdicts

| goal | frozen kill condition | verdict | scope |
|:--|:--|:--|:--|
| G1 | no cell with the gap | **killed on the scaled cells**: `gp`, `fp` and `coord` sit exactly at random in every cell (gap 0.000). **Not testable on the frozen cells** `(13,5)`, `(17,6)`, `(19,7)` (resource) | prime `n ∈ {5, 7, 11}`, `l = ⌈n/3⌉`, `m = 3`, class-0 targets |
| G2 | rise ≥ 1 against the plain chain | **survives at these sizes**: rise 0.000 in every arm (`s ∈ {1, 2, l}`, `second` and `ordered`, random and `fp`) | `(5,2)`, `(7,3)`, `(9,3)` |
| G3 | class-0 `D` flat in the panel to 0.3 | **killed**: panel spread 0.000 at every cell | `(5,2)`, 6 moduli; `(7,3)`, 7 moduli; `(11,4)`, 5 moduli |
| G4 | — | **not built** (budget) | — |
| G5 | class-0 `D*` flat in `β` to 0.3 | **killed** (the expected null): `β`-spread 0.06–0.19 | `m = 2`, `(7,3)`, `(9,4)`, `(11,5)`, `(13,6)` |
| G6 | — | **not run**: no survivor of G1 or G3. G2 "survives" a kill condition about degree, not about shape. Its question, whether the ordered tuple's yield gain carries to `n = 31`, has no G6 instrument here | — |

The positive control works: the subfield `F_8 ⊂ F_{2^9}` refutes at 5,
where random reads 6. The instrument can see a shape effect of one degree
at `m = 3`. The negative control (`random2` vs `random`) has gap 0.000 in
every cell.

## What ran (measured)

The host was a four-core container shared with four other lanes. Every run
used `nice -n 10` and one thread, with the address space capped at 2.5 GB.
Seed `20261015`. Sampled peak RSS and wall time are in `runs/*.meta`. They
are practicality notes, not metrics.

### G1: shape at `m = 3`, class-0 targets, mean refutation degree

| cell | unknowns | random | random2 | gp | fp | coord | subfield |
|:--|--:|--:|--:|--:|--:|--:|--:|
| `(5, 2)`, `d_max` 7 | 11 | 5.000 (8) | 5.000 (8) | 5.000 (8) | 5.000 (8) | 5.000 (7) | — |
| `(7, 3)`, `d_max` 7 | 16 | 6.000 (8) | 6.000 (8) | 6.000 (8) | 6.000 (8) | 6.000 (8) | — |
| `(9, 3)`, `d_max` 7 (control) | 18 | 6.000 (8) | — | — | — | — | **5.000 (8)** |
| `(11, 4)`, `d_max` 6 | 23 | 6.000 (6) | 6.000 (6) | 6.000 (6) | 6.000 (6) | 6.000 (6) | — |

- The number in brackets is the count of class-0 unsatisfiable targets
  measured.
- Every draw resolved by refutation. No draw was censored or hit the caps.
- **Within each cell every draw reads the same degree.** That holds for
  every shape and every instance, so the instance-mean spread is 0.
- `coord` at `(5, 2)` has 7 draws, not 8. Its window is fixed, and only 7
  class-0 targets of `F_32` are unsatisfiable for it, under sampling without
  replacement. That was the whole field.

**The frozen cells.**

- `(17, 6)` and `(19, 7)` are out of budget (derived, PROTOCOL §3). They
  have 35 and 40 unknowns.
- `(13, 5)` is not affordable either (measured). The feasibility probe, one
  random draw at `d_max = 6`, did not finish in 60 min under the cap on the
  shared host, and it was stopped. The protocol made the one-sided
  `(13, 5)` arm conditional on that probe, so `(13, 5)` was not run.
- The repository's earlier `(13, 5)` draws read `≥ 7` for random
  (`RESEARCH_DREG_MEASUREMENT.md` Result 7). They are not class-conditioned.

### G3: coordinate window over the modulus panel

| cell | moduli | class-0 mean `D` per modulus | spread |
|:--|--:|:--|--:|
| `(5, 2)` | 6, all the irreducibles of degree 5 | 5.000 each (6 draws each) | 0.000 |
| `(7, 3)` | 7: 4 trinomials, `z^7+z^3+z^2+z+1`, and the two heaviest | 6.000 each (6 draws each) | 0.000 |
| `(11, 4)` | 5: `z^11+z^2+1`, `z^11+z^9+1`, `z^11+z^4+z^2+z+1`, and the two heaviest | 6.000 each (4 draws each) | 0.000 |

### G2: τ-twisted chains (`d_max = 7`, 3 instances × 2 class-0 unsat targets)

**Degree (measured).** Every arm in every cell reads the plain chain's
degree exactly:

| cell | degree |
|:--|--:|
| `(5, 2)` | 5 |
| `(7, 3)` | 6 |
| `(9, 3)` | 6 |

The arms are `plain`, `second_s` and `ordered_s` for `s ∈ {1, 2, l}`, with
random and `fp` bases.

- The rise is 0.000 everywhere.
- At `(5, 2)`, `s = 2` and `s = l` coincide, so those arms carry 12 draws.

**Yield (derived).** This is the class-0 satisfiable fraction over the
targets drawn before the measured ones, so it is a coarse estimate from
2–127 drawn targets per instance.

| cell | shape | plain | `second_s1` | `ordered_s1` | `ordered_s2` | `ordered_sl` |
|:--|:--|--:|--:|--:|--:|--:|
| `(5, 2)` | random | 0.333 | 0.700 | 0.826 | 0.818 | (= s2) |
| `(7, 3)` | random | 0.571 | 0.941 | 0.952 | 0.951 | 0.945 |
| `(9, 3)` | random | 0.250 | 0.333 | 0.455 | 0.667 | 0.500 |
| `(5, 2)` | fp | 0.455 | 0.455 | 0.455 | 0.455 | (= s2) |
| `(7, 3)` | fp | 0.538 | 0.538 | 0.538 | 0.538 | 0.538 |
| `(9, 3)` | fp | 0.333 | 0.400 | 0.455 | 0.538 | 0.538 |

**What the yield table shows:**

- **Random bases.** The twisted arms decompose more class-0 targets than the
  plain chain at every cell. That is consistent with the hypothesis's
  "ordered tuple, yield × m!", but these small fields saturate, so the factor
  6 cannot be read off.
- **`fp` bases at `(5, 2)` and `(7, 3)`.** Every arm saw the identical
  satisfiable/unsatisfiable sequence on the same target stream. For example,
  instance 0 at `(7, 3)` measured targets 21 and 65 in all seven arms. This
  is an observation on a few dozen targets. It is not explained here, and it
  is not tested beyond these draws.

### G5: Möbius pole, `m = 2` (4 random `V` × 4 class-0 unsat targets per arm)

| cell | `x_subspace` | `β = 0` | `β = 1` | `β` random (pooled, 48) | `β`-spread | random-`β` arms' own spread |
|:--|--:|--:|--:|--:|--:|--:|
| `(7, 3)`, `d` 8 | 2.750 | 2.750 | 2.938 | 2.833 | 0.188 | 0.188 |
| `(9, 4)`, `d` 8 | 3.125 | 3.125 | 3.125 | 2.979 | 0.146 | 0.062 |
| `(11, 5)`, `d` 7 | 4.000 | 4.000 | 3.938 | 3.979 | 0.062 | 0.062 |
| `(13, 6)`, `d` 6 | 4.000 | 4.000 | 3.812 | 4.000 | 0.188 | 0.000 |

- **`β = 0` reproduces the `x_subspace` arm draw for draw**, in degrees and
  satisfiable fractions (derived identity, checked by hand). At `b = 1`,
  `(u₁u₂)²·S₃(1/u₁, 1/u₂, x₃) = S₃(u₁, u₂, x₃)`, so the pole at 0 is the same
  system. This serves as an internal identity check of the G5 builder.
- **`β = 1`, the A14 frame.** It is the lowest arm at `(11, 5)` and
  `(13, 6)`, by 0.06 and 0.19. Both are under the 0.3 threshold. At
  `(13, 6)` the difference comes from 3 draws of 16 reading 3, against 4 for
  all other arms.
- **Engine substitution.** PROTOCOL §3 replaces the frozen
  `coordinate_quotients.rs` engine, which measures no refutation degree.

## Build provenance and cross-check

The measured records were produced by a binary built in the **shared**
target dir `/home/user/cargo-target` (sha256 `bbc02561…74ef`).

**The hazard.** The coordinator flagged that with a shared target dir, cargo
can treat the release library as fresh from another worktree's build. A
binary built "successfully" here may therefore link another lane's library
code.

**The rebuild.** The example was rebuilt from this worktree, at this branch's
source, in the lane-private target dir `/home/user/cargo-target-shape-search`
(sha256 `54feb2b2…2506`). The binaries differ.

**The cross-check.** A subset was re-run with the private build
(`crosscheck.sh`, output in `crosscheck.out`), and each record was compared
line for line with the committed one, with only the wall-time field `secs`
removed. Every record was **identical**:

| re-run | measured lines |
|:--|--:|
| G1 `(5,2)` | 39 |
| G1 `(7,3)` | 40 |
| G1 `(9,3)` control | 16 |
| G2 `(5,2)` | 78 |
| G2 `(7,3)` | 77 |
| G3 `(5,2)` | 36 |
| G3 `(7,3)` | 42 |
| G5, all four cells | 4 × 96 |
| G1 `(11,4)`, `gp` arm | 6 |

That makes **718 measured lines plus the per-instance summaries**. The
summaries carry the bases, targets and satisfiable counts.

**Not re-run:** G2 `(9,3)`, the rest of G1 `(11,4)`, and G3 `(11,4)`. That
was for time; they cost about 60 min of the shared host. They used the same
binary as the re-run cells, which reproduced exactly.

## Status of the shape-search note's §5 ("The frozen scale-8 run is in progress", L239)

The run never landed. This was checked read-only from the main checkout and
through the GitHub commits API on 2026-10-15.

- `experiments/factor_base_closure_sweep.json` does not exist on
  `origin/main`. The commits API lists no commit that ever touched that path.
- The note's freezing commit `110117d8c` (2026-10-01) added the example and
  five notes, and no JSON.
- The note's §5 still reads "Pending". Its §6 readings, and the FFD-ledger
  row P3 that cites the JSON, therefore rest on the uncommitted half-scale
  smoke run (seed 3, scale 0.5) that §5 itself names.
- No open or merged pull request among the 100 most recent mentions the
  run.

The status is **not completed; the frozen evidence file is missing**. This
lane did not re-run it. It is the `m = 2` sweep and was out of this brief.

## Reading

- **At `m = 3` and these sizes, the refutation degree is set by `l` alone.**
  - It is 5 at `l = 2` and 6 at `l = 3, 4` (`(11,4)` reads 6 on all 30 G1 draws).
  - Shape (`gp`, `fp`, `coord`, random), modulus (G3) and Frobenius twist
    (G2) move it by exactly zero.
  - This matches `RESEARCH_DREG_MEASUREMENT.md` Result 10's `⌈l/2⌉ + 4`.
- **The subfield is the one shape that moves it**, by one degree at
  `(9, 3)`.
  - Its class-0 decomposable fraction there is 0.000 (0 of 8 class-0 targets drawn). At
    `b = 1` its points form the subgroup `E(F_8)`.
  - That is the yield collapse §6.5 of the note records at `m = 2`, now
    seen at `m = 3`.
  - It does not exist at prime `n`.
- **The degree is quantized and constant within a cell.** So "no gap" here
  means "the same integer on every draw". The statistic has no spread to
  hide a 0.5 effect at these `l`.
  - The open question is `l ≥ 5`, where the ladder reads `≥ 7` and could
    separate shapes. This host cannot reach it (frozen cells, above).
  - At `l = 2` every 2-dimensional subspace has `dim V·V = 3 = 2l − 1`. So
    `(5, 2)` cannot separate `gp`'s Kneser-extremal closure from random in
    any case.
- **Nothing here bears on ECC2K-130** beyond what the note's §2.1 and §8
  already say.

## Skipped, and why

| item | why skipped |
|:--|:--|
| G1 frozen cells | resources (above) |
| G4 | not built, for time; it needs a new base and a new system |
| G6 | no survivor |
| G2 at `(11, 4)` | budget |
| `KIC_*` engine path | the engineering path is documented as outcome-invariant and was not re-verified here. The 2.5 GB cap was respected through `KIC_SPARSE_DENSE_BUDGET_MB=1200` |

## Files

| path | contents |
|:--|:--|
| `PROTOCOL.md` | the pre-registration |
| `runs/*.jsonl` | one line per measured target, plus per-instance summaries. Each has an `.err` and a `.meta` (exit code, wall time, sampled peak RSS) |
| `summary.md` | the `score` output over every run file |
| `shape_search_degrees.svg`, `.png` | the graph, from `plot.py`, which plots `summary.md` values only |
| `report.typ`, `report.pdf` | the PDF report |
| `crosscheck.sh`, `crosscheck.out` | the private-build cross-check |

Reproduce:

```sh
cargo build --release --example factor_base_shape_ladder
B=target/release/examples/factor_base_shape_ladder
export F4_F2_MAX_ROWS=50000000 F4_F2_MAX_COLS=50000000 KIC_SPARSE_DENSE_FINISH=1 KIC_SPARSE_F5=1 KIC_SPARSE_DENSE_BUDGET_MB=1200
$B selftest
$B g1 --cell 7:3:7 --shapes random,random2,gp,fp,coord --instances 4 --targets 2 --seed 20261015
# … every command is in runs/*.meta ("args")
$B score runs/*.jsonl > summary.md
```
