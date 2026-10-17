# Is the `m = 4` tree-split set cut out at low degree on a structured factor base? Results

This file was written after the run. [PREREGISTRATION.md](PREREGISTRATION.md) is unchanged
since its registration commit `3dbbfa1d`, which was pushed before any registered cell ran.
The supporting records:

- run report: [reports/run-20261011T2105Z/](reports/run-20261011T2105Z/) (report, graph,
  PDF);
- readout: [runs/registered/readout.txt](runs/registered/readout.txt);
- records: [runs/registered/cells.jsonl](runs/registered/cells.jsonl);
- Lean: [lean/](lean/);
- cairn log: [cairn/log.jsonl](cairn/log.jsonl).

## Registered verdict: **CLOSED** for `poly`, `frob` and `rand`; validity gate **PASS**

**Structure in the factor-base subspace does not lower the sampled cut-out degree of `T2`.**
At every cell of every arm, `D_s(curve)` is within one degree of a random set of the same
size, and it climbs with `ℓ` at the random sets' rate.

| arm | series | uncensored cells | gaps (`D_s(curve) − D_s(R)`) | OLS slopes of `D_s(curve)` on `ℓ` | verdict |
|:--|--:|--:|:--|:--|:--|
| `poly` | 6 | 30/30 | all 0 | 1.0, 1.0, 0.9, 0.7, 0.7, 0.8 | **CLOSED** |
| `frob` (`n = 21`) | 4 | 20/20 | ~~14 × 0, 6 × −1~~ 15 × 0, 5 × −1 (corrected 2026-10-16) | 0.9, 0.9, 1.1, 1.0 | **CLOSED** |
| `rand` (`n = 21`) | 2 | 10/10 | all 0 | 0.7, 0.8 | **CLOSED** (replicates the predecessor) |
| `frob`, `n = 17`, supplementary | 2 | 2/2 | 0, 0 | — | not in the verdict |
| `subfield` gate | 2 | 2/2 | `T2 ⊆ GF(2^7)`; degree-1 rank 8, deficit 14 in both (gate needs ≥ 14) | — | **PASS** |

The ALIVE threshold was a slope of `≤ 0.15` together with a gap of `≤ −2` at the top `ℓ`.
No series comes near either. No cell was censored and none failed (`cells.log`: 64 `ok`).

> **Correction (2026-10-16, sheriff case SHF-20261016-3caa6f, finding L4).** The registered
> readout has **five** `frob` (`n = 21`) cells at gap −1 and **fifteen** at gap 0, not six
> and fourteen as first written here and in the run report
> ([`runs/registered/readout.txt`](runs/registered/readout.txt), lines 19, 26, 51, 56, 58:
> `a = 0` seed 20261011 `ℓ = 4`; `a = 0` seed 20261012 `ℓ = 6`; `a = 1` seed 20261011
> `ℓ = 4`; `a = 1` seed 20261012 `ℓ = 4` and `ℓ = 6`). Six `frob` cells have `V` inside
> `GF(2^7)`; one of them (`a = 0`, seed 20261012, `ℓ = 4`, `|T2| = 1`) is degenerate at gap 0,
> as the section below already says, so five, not six, are at −1. The miscount is in the
> prose only: the records, `Data.lean` and the registered verdict (`Verdict.lean`, CLOSED)
> are unchanged. The original wording is kept, struck through, where it was wrong.

### The ~~six~~ five `frob` cells at gap −1 are subfield leakage (post hoc; count corrected 2026-10-16)

The registered `frob` arm excluded `g = t^d − 1`, which is a subfield itself. It did not
exclude a proper Frobenius-stable subspace of a subfield. Three of the ten `g` the run used
divide `t^7 − 1`: `23`, `29` and `127`. So their `V` lies inside `GF(2^7)`. This is proved
in [`lean/TreeSplit/Leak.lean`](lean/TreeSplit/Leak.lean), by kernel evaluation
(`decide +kernel`; corrected 2026-10-16: it was first closed ~~with `native_decide`~~, which
trusts the compiler and is not a kernel proof, sheriff case SHF-20261016-3caa6f, finding L3).
For such a `V`, `F ⊂ E(GF(2^7))` and `T2 ⊆ GF(2^7)`, which is a 7-dimensional linear subspace.

- **They are exactly the gap −1 cells.** These are `ℓ = 4` for both seeds and `ℓ = 6` for
  seed 20261012 (`|T2|` = 29, 28, 70, 29, 35). The sixth subfield cell, `a = 0`, seed
  20261012, `ℓ = 4`, is degenerate: `|F| = 2`, `|T2| = 1`, gap 0.
- **Every `frob` cell whose `V` is not inside a subfield has gap 0.**
- **Why it is not a lever.** Inside a subfield, decompositions exist only for targets in
  `E(GF(2^7))`, which is the positive control's effect at smaller size, not a new
  mechanism.
- **The verdict does not depend on these cells.** It is CLOSED with them, and it would be
  CLOSED without them.
- **For a follow-up**, a `frob` arm should exclude every `g` that divides some
  `t^d − 1` with `d | n`, `d < n`.

**Run provenance.**
- The run used commit `dbd878f0` (`runs/registered/commit`). Its harness is the
  registration commit's, byte for byte.
- After the run, `examples/tree_split_structured.rs` was put through `rustfmt`, which only
  re-wraps lines (`git diff -w` shows wrapping only).
- The rebuilt binary regenerates `readout.txt` and `lean/TreeSplit/Data.lean` byte for byte
  from `cells.jsonl`. CI repeats that check.

## What the Lean and cairn checks establish

Lean has no clock: none of this is about time. It checks that the reported result is the
one the registered rule gives on the recorded data, and that part of the recorded data is
reproduced by a second implementation.

| check | result | trust |
|:--|:--|:--|
| `Rule.lean` is the registered rule | sha256 `714601c1…760d`, the same as in §7 of the registration and in cairn objective `sha256:f236a122…` (posted before the run) | hash |
| `Verdict.lean`: `report cells = ⟨gate := true, poly := closed, frob := closed, rand := closed⟩` and `cells.length = 64` | proved | kernel `decide`. Changing `frob := .closed` to `.alive` fails the build |
| `Replay.lean`: an independent implementation recomputes `|T2|`, `D_s(curve)` and `D_s(R)` for the 12 `poly` cells with `ℓ ≤ 5` | all 12 match the records exactly | `native_decide`, which trusts the Lean compiler. Changing one recorded `D_s` fails the build (checked during development) |
| `Leak.lean`: which `g` divide `t^7 − 1` | `[23, 29, 127]` | kernel, `decide +kernel` (2026-10-16; first ~~`native_decide`~~) |
| [`lean/AxiomAudit.lean`](lean/AxiomAudit.lean) (added 2026-10-16) | every theorem above depends only on `propext`, `Classical.choice` and `Quot.sound`, except `replay_ok`, which adds its own `native_decide` compiler-trust axiom and nothing else | the audit fails the run on any other axiom (`sorryAx` included) |
| cairn objective `sha256:2fb64f9a…` (result): the same theorem as `Verdict.lean`, over `Rule.lean` and `Data.lean` | **accepted** by cairn's Lean jail, kernel | `cairn audit` re-verifies it |
| cairn objective `sha256:d1fd6832…` (replay) | **rejected**: cairn accepts kernel proofs only, and refuses `native_decide`. Recorded as is | — |
| CI: [`.github/workflows/tree-split-structured-v-lean.yml`](../../.github/workflows/tree-split-structured-v-lean.yml) | Regenerates `Data.lean` and the readout from the records and diffs them, checks the rule hash, and runs `lake build` (all of the above) on every change to the study; since 2026-10-16 also the axiom audit (log kept as an artifact), `leanchecker --fresh`, and a check that `native_decide` closes `replay_ok` only | — |

**cairn log** (`cairn/log.jsonl`, 19 entries, cairn 6.1.0):
- It was built locally with `CAIRN_EPOCH_SECONDS=5`, so it audits with that setting.
  `CAIRN_EPOCH_SECONDS=5 cairn --log cairn/log.jsonl audit` reports "chain intact, every
  settled claim re-verified".
- It holds one rejected registration attempt. The proof text lacked its leading `:=`, and
  the accepted resubmission follows it.
- Rewards are 0. The log is a timestamped, re-verifiable record, not a payment.

## Scope and limits

- **What was tested.** Two Koblitz curves, `n ∈ {17, 19, 21}`, `ℓ ≤ 8`, `T2` only (not
  `Z`), with a sampled cut-out degree: 2,000 samples, capped at `D = 6` and 90,000
  monomials. The structured families were the polynomial-basis window and Frobenius-stable
  kernels at `n = 21`.
- **What was not tested.** Normal-basis windows, other low-expansion bases, `Z ∩ A_R`,
  linearisation with auxiliary variables, and larger `n` or `ℓ`.
- **Not timed.** The metric is a degree. The `secs` fields are logged, not measured.
- **Claims.** No exponent, speedup or deployed-curve claim is made. The tree-split
  reopening condition of
  [DECISION-20260929-stop.md](../ic_candidate_tournament_20260915/campaign_20260916/DECISION-20260929-stop.md)
  stays unmet for these structured families.

## Graphs

- [`reports/run-20261011T2105Z/ds_vs_ell.svg`](reports/run-20261011T2105Z/ds_vs_ell.svg)
  plots `D_s` against `ℓ`, curve against control, per series, for `poly`, `frob` and `rand`.
- **Canonical graphs checked and left unchanged:**
  - `docs/index-calculus-scoreboard.html`, `docs/ic/progress-timeline.json` and
    `docs/ic/LEADERBOARD.md`. This study prices no attack against rho and changes no
    boundary row. It is a degree screen.
  - The predecessor study is not on any canonical graph (`rg ic_tree_split_cutout_20260930`
    over `docs/`: no hits), and neither is this one.
