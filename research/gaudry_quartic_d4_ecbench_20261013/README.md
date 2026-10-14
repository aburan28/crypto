# Gaudry's `k = 4` index calculus with the `D₄` oracle, end to end in ecbench

2026-10-13.  Protocol: [`PROTOCOL.md`](PROTOCOL.md) (registered before the
session).  Spec: `specs/SPEC-k4-p29.json`; sealed session `sessions/k4-p29/`;
audit `audits/k4-p29.audit.json` (**OK: 3 records, 3 verified, 3 replays
reproduced, 0 problems**).  Report: [`REPORT.pdf`](REPORT.pdf).
**Class: engineering.**

## Result

The whole method, every phase charged, recovers the planted logarithm on
`E(F_{29⁴})`, `#E = 708,208 = 16 · 44,263` (its ICV1 id is in
`sessions/k4-p29/plan.json`; the curve is not in the current registry), one target:

| arm | method | S | isolation |
|:--|:--|--:|:--|
| `rho-neg` (reference) | `rho.negation` | 1.065 | L2 |
| `ic-k4-d4` (candidate) | `ic.pipeline`, `fghr-y-line`, `gaudry-d4:macaulay=14` | **7.17·10⁸** | L1 |
| `rho-neg-aa` (control) | `rho.negation` | 1.065 | L2 |

**Index calculus costs 6.7·10⁸× rho here**, within 2× of §11.16's model
with the measured inputs (`1.3·10⁹`) and inside the registered range.
`S` is a count: the oracle's `F_p` multiplications priced at the curve's 97
per addition (`priced_by = counted`), no wall-priced term.  It omits the
work of the 9 degenerate systems (under 1.4 %, below).

| phase | group-addition equivalents | share |
|:--|--:|--:|
| factor base (7 columns, 13 lifts) | 424 | 3·10⁻⁷ % |
| relations (633 trials, 631 `D₄` solves) | 1.509·10¹¹ | ≈ 100 % |
| linear algebra (8 × 8, rank 8, 126 row operations) | 1.3 | 9·10⁻¹⁰ % |
| verification | 17 | 1·10⁻⁸ % |

Per solve: `2.32·10¹⁰` `F_p` multiplications (elimination 76 %, normal forms
23 %, characteristic polynomial 0.6 %), the cost `../gaudry_quartic_d4_20261013`
measured.  Wall time 2.0 h on one reserved core (L1: not admitted as a
timing; `S` does not depend on it).

## Predictions (PROTOCOL.md)

- E1 (verified; audit replays identical): **held**, 3 of 3.
- E2 (600–4,000 trials): **held**, 633.  The rate the range was derived
  from did not: the protocol assumed 0.003–0.008 decompositions a trial
  (uniform targets), and the walk's subgroup targets decomposed at 11 in
  633 (0.017; 8 in 633 kept as relations).  The trial range held because
  the method needed fewer relations than assumed; the stated rate missed.
- E3 (`S_IC` in `[3·10⁸, 10¹⁰]`, ratio in `[10⁸, 10¹⁰]`, model `1.3·10⁹`):
  **held**, `7.17·10⁸`, ratio `6.7·10⁸`.
- E4 (relations ≥ 99 % of `S`): **held**.

## What the run also showed

- **Lost relations.** 11 trials decomposed; 3 of them (`lift_failures`)
  had a summand outside the base — most likely the point with `Y = 0`,
  which `fghr-y-line` leaves out (its `T`-translate is its negative, so
  `[h]P = O` and it contributes nothing to a row).  Dropping such a summand
  instead of the relation would have saved about a quarter of the trials.
- **Degenerate systems, and an accounting gap.** 9 of 631 (1.4 %) did not
  close at degree 14, as in the `p = 269` measurement (1 in 20); a trial is
  discarded.  The oracle as recorded charged nothing for those nine
  attempts, so **`S` is understated** by their cost: at most nine solves'
  worth, `≤ 9 · 2.32·10¹⁰` of `1.46·10¹³` `F_p` multiplications, under
  1.4 %.  `gaudry-d4:…,charge_unsolved=1` charges them; it is off by
  default so this session's replays keep their counts (found in review on
  #1792).
- **Isolation.** The two-hour IC run graded L1, the rho runs L2.

## Amendment 1: the `Y = 0` point back in the base

`fghr-y-line:zero=1` adds the order-4 point with `Y = 0` as its own column
(its logarithm is 0: `[h]P₀ = O`).  Measured first by exact enumeration on
24,000 uniform targets (`raw/hit_rate_zero.jsonl`): 136 decompose with it,
96 without, **1.42×**, 29 % ± 8 % of hits needing it.  Then end to end on
the same workload and seeds (`sessions/k4-p29-zero/`, audit **OK**, 2 of 2
replays identical):

| | first session | with `Y = 0` |
|:--|--:|--:|
| trials | 633 | **574** |
| relations (columns + target) | 8 (8) | 9 (9) |
| decompositions lost (`lift_failures`) | 3 | **0** |
| degenerate systems | 9 | 9 |
| `S` | 7.17·10⁸ | **6.49·10⁸** (−9.5 %) |
| ÷ rho | 6.7·10⁸ | 6.1·10⁸ |

Predictions (amendment): Z1 (verified, identical, no lost decompositions)
**held**; Z2 (fewer trials, 290–610) **held** at 574, above the ~450 point
estimate because the extra column needed one more relation; Z3 (`S` in
proportion to the trials) **held** in proportion, but the ~`5·10⁸` figure
built on 450 trials **missed**.  One target cannot resolve a 1.4× rate;
the 24,000-target count is the measurement.

## Scope

One curve, one target, `p = 29`: the method working and priced end to end
at a size where it was expected to lose by about nine orders of magnitude.
§11.16's projection (crossover near `2^{135}` with the `D₄` oracle,
`h = 2`) is supported here only in its relation-phase constant, which this
run confirms within 2×.

## Code

`ecbench`: curve construction `fp4_torsion`, `Instance::Fp4`,
`curve_id::extension` (equal to `scripts/curve_id.py`'s `extension_id`,
pinned by a unit test); `CountedGroup` for `E(F_{p⁴})`; `ic.pipeline` on it
with `fghr-y-line` and `gaudry-d4`; counted pricing of field
multiplications in `ic_framework` (`priced_by = counted`) and in the
report.  Existing methods' accounting is unchanged.

## Reproduction

```sh
cargo build --release --bin ecbench
target/release/ecbench run --spec research/gaudry_quartic_d4_ecbench_20261013/specs/SPEC-k4-p29.json --out OUT --cpus auto --wait
target/release/ecbench verify --dir OUT --replay-all
```

Environment: Claude Code cloud container, Intel Xeon @ 2.10 GHz, 4 vCPU,
15 GiB RAM, rustc 1.97.0; binary `7a69171c14f656d5…`.  Conductor's control
plane was unreachable: no task card.
