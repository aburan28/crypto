# ECC2K-130: what a growing solving degree means for the Weil-descent route

2026-09-30. Class: **accounting**. It is a derived reading of measurements
already committed plus one new measurement stage, the degree-7 study
(`research/dreg_degree7_ell5_20260930/`). It adds no algorithm, no `S` and
no speed claim, and it moves no ECC2K-130 cost. Every number at
`n = 131` below is **derived, not measured**, and is marked with the
exponents it rests on.

**Updated 2026-10-06** with the first `m = 4` measurement
(`research/dreg_m4_chain_20261006/`). It adds §5 and point 5 of the short
version, plus the new rows of §2's table. The cost model is now native
(`examples/descent_degree_cost`), and its first 16 rows reproduce the
2026-09-30 table exactly.

## The short version

1. **Measured.** The system is the `m = 3` chained `S₃` system over
   random `ℓ`-dimensional subspaces of GF(2^n), with `n ≤ 15`.
   - Its **first fall degree is 3 everywhere**. That is Kosters–Yeo's trace
     equation (`RESEARCH_KOBLITZ_SCALING_TARGET.md`, H1 attribution).
   - **The degree at which the Macaulay matrix refutes it grows** over
     `ℓ = 2, 3, 4, 5`: 5, then 6, then 6–7, then **exactly 7** at `(10, 5)`,
     one degree above 6. That is the degree-7 study.
   - At fixed surplus it grows in every pair of the registered ladder.
2. **It tracks a generic reference.** In every cell with exact values, the
   refutation degree is within one of the semi-regular degree of
   regularity of the system's shape, and equal to it in 10 of 14. The
   reference predicted the 7 at `(10, 5)` before it was measured, and the
   `ℓ = 5` lower bounds are consistent with it. That reference grows about
   linearly in the number of unknowns.
3. **ECC2K-130's verdict does not depend on it.**
   - At `m = 3` the route loses to rho even with a free oracle:
     `2^68.58` against `2^60.81`.
   - At `m = 4` it needs an oracle that costs under `2^5.71` operations a
     target. The descended system has more unknowns than that before any
     elimination starts.
   - This is `RESEARCH_ECC2K130_DECOMPOSITION.md` §5.3, and it stands.
4. **What the growth removes is the asymptotic case.**
   - At a constant degree the per-target oracle is polynomial in the
     unknowns. The published sub-exponential estimates rest on that, with a
     crossover near `n ≈ 1250` (Kousidis–Wiemers, as quoted in
     `RESEARCH_ECC2K130_IC_LITERATURE.md`).
   - If the degree tracks the semi-regular reference, the oracle is
     exponential in the unknowns, and that crossover does not exist.
   - At `n = 131`, `m = 3` then prices at `≥ 2^196.95`. That is above an
     exhaustive search of the `2^129` subgroup. This is an extrapolation
     from `n ≤ 15`.
5. **`m = 4`, measured (2026-10-06, §5).** This is the decomposition size
   whose free-oracle floor is below rho.
   - It refutes at `ℓ + 4` for `ℓ ≤ 2`, one degree above `m = 3` at the
     same `ℓ`, and at ≥7 for `ℓ = 3`.
   - Like `m = 3`, it is flat in `n` at fixed `ℓ`. So the measured degree
     is set by `m` and `ℓ`, and the semi-regular reference overstates its
     growth in `n`.
   - **Any degree that grows with `ℓ` collapses the route.** Under either
     small-`ℓ` fit the best cell sits at the smallest factor base in the
     range, where the route is a search over targets:
     - `m = 3`: `2^158.96`;
     - `m = 4`: `2^173.13`.
     Both are above an exhaustive search of the subgroup. These are
     extrapolations of three-rung fits.

## 1. What is measured

Refutation degree of the `m = 3` chained `S₃` system (`b = 1`), per cell
`(n, ℓ)`: four unsatisfiable draws each, `N = n + 3ℓ` unknowns, `n` cubic
and `n` quadratic equations, surplus `S = n − 3ℓ` (equations minus
unknowns). A value is exact where the degree-`D` Macaulay row space first
contains `1`; `≥D` is a lower bound (the full degree-`(D−1)` matrix does
not). Sources: `RESEARCH_DREG_MEASUREMENT.md` Results 4–9, the last being
the degree-7 study.

| `ℓ` | cell `(n, ℓ)` | `N` | `S` | measured, four draws | semi-regular `D_reg` |
|--:|---|--:|--:|---|--:|
| 2 | `(5,2)` `(7,2)` `(8,2)` `(10,2)` `(12,2)` | 11–18 | −1 to +6 | 5 in every draw | 5 |
| 3 | `(4, 3)` | 13 | −5 | 5 5 5 5 | 5 |
| 3 | `(5, 3)` | 14 | −4 | 6 6 6 6 | 5 |
| 3 | `(7, 3)` | 16 | −2 | 6 6 6 6 | 6 |
| 3 | `(9, 3)` | 18 | 0 | 6 6 6 6 | 6 |
| 4 | `(7, 4)` | 19 | −5 | 7 7 6 7 | 6 |
| 4 | `(8, 4)` | 20 | −4 | 6 7 6 7 | 6 |
| 4 | `(11, 4)` | 23 | −1 | 6 6 6 6 | 7 |
| 4 | `(13, 4)` | 25 | +1 | 6 6 6 6 | 7 |
| 5 | **`(10, 5)`** | 25 | −5 | **7 7 7 7** | 7 |
| 5 | `(11, 5)` | 26 | −4 | ≥7 ≥7 ≥7 ≥7 | 7 |
| 5 | `(13, 5)` | 28 | −2 | ≥7 ≥7 ≥7 ≥7 | 7 |
| 5 | `(15, 5)` | 30 | 0 | ≥7 ≥7 ≥7 ≥7 | 8 |

The last column is `D_reg`, the semi-regular reference: the first
non-positive coefficient of `(1+z)^N / ((1+z²)^n (1+z³)^n)`. It comes from
`research/dreg_degree7_ell5_20260930/semireg.py`, a formula, not a run.

- **The degree never falls as `ℓ` rises at a fixed surplus, and it rises
  from `ℓ = 3` to `ℓ = 5` in every column where both are measured.**
  - `S = −5`: 5 to 7. `S = −4`, `−2` and 0: 6 to ≥7.
  - At `S = −5` it is flat from `ℓ = 4` to `ℓ = 5`: 7 7 6 7, then 7.
- **At fixed surplus it grows in every registered pair**, two of them with
  `ℓ ≥ 3` at both ends (Results 7 and 8).
- **The first fall degree is 3 throughout.** So the gap between the fall
  degree and the refutation degree widens from 2 at `ℓ = 2` to
  4 at `ℓ = 5`, where the measured degree is 7.
  - The first-fall-degree assumption behind the sub-exponential estimates
    (Petit–Quisquater, *On polynomial systems arising from a Weil descent*,
    ASIACRYPT 2012) treats that gap as bounded.
  - On this system, at these sizes, it is not.
- **Against the reference, the measured degree is within one degree in
  every cell with exact values.**
  - It is equal in 10 of those 14 cells. One of them is `(10, 5)`, where the
    reference predicted 7 before the degree-7 run.
  - It is one above in `(5, 3)` and `(7, 4)`, and one below in `(11, 4)` and
    `(13, 4)`.
  - The lower bounds of the three unrun `ℓ = 5` cells are consistent with
    it.

## 2. The price at `n = 131`, by degree

`examples/descent_degree_cost` evaluates the cost model of
`RESEARCH_ECC2K130_DECOMPOSITION.md` §5.1 with the oracle priced as a
Macaulay elimination of the descended system:

- It is the native port of `descent_degree_cost_20260930.py`, which
  produced this table on 2026-09-30 and stays as that record.
- The port's first 16 rows are identical, and its output is
  `descent_degree_cost-output-20261006.txt`.

The model:

- **relations:** `2^ℓ`;
- **targets:** `m!·2^{n−(m−1)ℓ}`, at least one a relation;
- **linear algebra:** `m·2^{2ℓ}`;
- **oracle:** `C(N, ≤D)^w` a target, the column count of the degree-`D`
  matrix raised to `w`.
  - `N = (m−2)·n + m·ℓ` unknowns, since the chain has `m − 2` intermediate
    points.
  - `w = 1` prices one touch per column, which no elimination beats.
  - `w = 2` is the usual estimate at `ω = 2`.

Each row is minimised over integer `ℓ`. The rho reference is `2^60.81`
(§5.2 of that note, with `⟨−1⟩ × ⟨π⟩`).

| `m` | degree `D` | `w` | `ℓ` | `N` | oracle a target | total | vs rho |
|--:|---|--:|--:|--:|--:|--:|--:|
| 3 | free oracle (the floor) | — | 33 | 230 | — | `2^68.58` | `2^+7.77` |
| 3 | 2, any Macaulay matrix | 1 | 37 | 242 | `2^14.84` | `2^76.12` | `2^+15.31` |
| 3 | 6, constant | 1 | 43 | 260 | `2^38.59` | `2^88.05` | `2^+27.24` |
| 3 | 6, constant | 2 | 45 | 266 | `2^77.58` | `2^122.58` | `2^+61.77` |
| 3 | 7, constant | 1 | 44 | 263 | `2^43.90` | `2^90.53` | `2^+29.72` |
| 3 | 7, constant | 2 | 45 | 266 | `2^88.02` | `2^133.02` | `2^+72.21` |
| 3 | **semi-regular `D_reg(N)`** (28 at `ℓ = 26`) | 1 | 26 | 209 | `2^115.37` | **`2^196.95`** | **`2^+136.14`** |
| 3 | semi-regular `D_reg(N)` (17 at `ℓ = 4`) | 2 | 4 | 143 | `2^144.31` | `2^269.90` | `2^+209.09` |
| 4 | free oracle (the floor) | — | 27 | 370 | — | `2^56.46` | `2^−4.35` |
| 4 | 2, any Macaulay matrix | 1 | 30 | 382 | `2^16.16` | `2^62.88` | `2^+2.07` |
| 4 | 2, any Macaulay matrix | 2 | 33 | 394 | `2^32.50` | `2^69.64` | `2^+8.83` |
| 4 | 6, constant | 1 | 34 | 398 | `2^42.30` | `2^76.31` | `2^+15.50` |
| 4 | 6, constant | 2 | 34 | 398 | `2^84.59` | `2^118.59` | `2^+57.78` |
| 4 | 7, constant | 1 | 34 | 398 | `2^48.11` | `2^82.11` | `2^+21.30` |
| 4 | 7, constant | 2 | 34 | 398 | `2^96.21` | `2^130.21` | `2^+69.40` |
| 3 | `⌈ℓ/2⌉ + 4`, the `m = 3` fit (6 at `ℓ = 4`) — *new 2026-10-06* | 1 | 4 | 143 | `2^33.38` | `2^158.96` | `2^+98.15` |
| 3 | `⌈ℓ/2⌉ + 4`, the `m = 3` fit — *new* | 2 | 4 | 143 | `2^66.76` | `2^192.34` | `2^+131.53` |
| 4 | semi-regular `D_reg(N)` of the `m = 4` shape (55 at `ℓ = 27`) — *new* | 1 | 27 | 370 | `2^220.56` | `2^275.14` | `2^+214.33` |
| 4 | semi-regular `D_reg(N)` (38 at `ℓ = 5`) — *new* | 2 | 5 | 282 | `2^314.45` | `2^435.04` | `2^+374.23` |
| 4 | **`ℓ + 4`, the measured `m = 4` fit** (8 at `ℓ = 4`) — *new* | 1 | 4 | 278 | `2^49.55` | **`2^173.13`** | **`2^+112.32`** |
| 4 | `ℓ + 4`, the `m = 4` fit — *new* | 2 | 4 | 278 | `2^99.10` | `2^222.68` | `2^+161.87` |

- **The free-oracle rows reproduce §5.3.** At `m = 3`, `ℓ = 33` gives
  `2^68.58`, exactly. At `m = 4` the integer optimum `ℓ = 27` gives
  `2^56.46`, against §5.3's `2^56.40` at the continuous `ℓ = 26.83`.
- **Units.**
  - The oracle column counts monomial-column touches. Rho counts group
    operations, and one group operation over GF(2^131) is many such
    touches.
  - Converting would lower every row with an oracle by `log₂ K`, where `K`
    is touches per group operation. `K` is not measured here.
  - So the `m = 4`, `D = 2`, `w = 1` row, `2^+2.07`, is **not a separation**.
    Every other row with an oracle clears rho by 8 to 209 bits, and a
    reader should subtract `log₂ K` from those gaps.
- **The `m = 4` constant-degree rows are `m = 3`'s degrees applied to
  `m = 4`.** The chain has `2n` cubic equations. At equal `ℓ` it has
  `n + ℓ` more unknowns, one more 131-bit intermediate point and one more
  summand.
  - On 2026-09-30 these rows were unmeasured and called the optimistic
    case.
  - §5 now measures `m = 4`: 6 at `ℓ = 2` and ≥7 at `ℓ = 3`. The `D = 6`
    and `D = 7` rows are already exceeded at `ℓ = 3`, far below the
    `ℓ ≈ 34` they are priced at, so they are confirmed as optimistic.
- **The fit rows' minima sit at the bottom of the search range,
  `ℓ = 4`.**
  - Once the degree grows with `ℓ`, a larger factor base costs more than
    it saves in targets. The optimiser therefore shrinks the base until
    the target count, `2^125.6` here, is nearly the whole group.
  - The route has then turned into an expensive search over targets.

### What the table says

- **The verdict is degree-independent.**
  - `m = 3` is above rho at every degree, including a free oracle.
  - At `m = 4`, the lowest Macaulay degree there is (`D = 2`) already
    spends the whole margin a free oracle had. Its `2^−4.35` becomes
    `2^+2.07`.
  - That is the §5.3 budget of `2^5.71` a target, seen from the oracle's
    side. It holds at any degree.
- **The degree decides how far past rho the route is, and whether it has
  an asymptote at all.**
  - At a constant degree the oracle column is `N^D/D!`, polynomial in
    `n`. That is the regime of the published crossover estimates: `2^86.0`
    at `n = 131` by Kousidis–Wiemers' formula, "about `2^25` short"
    (`RESEARCH_ECC2K130_IC_LITERATURE.md` §2).
  - If the degree tracks `D_reg(N)`, which rises linearly in `N`, the
    column is exponential in `N`. The best `m = 3` cell is then `2^196.95`,
    above both rho and exhaustive search of the subgroup (`2^129`).
  - At that point the descent is not a slow sub-exponential method. It is
    worse than brute force.
- **Which regime the data support.** The measurements at `n ≤ 15` rule out
  constant degree 6, which is what `ℓ = 3` and most of `ℓ = 4` show.
  They do not yet separate a degree that has levelled off at 7 from one
  that follows `D_reg`. Both say 7 at `(10, 5)`.
  - Between those two, the `n = 131` price differs by 106 bits at `w = 1`.
  - `(15, 5)` is the cell that decides it.
  - *Added 2026-10-06:*
    - Both `m` show a degree that is flat in `n` at fixed `ℓ` and rises
      with `ℓ`.
    - That favours an `ℓ`-driven degree, which the fit rows model, over
      both "levelled off" and the `N`-driven `D_reg`.
    - The fit rows sit between those two prices, and above brute force.
    - Whatever the right model, any degree that keeps rising with `ℓ`
      gives the same verdict: the optimiser shrinks the factor base until
      the route is a target search.

## 3. What changes, and what does not

**Does not change.**

- The ECC2K-130 cost verdict. Index calculus by summation polynomials
  loses to rho at `n = 131` by the product law and the oracle budget
  (`RESEARCH_ECC2K130_DECOMPOSITION.md` §5.1–5.3). The scoreboard already
  says so.
- What remains open is also unchanged: an oracle that is not an
  elimination of the descended system at all (§5.4 of that note). This
  note prices Macaulay-type oracles, not every oracle.

**Changes.**

- **The constant-degree premise, for this system.**
  `RESEARCH_IC_BOUNDARY_LEDGER.md` records "the solving degree is flat in
  `n`" at `m = 2`. At `m = 3` it is not flat.
  - The chained `S₃` system's refutation degree rises with `ℓ`.
  - It rises at fixed surplus, across the four registered ladder pairs.
  - The first fall degree stays at 3.
  - H1 of `RESEARCH_KOBLITZ_SCALING_TARGET.md`, "the fall degree stays
    ≤ 3", therefore holds. It is Kosters–Yeo's trace equation.
  - H1's falsifier was a proxy for "the subspace-restricted systems are not
    as benign as the first data suggests". The fall degree does not show
    that. The refutation degree, which H1 does not measure, does.
- **The model to extrapolate with.**
  - A constant degree is contradicted at these sizes.
  - The semi-regular reference is the one model on hand that fits every
    measured cell within one degree, and it is not constant.
  - Any future ECC2K-130 figure for this route should state which degree
    model it assumes. The two differ by more than 100 bits at `n = 131`.
  - *2026-10-06:* §5 adds a third model, a degree set by `m` and `ℓ` and
    flat in `n`. It fits both `m` measured so far, and it is the one the
    data favour.

## 4. What would move it

- **A degree that stops growing.**
  - `(15, 5)` is the discriminating cell: `D_reg` is 8 there, and every
    other measured `ℓ = 5` cell has a reference of 7.
  - Its degree-7 matrix is 3.1M × 2.8M. At `(10, 5)`'s measured band-7
    rank ratio its dense block is about 110–120 GB, so it needs a host of
    about 128 GB.
  - A resolution at 7 there would put `ℓ = 5` at 7 wherever it is measured,
    and would favour the levelled-off reading.
  - A ≥8 would follow `D_reg` into a fifth consecutive rung.
- **A Frobenius-stable factor base.** It would divide the target count
  by `n`. At `n = 131` the only stable subspaces have dimensions
  0, 1, 130 and 131, because `ord₁₃₁(2) = 130`
  (`RESEARCH_ECC2K130_DECOMPOSITION.md` §6). No such base exists at the
  `ℓ` this route needs.
- **A different oracle.** SAT is measured and goes the wrong way (§5.4
  there). Anything else is unpriced.

## 5. The `m = 4` measurement (2026-10-06)

`research/dreg_m4_chain_20261006/` and Result 10 of
`RESEARCH_DREG_MEASUREMENT.md`. It was pre-registered, and its predictions
held.

| `ℓ` | `m = 4` cells, unknowns | measured, four draws each | semi-regular `D_reg` | `m = 3` at the same `ℓ` |
|--:|---|---|---|---|
| 1 | (4,1) (5,1) (6,1), 12–16 | 5 in every draw | 5 | — |
| 2 | (5,2) (6,2) (7,2) (8,2), 18–24 | 6 in every draw | 6, 7, 7, 7 | 5 in every draw |
| 3 | (9,3) (10,3), 30 and 32 | ≥7 ×4; M4_10_3 | 8, 9 | 6 (5 at surplus −5) |

- **`m = 4` sits one degree above `m = 3` at the same `ℓ`**, and the gap
  does not close as `n` grows.
- **It rises with `ℓ`:** 5, then 6, then ≥7. That is at least one degree
  per unit of `ℓ`, against about half a degree for `m = 3`.
- **It is flat in `n` at fixed `ℓ`.**
  - So the measured degree is set by `m` and `ℓ`.
  - The semi-regular reference, which counts unknowns and equations, is
    within one in every cell but overstates the growth in `n`.
- **What it does to §2.**
  - The `m = 4` constant-degree rows are now known to be optimistic: 7 is
    exceeded at `ℓ = 3`, against the `ℓ ≈ 34` they assume.
  - The new `ℓ + 4` row prices `m = 4` at `2^173.13` at `w = 1`, with its
    minimum at the smallest factor base. That is `2^+112` over rho, and
    above an exhaustive search of the subgroup.
  - It is an extrapolation of a fit to three rungs, one of them a bound.
    What it shows robustly is the shape: a degree that rises with `ℓ` makes
    the factor base shrink until the route is a search.
- **Unchanged:** the verdict against rho, which never depended on the
  degree (§2, "What the table says").

## Scope and disclosures

- **System.** `m = 3`, the chained `S₃`, with `b = 1`.
  - `b = 1` is the ECC2K-130 family's own coefficient: `E₀: y² + xy = x³ + 1`.
  - Random `ℓ`-dimensional subspaces, four unsatisfiable draws per cell.
  - The degree is the refutation degree of bounded Macaulay linear algebra
    (the constant `1` in the row space). It is neither an F4/F5 solving
    degree nor a certified first fall degree
    (`degree-reporting-correction-20260922` on the scoreboard).
- **Fields (`AGENTS.md` §8b).**
  - The measured degrees are `n = 4, 5, 7, 8, 9, 10, 11, 12, 13, 15`.
  - Several have proper intermediate subfields over GF(2): `n = 4`
    (GF(2²)), 8 (GF(2²), GF(2⁴)), 9 (GF(2³)), 10 (GF(2²), GF(2⁵)), 12
    (GF(2²), GF(2³), GF(2⁴), GF(2⁶)) and 15 (GF(2³), GF(2⁵)).
  - The subspaces are random, so no subfield structure is chosen or used.
  - The `m = 4` cells (§5) are `n = 4`–`10`, with GF(2^6) adding GF(2²)
    and GF(2³). Their other subfields are listed above.
  - The challenge field GF(2^131) has no proper intermediate subfields.
  - No result here depends on a subfield. No Frobenius action is used.
- **§8a.** Nothing here is an index-calculus improvement, so the `m = 83`
  gate is neither claimed nor discharged. No statement here is evidence
  at `m = 83` or `m = 131`. The `n = 131` rows are extrapolations.
- **Extrapolations, named.**
  - Every `n = 131` figure in §2 is derived.
  - The constant-degree rows assume the degree stops where `n ≤ 15` left
    it.
  - The semi-regular rows assume the one-degree agreement seen at `n ≤ 15`
    persists to `N ≈ 209`. For `m = 4` that is `N ≈ 370`.
  - The fit rows (2026-10-06) assume that `⌈ℓ/2⌉ + 4` (`m = 3`) and
    `ℓ + 4` (`m = 4`), fitted at `ℓ ≤ 5` and `ℓ ≤ 3`, hold at `ℓ = 4`.
    That `ℓ` is in range, but `n = 131` is not.
  - None of these is measured.
- **Units.** Oracle costs are in monomial-column touches and rho is in
  group operations, unconverted (§2).

## Reproducing

```sh
# native (AGENTS.md: no Python); the 2026-09-30 Python tools stay as that
# round's record and are reproduced exactly
cargo run --release --example dreg_score -- semireg 3 10:5 13:5 15:5     # the reference column
cargo run --release --example dreg_score -- score research/dreg_degree7_ell5_20260930/runs/cell-10-5-7.u*.jsonl
cargo run --release --example dreg_score -- score research/dreg_m4_chain_20261006/runs/cell-*.jsonl  # §5
cargo run --release --example descent_degree_cost                        # §2's table
```
