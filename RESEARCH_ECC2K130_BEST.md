# ECC2K-130: the best index-calculus algorithm this repository has, and the distance to rho

**Question.** Of everything measured here against ECC2K-130
(`K_0 : y² + xy = x³ + 1` over `F_{2^131}`, prime subgroup order
`r ≈ 2^{129}`), which index-calculus pipeline is cheapest in the one
unit every ECC2K-130 row shares, how far is it from Pollard rho, and
can an `F_2` SAT solver as the decomposition oracle close that gap?

**Answer.** The best measured pipeline is the G7e stack: type-II ONB ×
Hamming-weight-2 even-trace factor base × Frobenius-orbit collapse × a
stored pair table, at `|F| = 8384`. It costs `2^{116.84}` field
products, `2^{53.62}×` rho. Its counting identity is tight (ratio to
the pair-table identity is 1), so no constant-factor engineering on it
changes the exponent. The only way to parity is an oracle whose
per-call cost `Q = 2^{αl}` has `α ≈ 0.4–0.6` at `l = 45`, depending
on the frame (§2); anything at `α ≥ 1` stays `2^{21}×` rho or worse.
The `F_2` SAT route is measured in §4 on a five-rung ladder against
exactly that scale. WDSat, the published solver for this encoding,
fits **`α = 2.998` in conflicts** (`3.31` in product-equivalents) over
`l = 5..9`, 30/30 answers verified, one conflict per sorted triple
on every rejection: it enumerates, at a slower rate than the
pairs-and-solve oracle it was measured beside (`α = 2.13`), and the
gap widens from 215× to 5,100×. Pre-registered prediction `α ≈ 3`
confirmed; abandonment condition `α ≥ 1.9` met. **No, an `F_2` SAT
oracle does not make this pipeline capable of rho parity**, and the
distance is an exponent, not a constant. Class: measurement.
Scoreboard panel `#ecc2k130-best-20260919`.

This note re-runs nothing from earlier threads. §1–§2 cite frozen
artefacts; §3–§4 are one new measurement, protocol frozen before the run.

## 1. The best algorithm, highlighted

Cited from
[`ecc2k130/benchmarks/indexcalc-g7e/summary.json`](ecc2k130/benchmarks/indexcalc-g7e/summary.json)
and [`RESEARCH_ECC2K130_G7E_INDEX_CALCULUS.md`](RESEARCH_ECC2K130_G7E_INDEX_CALCULUS.md).
Unit: field products in `F_{2^131}`; affine addition priced at 10
products; rho at the measured 5.3125 products per iteration on the
G7e SKU. `S = products / √r`, `√r = 2^{64.5}`.

| Piece | What it is | Why it is in the stack |
|---|---|---|
| Curve and unit | `K_0/F_{2^131}`, `r ≈ 2^{129}`, rho `2^{60.809}` iterations = `2^{63.22}` products | reference, measured |
| Representation | type-II optimal normal basis | Frobenius is a coordinate rotation, squaring is free |
| Factor base | `{x : HW(x) = 2, Tr(x) = 0}`, 4,192 abscissae, `|F| = 8384` points | `σ`-stable for every weight; even trace is the odd-order filter; polynomial-basis weight is *not* `σ`-stable |
| Unknowns | `|F| / 131 = 64` Frobenius orbits | divides the relation count and the linear algebra by 131 |
| Oracle | stored pair table, `C(|F|, 2) = 35,137,344` sums, `|F|` probes per target | the `2^{13.03}` cut below streaming pairs at this `|F|` |
| Correctness | 8/8 planted decompositions recovered; 0 hits on the generator, matching the expected yield `3.6 × 10^{-29}` | verified |

| Variant | log₂ products | `S` | vs rho | Class |
|---|---:|---:|---:|---|
| Pollard rho `⟨−1⟩ × ⟨π⟩` *(reference)* | 63.22 | 0.411 | `2^{0}` | baseline |
| Product-law floor, streaming pairs, any `|F|` | 136.91 | `2^{72.4}` | `2^{+73.69}` | floor |
| + Frobenius collapse `/131` | 129.87 | `4.78 × 10^{19}` | `2^{+66.66}` | accounting |
| Streaming pairs, G7e, `|F| = 8384` | 129.87 | `4.78 × 10^{19}` | `2^{+66.66}` | engineering |
| **Pair table, G7e, `|F| = 8384`** | **116.84** | **`5.70 × 10^{15}`** | **`2^{+53.62}`** | engineering |

Everything else measured against this instance is either on the same
floor (GPU inner loops, SAT as a slower enumerator), does not recover a
logarithm (homogeneous relations), or is structurally closed at prime
`n = 131` (subspaces, covers, towers);
[`RESEARCH_ECC2K130_IC_SYNTHESIS.md`](RESEARCH_ECC2K130_IC_SYNTHESIS.md)
§3 and §6 hold the ledger.

## 2. What parity requires

Restated from the route-target note, before any new number is named.
With `m = 3` summands from a base of size `2^l`, a uniform target
decomposes with probability about `2^{3l} / (3! · #G)`, so charging
`|F| / 131` relations costs

```
G_IC ≈ (2^l / 131) · (3! · 2^{131} / 2^{3l}) · Q  =  6 · 2^{131 − 2l} · Q / 131
```

per-call oracle cost `Q` in field products, the Frobenius collapse
credited. Against rho's `2^{63.22}` products, at the natural cap
`l = 45` (where `2^{3l}/6 ≈ 2^{131}` and the yield saturates):

```
G_IC / G_rho  =  6 · 2^{131 − 2l} · Q / (131 · 2^{63.22})  =  Q · 2^{−26.67}   at l = 45
```

| Oracle | `Q` (products) | `α` in `Q = 2^{αl}` | `l` | `G_IC / G_rho`, this frame | route-target frame |
|---|---|---:|---:|---:|---:|
| Streaming pair enumeration | `C(2^l, 2) · 20` | 2 | any | `2^{+66.65}` | `2^{+64.74}` |
| Stored pair table, memory `2^{2l−1}` | `2^l · 10` | 1 | 45 | `2^{+21.65}`, `2^{89}` entries | `2^{+27.78}` |
| Measured pair table (§1) | `2^l · 10` | 1 | 13.03 | `2^{+53.62}` | — |
| **Parity** | **`2^{26.67}`** | **0.59** | 45 | `2^{0}` | `Q = 2^{17.2}`, `α = 0.38` |

The streaming row reproduces the frozen G7e `2^{66.66}` to the second
decimal, which is the check that this frame is the receipt's frame.
The right-hand column is
[`RESEARCH_ECC2K130_ROUTE_TARGETS.md`](RESEARCH_ECC2K130_ROUTE_TARGETS.md)'s
pre-registered scale table in its own unit (group operations, no
orbit collapse); the two frames differ by constants (`2^{7.03}` for
the collapse, `2^{2.41}` products per rho iteration, `2^{3.32}` per
affine addition), and a measured exponent must clear both.

So "capable of reaching rho parity" is one condition: an oracle on a
45-dimensional factor base that answers a decomposition query,
rejections included, in `2^{17}`–`2^{27}` field products,
`α ≈ 0.4–0.6`. An oracle with `α ≥ 1` cannot get within `2^{21}` of
rho at any `l`. Anything that lowers the constant in front of `2^{αl}`
without lowering `α` is engineering on a flat ratio.

## 3. X7 — the `F_2` SAT oracle's exponent, protocol frozen before the run

**What is measured.** WDSat (Trimoska–Ionica–Dequen, CP 2020), the
published solver for this exact Weil-descended symmetrised `S₄`
encoding, pinned at upstream commit `61c6ff3f` through the frozen
suite's manifest and built with widened static limits so `l ≤ 11`
fits (diff recorded in the receipt; the pinned build answers the
pinned `n19l6` instance identically: UNSAT, 44,050 conflicts). Flags
`-n -l -m 3 -b`, the frozen suite's `symmetry_on` configuration, which
is the fastest of the four (`-x` and `-x -b` are 2.9–17 s against
0.39 s on the pinned instance and still walk tens of thousands of
conflicts).

**Ladder.** `(n, l) = (15,5), (19,6), (21,7), (24,8), (27,9), (30,10)`,
`n ≈ 3l` as in `RESEARCH_SAT_SEMAEV.md`, so a uniform target
decomposes with probability near `1/3!` and both classes occur. Six
uniform targets per rung from seed `0x5A70001`, each an abscissa that
lifts to `K_0`. Irreducibles from `find_irreducible_sparse(n)`,
recorded. Timeout 3,600 s per solver process; a rung with a timeout
is reported and excluded from the fit.

**Reference on the same targets.** The in-repo pairs-and-solve oracle
`semaev_decomp::decompose`: loop over pairs, solve the quartic for the
third abscissa, keep subspace roots by `gcd(q, L_V mod q)`;
`O(2^{2l} · l)` field operations, `α = 2` by construction. It is the
`S₄` decision procedure for the same `b = 1` curve (the Semaev
polynomial does not see `a`).

**Unit.** Per-call cost, three ways: WDSat conflicts (a count, no
hardware); wall seconds; and wall seconds divided by the measured
nanoseconds per `F_{2^n}` multiplication on the same host in the same
run, giving a product-equivalent so the row can sit under §2's
columns. The fit uses rejection instances only. Satisfiable instances
are reported beside them, not fitted.

**Correctness.** Every WDSat `SAT` must carry an assignment whose three
abscissae zero `f₃`; every WDSat `UNSAT` must agree with the exhaustive
pairs-and-solve `None`; every pairs-and-solve witness must zero `f₃`.
A rung with a timeout or a disagreement is excluded from the fit and
reported.

**Fit.** Least squares of `log₂(median cost)` on `l` over complete
rungs, at least four. Extrapolated `log₂ Q` at `l = 45` reported as
an extrapolation.

**Prediction, written before the run.** At `l = 5, 6` WDSat spends
0.92–0.97 conflicts per sorted triple (`2^{3l}/6`), so the predicted
conflict exponent is **`α ≈ 3`**, and the wall exponent is above 3
because instance size grows with `l`; pairs-and-solve is predicted at
`α ≈ 2.1–2.3` in wall (the `l` factor). WDSat is predicted to remain
slower than pairs-and-solve on every rung, by a factor that grows.

**Success and abandonment.** The SAT route would be worth continuing
as an exponent lever if the WDSat exponent came in at `α ≤ 1.5` on
four rungs with every answer verified. It is closed as an exponent
lever if `α ≥ 1.9`. It reaches parity only at `α ≈ 0.38`; no outcome
of this run can establish that, and none is claimed.

**Inadmissible.** Reporting WDSat seconds as `S`; dropping rejections
from the fit; fitting satisfiable instances; changing flags after
seeing a rung; quoting the `l = 45` extrapolation as a measurement.

```bash
cargo run --release --example wdsat_oracle_alpha -- \
  --solver /tmp/wdsat-wide/wdsat_solver --max-l 10 --targets 6 \
  --out experiments/ecc2k130_wdsat_alpha_20260919
```

## 4. X7, measured

Receipt
[`experiments/ecc2k130_wdsat_alpha_20260919/`](experiments/ecc2k130_wdsat_alpha_20260919/):
`summary.json` (every row, every witness, the fit), `run.log`, the
thirty ANF instances, the WDSat `config.h` and its diff against the
pinned source, the build log, and `run-attempt-l10-stopped.log`. Host
`ip-172-31-19-103` (8 cores), code `1dee5d1e`, solver sha256
`bbe7c019…`. Six solver processes ran concurrently; the
pairs-and-solve arm ran alone. That contention is why WDSat's `l = 5`
wall varies 0.041–0.061 s across identical conflict counts, and it
is the reason the conflict column, not the wall column, carries the
headline. WDSat also warns that its widened `__MAX_ID__` makes
"running times not optimal"; on the pinned `n19l6` instance the
pinned and widened binaries take 0.38 s and 0.39 s.

**The `l = 10` rung was not completed.** At the fitted rate a
rejection costs about `180 × 10^6` conflicts and 75–90 minutes, past
the 3,600 s timeout the protocol fixed; the first attempt was stopped
at `l = 9` after all five lower rungs had finished, and rerun to
`l = 9` for a clean receipt. Five complete rungs exceed the protocol's
minimum of four. The stopped log is retained.

### 4.1 One table, one unit

Per-call cost of a **rejection** (median over the 5–6 rejection
targets per rung). Products are wall seconds divided by the measured
`F_{2^n}` multiplication time on this host (5.8–6.7 ns, in the
receipt), so WDSat and pairs-and-solve sit in one unit; a conflict is
not a product and is listed as a count.

| `n` | `l` | `\|F\|` | sorted triples | rejections / hits | WDSat conflicts | per triple | WDSat products | pairs-and-solve products | WDSat / pairs | verified |
|---:|---:|---:|---:|:--:|---:|---:|---:|---:|---:|---|
| 15 | 5 | 32 | 5,984 | 6 / 0 | 5,510 | 0.921 | `2^{23.33}` | `2^{15.58}` | 215× | 6/6 |
| 19 | 6 | 64 | 45,760 | 5 / 1 | 43,953 | 0.961 | `2^{26.39}` | `2^{17.70}` | 414× | 6/6 |
| 21 | 7 | 128 | 357,760 | 5 / 1 | 350,976 | 0.981 | `2^{29.61}` | `2^{19.82}` | 882× | 6/6 |
| 24 | 8 | 256 | 2,829,056 | 5 / 1 | 2,806,146 | 0.992 | `2^{33.29}` | `2^{21.95}` | 2,593× | 6/6 |
| 27 | 9 | 512 | 22,500,864 | 5 / 1 | 22,400,825 | 0.996 | `2^{36.41}` | `2^{24.09}` | 5,112× | 6/6 |

Every rung's satisfiable instance was found by both arms and its
witness zeroes `f₃`; hits cost WDSat 0.19–0.92 of the rejection
conflicts (trajectory luck) and are not fitted. Pairs-and-solve costs
15–20 products per `(pair, l)`, i.e. its constant is the quartic
solve, and WDSat costs 1,900–4,100 product-equivalents per conflict,
rising with `l`.

### 4.2 The fits

Least squares of `log₂(median)` on `l`, five rungs:

| Arm | unit | `α` | `log₂` intercept | `log₂ Q` at `l = 45` *(extrapolation)* | `G_IC / G_rho` at `l = 45` *(extrapolation, §2 frame)* |
|---|---|---:|---:|---:|---:|
| Pair enumeration *(identity, `2^{2l} · 10`)* | products | 2 | 3.32 | 93.3 | `2^{+66.65}` |
| Pairs-and-solve *(measured)* | products | **2.128** | 4.93 | 100.7 | `2^{+74.0}` |
| WDSat `-b` *(measured)* | conflicts | **2.998** | −2.56 | 132.3 | `≥ 2^{+105.7}` if a conflict were one product |
| WDSat `-b` *(measured)* | products | **3.307** | 6.66 | 155.5 | `2^{+128.8}` |
| Stored pair table *(identity)* | products | 1 | 3.32 | 48.3 | `2^{+21.65}` |
| Parity *(condition)* | products | 0.59 | 0 | 26.67 | `2^{0}` |

The conflict exponent is `3.00` to three figures with the per-triple
count converging on 1 from below: WDSat with symmetry breaking visits
every sorted triple once and refutes it, which is `2^{3l}/6` by
construction. That is the same finding the in-repo CDCL reached at
`l = 5..8` (`RESEARCH_SAT_SEMAEV.md`, `conflicts/triple` 0.97–1.02)
and the G7e CryptoMiniSat diagnostic reached at `n = 9` (228× slower
than a table lookup), now on the published solver for this encoding,
on a fifth rung, verified on every instance.

Footnote, not a fit (`xg_one_point_check.txt`): WDSat's
Gaussian-elimination mode `-x -b` spends 0.655, 0.761, 0.840 conflicts
per sorted triple at `l = 6, 7, 8`, rising toward 1, at 15–17× the
wall; two-point slope 3.15 in conflicts. It does not prune either.

### 4.3 What the measurement says

- **Prediction confirmed, abandonment condition met.** Predicted
  `α ≈ 3`; measured `2.998`. The falsifier for continuing the SAT
  route as an exponent lever was `α ≥ 1.9`.
- **Against the reference it is meant to replace, it loses on every
  rung and the gap widens**, 215× to 5,112×, because the exponent is
  higher (`3.0` versus `2.13`) and the constant per leaf is larger
  (a conflict is ~2,000–4,000 products; a quartic solve is ~15–20).
- **Against the parity condition it is not close in kind.** Parity
  wants `α ≈ 0.4–0.6`; a stored table already has `α = 1` and is
  `2^{21.65}×` rho at best. An oracle at `α = 3` is worse than not
  having an oracle: streaming pair enumeration at `α = 2` is
  `2^{39}` cheaper at `l = 45` even if a conflict were one product.
- **Class: measurement.** The ratio to the floor did not move; no
  constant in this pipeline changed. Nothing here is an advance, and
  nothing is engineering either, since the SAT arm was never the
  cheaper arm at any `l`.

### 4.4 What would have to be true

The three `F_2` SAT solvers this repository has now measured on the
Weil-descended Semaev systems — the in-repo CDCL, CryptoMiniSat,
WDSat with and without XOR-Gauss — all spend one conflict per
candidate tuple. The reason is structural rather than a matter of
solver tuning, and `RESEARCH_SAT_SEMAEV.md` named it: the descended
`S₄` constraints bind only once all three abscissae are nearly
determined, so no partial assignment is refutable early. An oracle
with `α < 2` needs a model that constrains *partial* tuples. The
candidates that survive in this repository's ledger are algebraic,
not SAT: a first-fall-degree bound that stays flat as `n` grows (X5's
`m = 4` chained symmetrised FFD is flat at 4 through `n = 15`, which
is evidence in that direction but not a solve), or a Crossbred
`(D, k)` rule with a determining space past the toy rungs (X1–X3
have none). Both are open; neither has an `α`; and even a linear
oracle leaves `2^{21}×` rho. That is where the question "can this
pipeline reach rho parity" stands: closed for SAT, open in theory for
Gröbner-type oracles, with no measured route to the required
`α ≈ 0.4–0.6`.

## 5. Scoreboard

Panel `#ecc2k130-best-20260919` on
[`docs/index-calculus-scoreboard.html`](docs/index-calculus-scoreboard.html):
the §1 best-algorithm rows, the §2 parity scale, and the §4.2 fits,
with the `l = 45` column marked as extrapolation. The page's verdict
sentence is unchanged by this round: the `n = 131` product-law floor
is untouched and every SAT row sits above pair enumeration.
