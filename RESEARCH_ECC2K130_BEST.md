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
per-call cost `Q = 2^{αl}` has `α ≈ 0.38` at `l = 45` (§2). The `F_2`
SAT route is measured in §4 on a six-rung ladder against exactly that
scale: WDSat's fitted exponent is **[to be measured]**, against `2` for
pair enumeration and `0.38` for parity. Scoreboard panel
`#ecc2k130-best-20260919`.

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

per-call oracle cost `Q`. Against rho's `2^{63.22}` products:

| Oracle | `Q` | `α` in `Q = 2^{αl}` | best `l` | `G_IC / G_rho` |
|---|---|---:|---:|---:|
| Streaming pair enumeration | `C(2^l, 2) · 20` | 2 | any | `2^{+66.66}` (after `/131`) |
| Stored pair table | `2^l · 10`, memory `2^{2l−1}` | 1 | 45 | `2^{+27.78}`, memory `2^{89}` entries |
| Measured pair table | as above at `l = 13.03` | 1 | — | `2^{+53.62}` |
| **Parity** | `≈ 2^{17.2}` products at `l = 45` | **0.38** | 45 | `2^{0}` |

So "capable of reaching rho parity" is one condition: an oracle at
`l = 45` that answers a decomposition query in about `2^{17}` field
products, on a 45-dimensional factor base, and does so for the
rejections as well as the hits. An oracle with `α ≥ 1` cannot get
within `2^{27}` of rho at any `l`. Anything that lowers the constant in
front of `2^{αl}` without lowering `α` is engineering on a flat ratio.

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
recorded.

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

[to be filled from `experiments/ecc2k130_wdsat_alpha_20260919/summary.json`]
