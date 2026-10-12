# lit2026 W10: can Coppersmith be the prime-field `S₅` small-root oracle?

**Date.** 2026-10-11. **Workstream.** W10 of [`../lit2026_speedup_program/`](../lit2026_speedup_program/README.md)
(kind: feasibility). **Branch.** `lit2026/w10-s5-coppersmith`, stacked on the program commit
`3c239129d` (PR aburan28/crypto#1681). **Code.** [`examples/lit2026_w10_coppersmith_bound.rs`](../../examples/lit2026_w10_coppersmith_bound.rs)
(Rust; Python only plots). **Run.** [`reports/run-20261011T2118Z/`](reports/run-20261011T2118Z/)
(`report.pdf`, `run.md`, `run.json`). **Results PDF.** [`W10_results.pdf`](W10_results.pdf).

**Verdict.** For the strategy family defined in §3, the asymptotic reach is
`δ*(S₃) = 1/6`, `δ*(S₄) = 1/24`, `δ*(S₅) = 1/80`, and no lattice in the family can exceed
`1/(m·d)` = `1/4`, `1/12`, `1/32` at any multiplicity. Arity `m = 4` needs `δ ≈ 1/5`.
`δ*(S₅) < 1/5` by a factor of 16, and the family-wide ceiling `1/32` is 6.4× short. **The
Coppersmith route to the `m = 4` oracle is closed for full-rank triangular single-polynomial
lattices over `(p, f)^ℓ`:** Jochemsz–May basic and extended shifts, any monomial set, any
leading vertex, in both the `x`-encoding and the elementary-symmetric `e`-encoding. Some routes
are not covered (§7), so this is not a proof that no small-root method exists.

## 1. Question and sources

Prime-field EC index calculus first beats rho at arity `m = 4`. That arity needs an oracle that
finds roots of Semaev's `S₅(x₁, …, x₄, x_R) ≡ 0 (mod p)` with `|x_i| < p^δ`, `δ ≈ 1/5`
(`RESEARCH_PRIME_FIELD_BREAKTHROUGH_PROGRAM.md` §2.4; §2.1's relation-count bound gives `1/6`;
`RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md` §0, §7–§8). W10 asks whether the 2024–26
automated-Coppersmith tools move the lattice reach to that range.

Papers (`research/literature_search_20261011/papers.csv` rows 35–37):

- Ryan, ePrint 2024/1577 (EUROCRYPT 2025, "cuso"). Gröbner-basis-optimal shift polynomials for a
  given monomial set, a graph heuristic for the lattice rank, and polytope-based automatic
  asymptotic bounds.
- Ding–Dai–Wu–Zhang–Zhang, ePrint 2026/1027. Asymptotic bounds for the Meers–Nowakowski
  automated Coppersmith method, posed as a linear program.
- Ding et al., ePrint 2026/423. Determinant-based elimination for generating shifts (MIHNP, ICG).

**Access limit.** `eprint.iacr.org` returned HTTP 403 to the PDF fetch of 2026/1027 (bot
protection). We did not try to bypass it, and no arXiv or author-hosted copy turned up in a
search. Only the three abstract pages were read. The LP below is therefore **our own**
formulation, built from the abstracts' description: exact optimisation over monomial sets for the
Meers–Nowakowski/cuso shift structure. It is not a reimplementation of 2026/1027's LP, and we
could not check our numbers against that paper's tables. The known bounds we reproduce instead
are listed in §5.1.

**The "measured 1/80" record.** The repo's `1/80` comes from `RESEARCH_PRIME_FIELD_BREAKTHROUGH_PROGRAM.md`
§2.2. It is the *derived* Jochemsz–May basic-lattice limit `2/(m(m+1)d)` for a box support, not a
measurement. Only `S₃` was measured (probes B/C, 40- and 56-bit: Howgrave-Graham holds to
`δ = 0.10`, the heuristic criterion fades at 0.12–0.16). `RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md`
lines 55–56 and 384–385 call it a "measured Coppersmith reach". This study confirms the
derivation (JM basic at `ℓ = 10⁴`: 0.01250). It also shows the value is the optimum of the whole
box family, not just of JM basic (§5.3).

## 2. The small-root problem, exactly

Data: a prime `p`; `E: y² = x³ + ax + b` over `F_p`; a target x-coordinate `z = x(R)`.
Unknowns: `x₁, …, x_m` (`m = 2, 3, 4`), with `|x_i| < X = p^δ`.
Equation: `f(x₁, …, x_m) := S_{m+1}(x₁, …, x_m, z) ≡ 0 (mod p)`, where
`S₃` is Semaev's closed form, `S₄ = Res_Y(S₃(x₁,x₂,Y), S₃(x₃,x₄,Y))` and
`S₅ = Res_Y(S₄(x₁,x₂,x₃,Y), S₃(x₄,x₅,Y))`.

**Monomial support** (computed exactly by tensor interpolation over a 61-bit prime field with
random `a, b, z`; run §2):

| poly | m | degree per variable | support | Newton polytope | total degree | `e`-encoding monomials (weighted degree) |
|:--|---:|---:|---:|:--|---:|---:|
| `S₃` | 2 | 2 | 9 of 9 | box `[0,2]²` | 4 | 6 (4) |
| `S₄` | 3 | 4 | 125 of 125 | box `[0,4]³` | 12 | 35 (12) |
| `S₅` | 4 | 8 | 6561 of 6561 | box `[0,8]⁴` | 32 | 495 (32) |

After substitution every polynomial is **dense on its box** and symmetric in `x₁..x_m`. Checks:

- the interpolant equals the direct resultant at 8/8 off-grid points, so the degree bounds hold;
- `S_{m+1}` vanishes at 8/8 random genuine relations `P₁ + … + P_m = R`;
- the interpolant vanishes at 4/4 planted decompositions of `z`;
- the symmetric reduction to `G(e₁, …, e_m)` ends with a zero remainder, which is the exact
  re-expansion check.

`S₃` is the sanity case, where the answer is known (`1/6`, and repo probes B/C). `S₄` is the
calibration case.

## 3. Strategy family and the exact optimiser

**Family `F`** (full-rank, triangular, single polynomial, known modulus). The family has four
ingredients:

- a vertex `a` of `Newt(f)` as leading monomial;
- a multiplicity `ℓ`;
- a finite monomial set `M`;
- for each `β ∈ M`, the shift `p^{ℓ−k} x^{β−ka} f^k` with `k = k(β) ∈ [0, ℓ]`, whose whole support
  lies in `M`.

The basis is triangular in any order that `a` maximises over `Newt(f)`. As `p → ∞` at fixed `ℓ`,
Howgrave-Graham holds iff

  `Σ_{β∈M} k(β) > Σ_{β∈M} λ·β`,  `λ_i = log_p X_i`.

`F` contains Jochemsz–May basic (`M = Newt(f^ℓ)`) and extended strategies. It also contains every
monomial-set choice the Meers–Nowakowski automation or cuso's graph heuristic could make for one
polynomial over `(p, f)^ℓ`.

**Exact LP (closure).** Fix `a`, `ℓ` and a bounding box. Define node `(β, j)` to mean "the
level-`j` shift at `β` has its support in `M`". Then choosing `(M, k)` is a maximum-weight-closure
problem:

- a unit gain per level `j ≥ 1`;
- a cost `λ·β` per selected monomial;
- implications `(β, j) → (β + s − a, j−1)` for `s ∈ supp f`.

For a box support the implications become 2^m sub-box edges, and a unit test checks that they
match the sumset form. The constraint matrix is totally unimodular, so the LP optimum is integral.
We solve it exactly by min-cut (Dinic) inside Dinkelbach iteration on the ratio, with
`δ_ℓ = max_M (Σk − Σe·β)/Σc·β` returned as an exact rational. Under the permutation symmetry the
maximal optimal closure is invariant, so we solve on orbits.

**Ceiling for the whole family, every `ℓ` (charging argument).** Take a vertex `c` with
`(k−j)a + jc ∈ supp(f^k)` for `1 ≤ j ≤ k ≤ ℓ`. This holds for `c = a` always, and for other `c`
whenever no coefficient cancels in `f^k`. That is generic, and the closure model assumes it by
taking full sumsets.

- The shift at `β` puts `u_j = β + j(c−a)` in `M` for `j = 1..k(β)`. Charge one unit of `k(β)` to
  each `u_j`.
- A receiver `u` gets at most one unit per `j`, and every sender satisfies `u ≥ jc`. So
  `charge(u) ≤ J(u) := min_{c_i>0} ⌊u_i/c_i⌋ ≤ λ·u / λ·c`.
- Summing, `Σk ≤ Σ_M λ·u / λ·c`. The success condition therefore forces **`λ·c < 1`**: every
  monomial of `f` must be smaller than `p` at the bounds.

For `S_{m+1}` with equal bounds and `c = (d, …, d)`, that gives `δ < 1/(md)`: `1/4`, `1/12`, `1/32`.
In the `e`-encoding the weighted degree of `G` is the same `md`, so the ceiling is the same. The
engine checks `Σk ≤ Σ J(u)` on every optimum it returns (it held with equality in all 70 cells).
This is an inequality, not a polynomial identity. The identity-certificate format
(`docs/identity-certificates/README.md`) certifies polynomial equalities and has no slot for it,
so no certificate is shipped. The proof is the paragraph above.

**Box family in closed form.** `M = [0, L]^m`, `L = λ·dℓ`, any vertex. As `ℓ → ∞`,
`δ(λ) = 2[1 − (1 − 1/λ)^{m+1}] / (m(m+1)d)` for `λ ≥ 1` and `2/(m(m+1)d)` for `λ ≤ 1`. The
supremum over `λ` is therefore exactly `2/(m(m+1)d)`. The best finite-`ℓ` extension `t = L − dℓ`
only speeds up convergence (part `families`, `ℓ ≤ 128`).

## 4. Commands

```bash
export CARGO_TARGET_DIR=/Volumes/SSD990-2/crypto/targets/lit2026-w10
cargo test  --release --example lit2026_w10_coppersmith_bound          # 6 tests
cargo build --release --example lit2026_w10_coppersmith_bound
./lit2026_w10_coppersmith_bound --part all --max-ell 3 --json run.json > run.md   # wrapper: reports/run-*/run_w10.sh
PY=<venv>/bin/python
$PY -I plot_w10.py reports/run-20261011T2118Z/run.json reports/run-20261011T2118Z <label>
$PY -I -c "import typst; typst.compile('reports/run-20261011T2118Z/report.typ', output='reports/run-20261011T2118Z/report.pdf')"
$PY -I -c "import typst; typst.compile('W10_results.typ', output='W10_results.pdf')"
```

Engine commit `a21b5de50`. The run started 2026-10-11 21:18 UTC and took 2564 s wall (301 s
user) on the shared M4 Pro at load 108–177. It is deterministic, so no isolation is claimed and
the `ms` columns are informational.

## 5. Results

### 5.1 Known bounds reproduced (validation)

| case | this engine | literature |
|:--|:--|:--|
| univariate degree `d ∈ {2,3,5}`, `ℓ ≤ 16` | closure optimum = `ℓ/(dℓ+d−1)` in 30/30 cells; → `1/d` | Coppersmith `1/d`; classic finite-`ℓ` lattice |
| Boneh–Durfee, `X = e^δ`, `Y = e^{1/2}`, `ℓ ≤ 32` | closure optimum = BD x+y-shift lattice with best `t`, every `ℓ` (0.27949 at `ℓ = 32`) | — |
| Boneh–Durfee x+y shifts, `ℓ = 2000`, `τ = 0.215` | 0.28466 | `(7−2√7)/6 = 0.28475` |
| JM basic, box `[0,2]²`, `ℓ = 1, 2, 3` | 1/18, 1/10, 0.11905 | repo probe-B theory column 1/18, 1/10, 0.119 |
| JM basic box limits `m = 2, 3, 4` | 0.16665, 0.04166, 0.01250 at `ℓ = 10⁴` | `2/(m(m+1)d)` = 1/6, 1/24, 1/80 |

Boneh–Durfee's 0.2929 (`1 − 1/√2`) is outside `F`. It needs a non-full-rank sublattice, or
Herrmann–May's unravelled linearisation `u = 1 + xy`, and the full-rank optimum here stops at the
0.2847 lattice, as expected. The closure optimiser is exact (min-cut), so the univariate and BD
rows show it finds the known optimal shift sets on its own.

### 5.2 Reach for `S₃`, `S₄`, `S₅`

| | `S₃` (m=2) | `S₄` (m=3) | `S₅` (m=4) |
|:--|---:|---:|---:|
| `δ*` in `F` (box limit, `ℓ → ∞`) | **1/6 = 0.1667** | **1/24 = 0.0417** | **1/80 = 0.0125** |
| ceiling of `F`, every `ℓ`: `1/(md)` | 1/4 | 1/12 | 1/32 = 0.03125 |
| exact closure optimum, largest `ℓ` run | 0.154789 (`ℓ = 12`) | 0.032955 (`ℓ = 5`) | 0.006334 (`ℓ = 3`, `L = 32`) |
| `e`-encoding closure optimum | 0.136111 (`ℓ = 3`) | 0.024903 (`ℓ = 2`) | 0.004061 (`ℓ = 1`) |
| `δ` needed (repo §2.4) | 1/2 (never wins) | 1/4 (ties) | **1/5** (first win) |
| repo measurement | 0.10 (HG) / 0.12–0.16 (heuristic) | none | none |

### 5.3 The closure optimum is always a box

The exact optimum over **all** monomial sets in the bounding box equals the best box
`[0, L′]^m` (`L′ ≤ L`) in **70 of 70** cells. That is every `(poly, ℓ, L)` run and every
leading-vertex class (`a = d·1_S`, `|S| = 0..m`), and all classes give the same rational. The cells:

- `S₃`: `ℓ = 1..12`, all 3 vertex classes; the box never binds.
- `S₄`: `ℓ = 1..5`, all 4 classes, `L = 4ℓ + t_opt + 4`; the box never binds.
- `S₅`:
  - `ℓ = 1`: corner and constant vertices at `L = 44`, where the box does not bind and the optimum
    is `L′ = 38`, `δ₁ = 0.005253`, dimension 2 313 441; all 5 classes at `L = 24`.
  - `ℓ = 2`: corner and constant vertices at `L = 32`; mixed vertices at `L = 24`.
  - `ℓ = 3`: corner vertex at `L = 32`.

  The `ℓ ≥ 2` boxes bind, so those values are lower bounds on `δ_ℓ`.

Up to these multiplicities, then, the 2026-style automation (all monomial sets, every leading
monomial) adds nothing over Jochemsz–May extended boxes for Semaev polynomials. Since the box
family's supremum is exactly `2/(m(m+1)d)` (§3), the extrapolated reach is `1/80` for `S₅`. The
charging argument caps it at `1/32` whatever happens at higher `ℓ`.

### 5.4 Finite lattices (repo LLL, 32-bit curve `a3-fast-32bit-p256class`)

`p = 4294967291`, `a = −3`, `b = 31`, prime order `n = 4294959973`. Planted `x(P_i) ≤ B = ⌊p^δ⌋`.
Success means at least `m` LLL-reduced rows vanish over Z at the root.

| poly | ℓ | dim | `≥ m` rows vanish | Howgrave-Graham on first `m` rows |
|:--|---:|---:|:--|:--|
| `S₃` | 1 | 9 | all trials to δ = 0.10 (B = 9); 4/8 at 0.12; 3/8 at 0.14; 0/7 at 0.16 | to 0.06 (6/7); 3/6 at 0.08–0.10; 0 from 0.12 |
| `S₃` | 2 | 25 | all trials to δ = 0.12 (B = 14); 6/8 at 0.14; 4/8 at 0.16; 0/8 at 0.20 | all to 0.10; 1/7 at 0.12; 0 from 0.14 |
| `S₄` | 1 | 125 | **0/36** at δ = 0.035–0.125 (B = 2 … 15) | 0/36 |

`S₃` reproduces the repo's probe-B picture at 32 bits. `S₄` at `ℓ = 1`, the JM basic lattice of
dimension 125, recovered **no** planted root at any tested bound, including the smallest, `B = 2`.
That matches its `p → ∞` reach `δ₁ = 0.00133`, which means `X < 1.03` at 32 bits. The next level,
`ℓ = 2`, has dimension 729 and reach `δ₂ = 0.0144`, i.e. `X < 2^{0.46}` at 32 bits. Even that
lattice predicts no nontrivial root at this size, and at this dimension the repo's BigInt LLL is
too slow to run (about 1 s per dimension-125 reduction here). It was not run.

## 6. Verdict (task item 4)

`δ*(S₅) = 1/80 < 1/5`, and the family ceiling is `1/32 < 1/5`. **Closed:** the Coppersmith route
to the `m = 4` oracle, for full-rank triangular lattices generated by `p^{ℓ−k} x^γ f^k` from the
single polynomial `f = S₅(·, z)`. This covers:

- Jochemsz–May basic and extended shifts;
- every monomial set: exact optimum for `ℓ ≤ 3` inside the bounding boxes, ceiling for all `ℓ`;
- every vertex of `Newt(f)` as leading monomial;
- the `x`-encoding and the elementary-symmetric `e`-encoding with bounds `X^i`.

There is no lattice dimension to report, because no strategy in the family reaches `1/5`.

## 7. Limits: routes this does not close

- **Non-full-rank and sublattice strategies.** These include Boneh–Durfee geometrically
  progressive sublattices, Herrmann–May unravelled linearisation and other relation-based
  rewritings. The `e`-encoding is one such rewriting and does no better. We have no general
  ceiling for this class.
- **cuso's Gröbner-optimal shifts for a fixed `M`.** In principle they can beat the obvious shift
  of `β` by cancelling out-of-`M` monomials. `F` and the charging argument assume the shift's own
  support lies in `M`. Our closure covers that case only.
- **More than one polynomial.** Examples: adding the curve equations `y_i² = x_i³ + ax_i + b`, or
  sharing unknowns across several targets. Multi-polynomial Coppersmith over the same `p` was not
  analysed.
- **Unequal or non-box bounds.** Examples: factor bases that are not intervals, or lopsided
  bounds per variable.
- **Finite-`p` effects.** Heuristic successes above the determinant bound (seen for `S₃`) do not
  reach `1/5` there either. For `S₄` and `S₅` they were not measurable at a tractable dimension.
- **The 2026/1027 LP itself.** It could not be read (HTTP 403), so it was not reproduced.

## 8. Graphs, canonical graphs, journey

- [`reports/run-20261011T2118Z/delta_vs_arity.svg`](reports/run-20261011T2118Z/delta_vs_arity.svg):
  δ against arity, with the needed-δ line and the repo's `S₃` measurement.
- [`delta_vs_multiplicity.svg`](reports/run-20261011T2118Z/delta_vs_multiplicity.svg):
  convergence in `ℓ`.
- [`finite_lattice.svg`](reports/run-20261011T2118Z/finite_lattice.svg): the 32-bit lattices.
- Graph source: `plot_w10.py`; data: `run.json`.
- Canonical graphs checked and left unchanged, because this is a feasibility negative with no
  speed or IC-versus-rho result: `docs/ic/progress-timeline.json`,
  `docs/index-calculus-scoreboard.html`, `docs/performance-gains/` (none mentions the prime-field
  Coppersmith reach).
- Journey: `../lit2026_speedup_program/journey/entries/20261011-W10-*.json` (kinds `negative` and
  `diagnostic`, no arms; a feasibility computation has no original/candidate binary pair).
