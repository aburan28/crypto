# Combining index-calculus ideas on ECC2K-130

**Question.** After the threads in this repository have each priced one
lever — pair enumeration, a stored pair table, SAT, Frobenius, homogeneous
meet-in-the-middle, Weil descent, residual-walk linear algebra, symmetrised
Semaev, Crossbred — which *combinations* survive the counting, oracle, and
linear-algebra constraints at once, and could any of them move the exponent
rather than a constant?

**Answer.** None of the measured combinations is an advance. The only
coherent multi-structure stack on this curve is already the one G7e ran:
type-II ONB × Hamming-weight-2 even-trace × Frobenius-orbit collapse ×
pair enumeration or a stored pair table. GPU search, SAT, and the table
compose with that stack and still sit on the product law. Invariant
subspaces, quasi-subfield polynomials, polynomial-basis Hamming weight,
GHS covers, and extension-field base change do not compose with it —
they need divisor degrees of `t^{131} − 1` that are not in `{1, 130}`.
Homogeneous relations that determine nothing, Wagner on curve addition,
and Joux–Vitse at prime degree are closed, not slow. Five solver
experiments remain compatible with those identities and could still move
the oracle exponent `α`; none of them has a fitted `α` yet, so none is a
result.

This note does not rerun anything. Every number is cited from a frozen
artefact or from a thread that already stated its boundary. The scoreboard
panel is `#ecc2k130-ic-synthesis-20260918`.

Companion map: [`RESEARCH_ECC2K130_ROUTES.md`](RESEARCH_ECC2K130_ROUTES.md)
and [`RESEARCH_ECC2K130_ROUTE_TARGETS.md`](RESEARCH_ECC2K130_ROUTE_TARGETS.md)
(what to try next, as X1–X6). Measured G7e oracles:
[`RESEARCH_ECC2K130_G7E_INDEX_CALCULUS.md`](RESEARCH_ECC2K130_G7E_INDEX_CALCULUS.md).

## 1. The boundary, stated before combining

Per `AGENTS.md` §1 this note has both kinds of boundary, written before the
combination table. The instance is `K_0 : y² + xy = x³ + 1` over `F_2^131`,
cofactor 4, prime subgroup order

```
r = 680564733841876926932320129493409985129,   log2 r = 129.0
#E = 4r ≈ 2^131
```

**The reference** is Pollard rho with the `⟨−1⟩ × ⟨π⟩` speed-up on the same
subgroup: `2^60.809` walk iterations. On the G7e SKU the shipping packed
client charges **5.3125 field products per iteration**
(`ecc2k130/THROUGHPUT-29B.md`). In the product unit that is

```
rho products = 5.3125 × 2^60.809 ≈ 2^63.22
S_rho        = (5.3125 × 2^60.809) / 2^64.5 ≈ 0.411
```

cited from `ecc2k130/benchmarks/indexcalc-g7e/summary.json`.

**The floor** is the product law for an `m`-summand factor base `F`. A
random target decomposes with probability `C(|F|, m) / #G`. Streaming pair
enumeration spends `C(|F|, 2)` affine additions on every trial. Charging
`|F|` relations, at leading order in `|F|`, costs

```
|F| · (#G / C(|F|, 3)) · C(|F|, 2) · 20   field products
    ≈  60 · #G
```

independent of `|F|`. For `G = E(F_2^131)` that is `60 · 2^131 ≈ 2^136.91`
field products, `2^73.69` times the rho product count. The Frobenius
collapse (`|F|/131` orbit unknowns) divides by `131 ≈ 2^7.03` and leaves
`2^66.66` times rho. Affine addition is priced at 10 field products, as in
the G7e note. Squares and Frobenius permutations are not charged.

**The oracle inequality** that ranks every algebraic combination, restated
from the route-target note before any Crossbred number is named. Substituting
a per-call oracle cost `Q` for the `C(|F|, m−1)` term:

```
Λ  =  m! · Q / 2^{(m−1)l}        (m = 3:  Λ = 6Q / 2^{2l})
```

`Λ` falls below `m` exactly when `Q < C(|F|, m−1)`. An oracle moves the
ratio exponentially in `l` only if `Q = 2^{α l}` with **`α < 2`**. A
constant-factor speed-up, however large, leaves `Λ` flat and is
**engineering** by `AGENTS.md` §3.

At `n = 131`, `m = 3`, `l = ⌈(n + log₂ m!)/m⌉ = 45`, even an oracle
*linear* in `|F|` — already a major theoretical result — leaves
`2^27.78×` rho in the route-target's group-operation unit. Matching rho
would need `Q ≈ 2^{17.2}`, i.e. `α ≈ 0.38` at `l = 45`. That scale is
cited, not re-derived; it is the numeric condition this note uses for
"advance".

**The linear-algebra caution**, from a different curve and therefore not a
row of the ECC2K-130 table: on `E(F_{p³})` residual-walk relations grow as
`n^{1/3}` and Wiedemann as `n^{0.68}`, so `S` bottoms near `200×–265×` rho
around `2^{50}` and the method's crossover with rho sits past `2^{230}`
([`RESEARCH_RESIDUAL_WALKS.md`](RESEARCH_RESIDUAL_WALKS.md) §11.7). An
exponent claim that prices only the oracle is the §5 mistake, drawn.

**Falsification target, stated in advance.** Combining ideas would be a
success if some compatible stack produced a fitted oracle exponent
`α ≤ 1.5` over four or more rungs, with every decomposition cross-checked
against an independent oracle, every completed discrete log satisfying
`[k]P = Q`, linear algebra priced, and total field products below
`2^63.22`. It is abandoned if every surviving combination has `α ≥ 1.9`
or is structurally inapplicable at prime `n = 131`.

Inadmissible: quoting GPU pairs/second or SAT seconds as `S`; dropping
linear algebra; counting homogeneous relations as a logarithm; moving the
boundary by changing `m`, `n`, or the product-per-add conversion after the
fact; presenting the conditional first-fall-degree formula at `n = 131`
as a measurement.

## 2. Method families, and what they are allowed to compose with

Each family is one idea this repository has already built, surveyed, or
closed. "Compose" here means the output of one is a legal input of another
without changing the instance or the unit.

| Family | What it actually does | Moves `α`? | Composes with | Does not compose with |
|---|---|---|---|---|
| Streaming pair enum | Exhaustive 3SUM, `Q = C(\|F\|, 2)` | No (`α = 2`) | Frobenius orbits, GPU inner loop, planted checks | SAT as a cheaper `Q` at `n = 131` (measured slower) |
| Stored pair table | Pay `C(\|F\|, 2)` once, then `Q = \|F\|` probes | No (`α = 1`, still `2^27.78×` rho if it held at `l = 45`; at the measured `\|F\| = 8384` it is `2^53.62×`) | GPU fill, host sort, Frobenius | A `2^{70}`-below-pairs claim; weight-4 tables that do not fit in 96 GiB |
| SAT / CryptoMiniSat | CDCL on Weil-descended Semaev | Unknown: conflicts are not products | Type-II ONB degrees `{5,9,11,23}`, toy DLP | The product table; GPU pair-enum as if SAT were faster |
| F4 / first-fall degree | Gröbner on chained `S₃` | Conditional `2^{86}` at `n = 131` is still `2^{25}` short of rho, and `D_reg = D_ff + o(1)` is a stalemate | Symmetrised systems (smaller, higher degree) | A rigorous `D_reg` bound (Caminata–Ceria–Gorla is vacuous at 263) |
| Crossbred | Algebra to degree `D`, enumerate `k` variables | **Open.** Premise: cost tracks the *system*, not `\|F\|` | Existing `F2BoolPoly` oracles, GPU search phase, symmetrised frame | Factor-base constructions that need a subfield (`131` prime does not touch it) |
| Symmetrised `u`-frame | Quotient by 2-torsion; ~350× SAT at `n = 15` | Predicted **engineering** (constant) until an e2e `Λ · n / m` slope is fitted | Crossbred, Frobenius view, chained `S₃` at `m = 4` | A claim that 350× is an advance |
| Homogeneous MITM | Relations among base points | Honest logarithm `2^{+6.48}` vs rho *with memory free*; `m = 18` "below rho" is relabelling | Frobenius quotient (memory and unknowns, **not** hit rate) | Wagner `k`-tree (no prefix on elliptic addition) |
| Residual walks | Subspace FB on `E(F_{p³})` | Exponent is real (`n^{4/9}`); constant is hopeless | The bookkeeping rule "price LA" | ECC2K-130 rows (different group, different unit) |
| GHS / Joux–Vitse / extension | Change the field or the cover | Structurally empty or `≥ 2^{70}×` rho | Nothing on this instance | Prime `n = 131`, `ord_{131}(2) = 130` |
| CM / `τ`-adic IC | Koblitz structure beyond Frobenius | **Unsearched** (literature item 5) | Unknown | Rho/`τ`-adic scalar multiplication, which is a different algorithm |

The hard arithmetic facts that pre-filter the closed families, and do **not**
pre-filter Crossbred or CM-for-IC:

- `131` is prime, so there is no intermediate field for a tower.
- `2` is a primitive root mod `131`, so `t^{131} − 1 = (t+1)·Φ_{131}` has
  divisor degrees `{0, 1, 130, 131}` only. Every Frobenius-stable
  `F_2`-subspace, quasi-subfield kernel, and GHS magic number is one of
  those four; there is no intermediate dimension.
- Type-II ONB exists only when `2m+1` is prime and `ord_{2m+1}(2) ∈ {m, 2m}`;
  the SAT ladder is `5, 9, 11, 23`, not `13, 15, 17`.
- The set of `n`-subset sums of a `σ`-stable support is itself `σ`-closed, so
  Frobenius does not multiply hit rate
  ([`RESEARCH_ECC2K130_RELATION_SWEEPS.md`](RESEARCH_ECC2K130_RELATION_SWEEPS.md)
  §5.4).
- The 262 automorphisms `⟨−1⟩ × ⟨π⟩` are already inside the rho reference.
  Quotienting the factor base by `σ` again does not buy a second `√262`.

### The factor-base stack that actually exists

A Frobenius-stable subspace of `F_2^{131}` is a binary cyclic code of
length 131. The 2-cyclotomic cosets are `{1}` and `{130}`, so the
achievable dimensions are `{0, 1, 130, 131}`: `E(F_2)` (4 points) or
about `2^{130}` abscissae. Quasi-subfield polynomials reduce to the same
divisor criterion; Euler–Petit put the linearized class at `β ≥ 3/4`
against the `β < 0.103` needed to beat generic `O(2^{n/2})`
([`RESEARCH_QUASI_SUBFIELD.md`](RESEARCH_QUASI_SUBFIELD.md)). Hamming
weight in a type-II ONB is a different predicate: `{x : HW(x) ≤ w}` is
`σ`-stable for every `w` because Frobenius is a coordinate shift, and it
is **not** a low-degree polynomial, so Semaev descent cannot treat it as
a subspace.

That is why the G7e factor base is the live stack rather than a subspace
one. Even trace is the odd-order filter `Tr(x) = 0`; on this curve every
`m ≥ 2` is admissible, unlike `K_1/F_2^7` where odd `m` fails
([`RESEARCH_ECC2K130_DECOMPOSITION.md`](RESEARCH_ECC2K130_DECOMPOSITION.md)
§3). A `ker Tr` base lives in an index-2 subgroup and buys about `2×`
yield — one bit, compatible with the orbit quotient, not an exponent.
Four-point sums over that support land in `H` of index 2, so the usable
translates are `E[4] ∩ H` of size 2, not 4
([`RESEARCH_ECC2K130_RELATION_SWEEPS.md`](RESEARCH_ECC2K130_RELATION_SWEEPS.md)
§4.4).

Polynomial-basis Hamming weight is **not** in this stack. Squaring sends
low weight out of the set, so a poly-basis weight-2 support cannot be
orbit-quotiented (relation-sweeps §4). That incompatibility is easy to
miss and is why the packed ONB is not an optional representation.

C1–C6 below change the *oracle* on top of this stack. They do not reopen
a subspace, a quasi-subfield, or a cover.

## 3. One table, one unit

Unit: **field products**. Conversion from rho iterations: 5.3125
products/iteration, measured on the G7e SKU. `S = products / √r` with
`√r = 2^{64.5}`. Every cost a method incurs belongs in `S`. Rows whose
native unit is not field products are excluded from this table, not
converted by guesswork; they are listed in §4.

| Variant | log2 products | S | vs rho | vs floor | Correctness | Class |
|---|---:|---:|---:|---:|---|---|
| Pollard rho `⟨−1⟩×⟨π⟩` *(reference)* | 63.22 | 0.411 | `2^+0.00` | — | shipping client | baseline |
| Product-law floor, `m = 3`, any `\|F\|` | 136.91 | `2^{72.4}` | `2^+73.69` | `2^+0.00` | derived | floor |
| + Frobenius collapse `\|F\|/131` | 129.87 | `4.78×10^{19}` | `2^+66.66` | `2^+0.00` | derived | accounting |
| Streaming pairs, G7e, `\|F\|=8384` | 129.87 | `4.78×10^{19}` | `2^+66.66` | `2^+0.00` | 8/8 planted; 0 generator hits, yield `3.6×10^{-29}` | engineering |
| Pair-table identity, `\|F\|=8384` | 116.84 | `5.70×10^{15}` | `2^+53.62` | `2^{-13.03}` vs streaming | derived, generic 3SUM with memory | accounting |
| Pair table, G7e, `\|F\|=8384` | 116.84 | `5.70×10^{15}` | `2^+53.62` | `2^{-13.03}` vs streaming | 8/8 planted; 0 generator hits; 35,137,344 sums | engineering |
| Honest MITM logarithm, memory free | 69.70 | 36.7 | `2^+6.48` | — | 6/6 planted logs on analogues; 1,011 EB uncharged | accounting |
| Homogeneous relations, `m = 28` | 59.17 | 0.025 | `2^{-4.05}` | — | relations determine nothing | **relabelling** |

Cited from `ecc2k130/benchmarks/indexcalc-g7e/summary.json` (rho, floor,
streaming, table) and from
[`RESEARCH_ECC2K130_RELATION_SWEEPS.md`](RESEARCH_ECC2K130_RELATION_SWEEPS.md)
§5.6 (MITM and homogeneous). The MITM row was published in rho-step
equivalents as `2^{67.287}` operations; multiplying by the measured
5.3125 products/iteration converts the unit and leaves the ratio to rho
unchanged. Memory is still charged at zero. Homogeneous `S` is
`2^{59.17 − 64.5}` from that converted `log2`, not the `0.018` printed
in §5.6, which does not equal `2^{56.76 − 64.5}` in the iteration unit
either. The ratio `2^{-4.05}` is the column that note's own rho
reference supports.

Reading the table:

- GPU pair enumeration and the stored table are the same algorithm as their
  counting identities. Ratio to those identities is 1. Class: engineering.
- Combining them (table on the GPU) is the `|F|` cut already priced, not a
  new exponent. Still `2^{53.62}×` rho at this `|F|`.
- The friendliest honest logarithm this repository has written down —
  Frobenius-quotiented MITM, table built once, probes and `E[4]` priced —
  is `2^{6.48}×` rho **and asks for a thousand exabytes**. Charging memory
  makes the gap larger.
- The only row below rho is the one that does not recover `log_P(Q)`.

## 4. What does not enter the product table

| Item | Why it is not a row | What it is | Class if forced |
|---|---|---|---|
| SAT growth `m=5,9,11` and toy DLP `n=5,9` | Conflicts have no measured conversion to field products | Diagnostic; pair lookup 228× faster on 16 planted triples at `n=9` | engineering |
| Conditional FFD formula at `n=131` (`2^{86}`) | Literature bit-complexity, gated on `D_reg = D_ff + o(1)` | Extrapolation, `2^{25}` short of rho even if granted | extrapolation |
| Residual-walk `S ≈ 200×–265×` at `2^{50}` | Different curve `E(F_{p³})` | Phase-pricing lesson: LA can dominate | accounting |
| Crossbred `word_ops` | Conversion factor to group ops is **unmeasured** (`AGENTS.md` §2) | X1 frozen 2026-09-18: no `α` fit (two usable rungs); X2 frozen | unmeasured |
| Symmetrised 350× SAT at `n=15` | Wall-clock on one size, not `Λ` | Size gain; predicted engineering | engineering |
| Crossbred kernel frontier | Bit ops, not field products | Space exists through `n=9, m=3` and `n=13, m=2`; `filters=0`; `n=13, m=3` is not a result | measurement |

SAT receipt: `ecc2k130/benchmarks/indexcalc-g7e/sat.json`. Putting any of
these into `S` without a measured conversion is relabelling.

## 5. Combinations that survive the constraints

These are the stacks whose pieces are compatible and whose obstruction is
empirical rather than structural. They change the *oracle* on top of the
ONB Hamming-weight stack in §2; they do not reopen a subspace or a cover.
They are the route-target experiments X1–X6, restated as combinations.
Predicted class is engineering unless a fitted `α ≤ 1.5` is produced.
None has been run as a combination in this note.

### C1 — Crossbred on the chained `m = 3` systems (X1 + X2)

**Why it survives.** The solver is orthogonal to the factor-base
construction. Prime `n` and primitive `2` do not touch it. Half of it is
already in `src/cryptanalysis/crossbred.rs`, correctness-gated against
matrix-F4 at `m = 2` and `m = 3`. The literature survey found no
publication combining Crossbred with binary ECDLP decomposition systems
([`RESEARCH_ECC2K130_IC_LITERATURE.md`](RESEARCH_ECC2K130_IC_LITERATURE.md)
§5).

**X2, frozen.** A determining space exists on every `agree = yes` rung
tried: `m = 3` through `n = 9` (`ℓ = 6`, `v = 27`) and `m = 2` through
`n = 13` (`ℓ = 12`). The X2 falsifier is not met; Route 1 stays open.
Every printed cell has `filters = 0`. The first larger `m = 3` rung,
`n = 13`, `v = 49`, extracts a kernel but fails the correctness gate
(`agree = NO`, `xb/F4 = 1150`). Receipt:
`experiments/ecc2k130_crossbred_kernel_20260918/`. **X1, frozen.**
Receipt `experiments/ecc2k130_crossbred_x1_20260918/`. T4 is applied.
No fit: only two agreeing `m = 3` rungs have `|F| ≥ 3` (`n = 5`,
`Q/C = 83.248`; `n = 9`, `Q/C = 319.222`). `n = 13` and `n = 23` have
no determining space. The two-point sketch slope is not a result.
Chained Route 1 cannot produce `α`. The remaining `α` measurement is
C2 / X3.

**The combination.** Existing Semaev systems → Crossbred `(D, k)` →
`DecompositionStrategy`, with the search phase (`2^k` independent points)
the one piece that maps onto a GPU without a new argument.

**What would make it an advance.** Least-squares slope of
`log₂(Q / C(|F|, m−1))` against `l` over ≥4 rungs ≤ −0.5 (`α ≤ 1.5`),
every call cross-checked against matrix-F4, Macaulay extraction counted
per target, `(D, k)` rule fixed before the run.

**Falsifier.** Slope ≥ −0.1 (`α ≥ 1.9`), or `kernel_dim = 0` past the toy
rungs so there is no Crossbred space to fit.

**Predicted class.** Unknown *a priori*; **blocked in practice** — the
chained ladder cannot supply four usable rungs, so C1 cannot produce `α`.

### C2 — Crossbred on the symmetrised `u`-frame (X3)

**Why it survives.** `build_symmetrised_system` already returns the
`F2BoolPoly` shape `extract_crossbred` consumes. The `u`-frame is smaller
and higher-degree than the chained `x`-system (13 unknowns at Boolean
degree 4 against 30 at degree 3, `n = 15`;
[`RESEARCH_EXOTIC_COORDINATES.md`](RESEARCH_EXOTIC_COORDINATES.md) §8.2).
Crossbred likes few variables and dislikes high degree, so the direction
is uncertain — which is why it is a run, not an argument.

**The combination.** Symmetrised system → Crossbred, same metric as C1.

**Falsifier.** As C1, plus `kernel_dim = 0` throughout (systems too small
for the technique). That outcome is still worth recording.

**Protocol, frozen 2026-09-18 before the run.** T4 selection, paired-oracle
divisor `divisor_for_dimension(n, (n+1).div_ceil(m))`, `|F| = |F_u|`,
axis `ℓ = dim V`, same slope falsifier as C1. Wired as
`cargo run --release --example crossbred_bench -- --sym`.

**X3, frozen.** Receipt `experiments/ecc2k130_crossbred_x3_20260918/`.
No fit: only two usable rungs (`K_1` `n = 7`, `Q/C = 4.653`; `K_1`
`n = 15`, `Q/C = 6.409`). `K_0` is degenerate at small `n`. Larger
rungs with a real `F_u` have no determining space. `filters = 0`.
C2 cannot produce `α` on this divisor convention.

**Predicted class.** Unknown *a priori*; **blocked in practice**.

### C3 — Symmetrised oracle through end-to-end collection (X4)

**Why it survives.** The 350× is a measured SAT/F4 constant at `m = 3` on
`K_0/F_2^{15}`, and `frobenius_view_of_symmetrised` supplies the orbit
structure the attack needs. Nothing about 350× contradicts the product
law; it is a cheaper constant in front of the same `C(|F|, m−1)` if the
oracle remains enumerative, or a different `Q(l)` if it does not.

**The combination.** `DecompositionStrategy::Symmetrised` → E1 ladder with
**every** phase priced: setup, encoding, failed attempts, solving,
lifting, verification, relation-matrix work, scalar recovery.

**Primary metric.** Slope of `log₂(Λ · n / m)` against `n`. A method that
tracks the Frobenius floor exactly already has `log₂(total)`-versus-`n`
slope ≈ 0.95 on the E1 ladder; that slope is not a falsifier. `Λ · n / m`
is 1 on the floor at every rung.

**Falsifier for "advance".** Slope consistent with zero → class
engineering. An advance requires the ratio to *decrease* with `n`.

**Predicted class.** Engineering. The prediction is pre-registered so a
later 350× e2e cannot be reported as progress if the ratio is flat.

### C4 — Chained symmetrised `S₃` at `m = 4` (X5)

**Why it survives.** Conditional theory wants `m ≈ n^{1/3} ≈ 5.1` at
`n = 131`; the harness reaches `m = 3`. Chaining the *symmetrised* `S₃`
gives `4(ℓ−1) + 1 + n` unknowns with bilinear links, named as missing in
the exotic-coordinates note §8.5. First-fall degree on those systems is
H1 of the scaling target, and nobody has a rigorous answer in either
direction ([`RESEARCH_ECC2K130_IC_LITERATURE.md`](RESEARCH_ECC2K130_IC_LITERATURE.md)
§1).

**The combination.** Symmetrised `S₃` chained for `m = 4` → FFD over 16
draws per instance, plus a solve under 64 unknowns on the scaling metric
`m·ℓ + (m−2)·n` at `m = ⌈n/ℓ⌉`.

**Falsifier.** FFD growing with `n`. That closes the route *and* is a
result: it is evidence on the FFD stalemate, which this repository's
small-`n` ladder currently holds at FFD = 3 for chained `m ≥ 3`.

**Predicted class.** Measurement, not an attack. A growing FFD is a
negative that still counts.

**AutoLab.** Local control plane
[`research/ecc2k130_crossbred_autolab_20260918/`](research/ecc2k130_crossbred_autolab_20260918/)
replays the X1/X3 freezes, refuses an `α` claim without four rungs in
one frame, and launches the chained-`x` FFD smoke. The chained
symmetrised `S₃` at `m = 4` is still unbuilt.

### C5 — GPU search phase of Crossbred, only after C1 has an `α`

**Why it survives.** Crossbred's search is `2^k` independent bitwise ANDs
against a precomputed table. That is the one IC inner loop this
repository's CUDA stack is already shaped for. It does **not** survive as
a substitute for C1: a faster inner loop at `α = 2` is the G7e pair-enum
story again.

**Blocked on the X2 freeze.** Every printed cell has `filters = 0`, so
there is currently nothing for that AND to test. Revisit only if a later
`(D, k)` rule produces `filters > 0` and C1 has an `α`.

**The combination.** C1's fitted `(D, k)` rule → GPU search, conversion
factor from word-ops to field products measured on the same SKU.

**Falsifier.** C1 already closed, or `k` growing like `2l`.

**Predicted class.** Engineering, unless C1 itself is an advance.

### C6 — Second literature pass on CM / `τ` for index calculus (X6)

**Why it survives.** Survey items 4 and 5 returned zero surviving claims,
flagged as absence of evidence. Item 5 is Koblitz structure *beyond*
Frobenius used for index calculus rather than rho: `τ`-adic expansion, CM
by `(1 ± √−7)/2`, class group of `Z[τ]`. The structural constraints do
not obviously touch it.

**The combination.** None until a source or a construction exists.
Composing `τ`-adic scalar multiplication with an enumeration oracle is
rho, not IC.

**Falsifier.** A second independent pass also empty → treat as genuinely
unexplored, then decide whether to build rather than whether to read.

**Predicted class.** Survey, not a row.

## 6. Combinations that look like they compose and do not

Recorded so they are not rebuilt.

| Combination | Why it fails | Class of the failure | Source |
|---|---|---|---|
| GPU pairs + SAT as the `n = 131` oracle | Unbounded SAT is 228× slower than pair lookup at `n = 9`; SAT has no product conversion | engineering, not an attack | G7e `sat.json` |
| Pair table as a `2^{70}`-below-pairs oracle | Generic 3SUM with memory; at `\|F\| = 8384` the cut is `\|F\|`, leaving `2^{53.62}×` rho | engineering | G7e table identity |
| GPU + table + Frobenius + SAT | Product of engineering constants on the same floor | engineering | this table, §3 |
| Homogeneous `m = 18` "below rho" + anything | Relations among unknown logs; no `log_P(Q)` | relabelling | relation-sweeps §5.2 |
| Wagner `k`-tree on curve points | No prefix: `x(P+Q)` is not a truncation of `x(P), x(Q)` | closed | relation-sweeps §5.5 |
| Joux–Vitse cover-and-decomposition | Needs a composite-degree tower; `131` is prime | closed | literature §4, routes |
| GHS / hyperelliptic transfer to `F_2` | Genera `1, 2^{129}, 2^{130}` only; the cheap window is `[130, ~300]` | closed | hyperelliptic note |
| Base change to `F_2^{131e}` | Subspaces are `dim ≤ e` (empty of relations) or `dim ≥ 130` (`≥ 2^{70}×` rho) | closed | extension note |
| MOV / Frey–Rück | Embedding degree `> 10^7` | closed | extension §5 |
| Pairs-and-solve as a Koblitz `m = 3` speed-up | `\|F\|² · ℓ` against enumerate's `\|F\|²`; wins only vs a triple loop never used | closed | routes |
| Yokoyama et al. as a lower bound that ends IC | Conditional, prime fields, naive Semaev, carves out Diem and QSP | closed | literature §4 |
| Frobenius ×131 on hit rate | `n`-subset sums are `σ`-closed | accounting | relation-sweeps §5.4 |
| Weight-3 even-trace base as larger than weight-2 | Odd Hamming weight is skipped by the odd-order trace filter | accounting | G7e note |
| Residual-walk `n^{4/9}` transplanted to ECC2K-130 | Different group; on its own curve it still bottoms ~200× rho | closed as a transplant | residual-walks §11.7 |
| Relation-phase crossover quoted as the method | Residual-walk `2^{96}` was one phase; LA is `n^{0.68}` | accounting | `AGENTS.md` §5 |
| Invariant subspace / cyclic code at `n = 131` | Dimensions `{0,1,130,131}` only | closed | quasi-subfield §5b, extension §2 |
| Quasi-subfield polynomial + this curve | Same divisor gap; linearized `β ≥ 3/4` vs need `< 0.103` | closed | `RESEARCH_QUASI_SUBFIELD.md` |
| Polynomial-basis HW-2 + `σ`-orbit quotient | Squaring leaves the weight class | closed | relation-sweeps §4 |
| GHS isogeny walk to a cheaper magic number | Magic `{1,130,131}` is a property of the field, not the curve | closed | hyperelliptic note |
| Construction equations `σ² + σ + 2 = 0` as DLP relations | Real rank, zero information about `log_P(Q)` | relabelling | relation-sweeps §4.7 |
| Counting `⟨−1⟩ × ⟨π⟩` a second time on the IC side | Already inside `S_rho` | accounting | this §2 |
| Hypothetical genus-130 Jacobian of `A/F_2` with `#A(F_2) = r` | Would be `≈ 2^{37}` vs rho, Torelli codimension 8128, no construction | open, not runnable | hyperelliptic §6 |

## 7. What would count as finishing

Any one of:

- **`α` fitted over ≥4 rungs** on C1 or C2, whatever its value. A number
  closes Crossbred either way; X2 is now frozen and the remaining gap is
  a fourth agreeing `m = 3` rung plus `Q / C(|F|, 2)`.
- **C3 classified** by the `Λ · n / m` slope, not by the 350×.
- **H1 falsified at `m = 4`** — a first fall degree that grows.

None of these threatens a deployed curve. The `α ≤ 0.38` row of the
route-target scale table is what that would take, and nothing on §3 is
within `2^{27}` of it. The 500×-or-worse rho verdict on
`docs/index-calculus-scoreboard.html` is unchanged.

## 8. What this note does not settle

- It is a map, not a measurement of C1–C6. Combining ideas on paper does
  not substitute for the word-op conversion factor, the `(D, k)` kernel
  frontier, or an end-to-end ladder with linear algebra priced.
- The MITM row still charges 1,011 exabytes at zero. That is inherited
  honesty from the relation-sweep note, not a new discount.
- Small-characteristic ECDLP remains a literature stalemate on `D_reg`.
  This repository's FFD ladder is evidence, not a proof.
- Prime-field structured factor bases, Nagao function-first solvers, and
  the WDSat regression suite are out of scope: they are scored on their
  own panels and do not become ECC2K-130 rows by juxtaposition.
- A simple abelian variety `A/F_2` of dimension 130 with `#A(F_2) = r`,
  isogenous to a Jacobian, is the one cover that would undercut rho
  (`≈ 2^{37}`). It is not a factor-base combination and has no
  construction in this repository.

## Sources (the threads this map is of)

- [`RESEARCH_ECC2K130_G7E_INDEX_CALCULUS.md`](RESEARCH_ECC2K130_G7E_INDEX_CALCULUS.md)
  — GPU streaming and pair table, SAT diagnostic.
- [`RESEARCH_ECC2K130_DECOMPOSITION.md`](RESEARCH_ECC2K130_DECOMPOSITION.md)
  — product law, existence vs finding.
- [`RESEARCH_ECC2K130_RELATION_SWEEPS.md`](RESEARCH_ECC2K130_RELATION_SWEEPS.md)
  — homogeneous MITM, Frobenius hit-rate, Wagner.
- [`RESEARCH_ECC2K130_EXTENSION.md`](RESEARCH_ECC2K130_EXTENSION.md),
  [`RESEARCH_ECC2K130_HYPERELLIPTIC.md`](RESEARCH_ECC2K130_HYPERELLIPTIC.md),
  [`RESEARCH_QUASI_SUBFIELD.md`](RESEARCH_QUASI_SUBFIELD.md)
  — prime `n`, primitive `2`, cyclic-code dimensions `{0,1,130,131}`.
- [`RESEARCH_ECC2K130_IC_LITERATURE.md`](RESEARCH_ECC2K130_IC_LITERATURE.md),
  [`RESEARCH_ECC2K130_ROUTES.md`](RESEARCH_ECC2K130_ROUTES.md),
  [`RESEARCH_ECC2K130_ROUTE_TARGETS.md`](RESEARCH_ECC2K130_ROUTE_TARGETS.md)
  — FFD stalemate, Crossbred as open ground, X1–X6.
- [`RESEARCH_EXOTIC_COORDINATES.md`](RESEARCH_EXOTIC_COORDINATES.md),
  [`RESEARCH_KOBLITZ_SCALING_TARGET.md`](RESEARCH_KOBLITZ_SCALING_TARGET.md),
  [`RESEARCH_KOBLITZ_INDEX_CALCULUS.md`](RESEARCH_KOBLITZ_INDEX_CALCULUS.md)
  — symmetrised frame, FFD = 3 on the measured chained ladder, which
  invariant bases exist.
- [`RESEARCH_RESIDUAL_WALKS.md`](RESEARCH_RESIDUAL_WALKS.md) §11.7
  — price every phase.
- Frozen receipts:
  `ecc2k130/benchmarks/indexcalc-g7e/summary.json`,
  `ecc2k130/benchmarks/indexcalc-g7e/sat.json`.
