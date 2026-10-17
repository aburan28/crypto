# PS-1: algebra spec, bit-degree and root-test formulations

Protocol: [`../protocols/PS-1.md`](../protocols/PS-1.md).  Plan: [`../PLAN.md`](../PLAN.md).  Baseline: [`../baseline_v1/`](../baseline_v1/README.md).
Run report with graphs and PDF: [`../reports/run-20261011T2217Z/`](../reports/run-20261011T2217Z/).  Class: **accounting**
(no timing here; the quantity is gates per pair).  Conductor task T-4.  Raw records: [`raw/`](raw/), commands:
[`run_ps1.sh`](run_ps1.sh).

## 1. Requested and executed

| PS-1 requirement | status | evidence |
|:--|:--|:--|
| H1: exact bit-degree of the quartic's coefficients | **verified complete** (measured, `l = 3…8`) | §2 |
| H2: a root test cheaper than `gcd(q, L_V mod q)` in counted gates | **verified complete** for one formulation (group law); the "trace-based affine quadratic" and "resolvent cubic" candidates were **not attempted** | §4 |
| State each candidate as a straight-line program and count its gates | **partial**: the group-law test is a program (run on 64 lanes and on a gate-counting lane, same code); the baseline formulation's count is **derived** from the counted operations of PS-0, not built as a program | §4 |
| Check equality with direct evaluation of `S₄` and with `gcd(q, L_V mod q)` | **verified complete**, exhaustive for `n ≤ 9`, sampled for `n = 11, 13`, and on the PS-0 corpus | §3 |
| Exhaustive agreement at `n ≤ 9` over all triples | **verified complete**: every field of degree 5…9, every `l ≤ n`, every rational target | §3.2 |
| Sampled agreement at `n = 11, 13` over `10^6` random pairs "against the PS-0 baseline" | **verified complete with one change**: the comparison is against the generic summation-polynomial enumerator (the PS-0 root-set generator), not against `decompose`, because the pairs are drawn at `(n, l)` that are not the PS-0 rungs (`n ≠ 3l`); the PS-0 rungs themselves are covered by the corpus check | §3.2 |
| Every identity relied on gets an `identity.certificate/v1` record | **blocked**: the schema cannot express a characteristic-2 identity | §5 |
| Cheaper formulation and its gate count named for PS-2 | **verified complete** | §4 |
| Stop: bit-degree ≥ 5; no formulation under `1.8·10^5` gates per pair | **not triggered** (maximum bit-degree 4; the group-law test is `7.4·10^3 … 1.5·10^4`) | §2, §4 |

## 2. Bit-degree of the quartic's coefficients (H1)

`S₄(x₁, x₂, t, x_R)` as a quartic in `t`: `c₀ … c₄` below are the coefficients `quartic_in_x3` returns (`b = 1`).
Each output bit of each `cᵢ` is a Boolean function of the `2l` bits of `(x₁, x₂)`; its algebraic degree is the
largest monomial in its ANF (Möbius transform of the truth table).  Measured maximum over all `n = 3l` output
bits and over six random targets `x_R` (`raw/degrees.jsonl`):

| `l` | `n` | joint degree in the bits of `x₁, x₂`: `c₀ c₁ c₂ c₃ c₄` | degree in the bits of `x₂` alone, `x₁` fixed: `c₀ … c₄` |
|--:|--:|:--|:--|
| 3 | 9 | 2 4 2 4 2 | 1 2 1 2 1 |
| 4 | 12 | 2 4 2 4 2 | 1 2 1 2 1 |
| 5 | 15 | 2 4 2 4 2 | 1 2 1 2 1 |
| 6 | 18 | 2 4 2 4 2 | 1 2 1 2 1 |
| 7 | 21 | 2 4 2 4 2 | 1 2 1 2 1 |
| 8 | 24 | 2 4 2 4 2 | 1 2 1 2 1 |

The pattern is the one the formulas predict: with `s = x₁ + x₂` and `p = x₁x₂`, `s` and `s²` are linear, `p`, `p²` and
`s²p`, `p²s²` are bilinear in the bits (degree 2), and only `p³ = x₁³x₂³` reaches degree 4 (it appears in `c₁` and
`c₃`).  With `x₁` fixed, `p` is linear in the bits of `x₂`, so `c₀, c₂, c₄` are affine and `c₁, c₃` quadratic.
**H1 holds**: every coefficient has bit-degree at most 4, at most 2 for `c₀, c₂, c₄`, and at most 2 in the bits
of the inner variable.  The stop condition (degree ≥ 5) is not met.  What it implies for a Gray-code
incremental update is *measured later* (PS-2): this task does not show that forming the quartic is where the
baseline's time goes.

## 3. The group-law formulation

**Reference (what PS-0 times).** `S₄ = 0` over `V³` by forming the quartic for each pair and `gcd(q, L_V mod q)`.

**Candidate F3, the group-law test.**  A rational base point is `(x, y)` with `x ∈ V`.  For a target `R` and a
base point `P₁` let `U₊ = P₁ + R`, `U₋ = P₁ − R` (two points per row, computed once).  The four candidate
third abscissae of the pair `(P₁, P₂)` are `x(U₊ ± P₂)` and `x(U₋ ± P₂)`, and with `λ = (y_U + y₂)/(x_U + x₂)`

```
  x(U + P₂) = λ² + λ + (x_U + x₂) + 1
  x(U − P₂) = λ̃² + λ̃ + (x_U + x₂) + 1,        λ̃ = λ + x₂/(x_U + x₂)
```

(`a = b = 1`).  A pair hits when some candidate lies in `V`, which for the low-order subspace is a test of the top
`n − l` bits.  Per pair that is two inversions (taken a row at a time with Montgomery's trick), four multiplications and
four squarings; no polynomial is formed, no `L_V` is reduced.  It decides the question over the *rational* points of
`V`; §3.1 says when that is the baseline's question.

## 3. What the group-law formulation decides, and when it is the baseline's question

Notation.  `E: y² + xy = x³ + x² + 1` over `F = F_{2ⁿ}`, `σ` the Frobenius,
`T = (0, 1)` the rational point of order 2, `E[2] = {O, T}`.  `V ⊂ F` is the
factor-base subspace.  `R` is a rational point with `x_R ≠ 0`.  An
abscissa `x ∈ F` is *liftable* if some `y ∈ F` has `(x, y) ∈ E`.

**Lemma 1 (Semaev).**  For finite `x₁…x₄` in an algebraic closure,
`S₄(x₁,…,x₄) = 0` iff there are `y_i` with `(x_i, y_i) ∈ E` and
`P₁ + P₂ + P₃ + P₄ = O`.  (The polynomial is the one the baseline evaluates;
`binary_semaev_s4` and the baseline's tests tie it to the group law.)

**Lemma 2.**  A non-liftable `x ∈ F` is the abscissa of two points `±P` over
`F_{q²}` with `σP = −P`.  A liftable one has two rational points `±P`.
*Proof.* `y` and `y^q` are the two roots of `Y² + xY + (x³+x²+1)`, so
`y^q = y + x`, i.e. `σP = (x, y + x) = −P`. ∎

**Theorem (roots are rational).**  Let `x_R ≠ 0` be the abscissa of a
rational point `R`, and `x₁, x₂, x₃ ∈ F` with `S₄(x₁,x₂,x₃,x_R) = 0`.  Let `N`
be the set of indices `i` with `x_i` non-liftable.  Then `|N| ∈ {0, 2}`, and
if `|N| = 2`, with `N = {a, b}` and `c` the other index, either

* (F1) `x_a = x_b` and `x_c = x_R`, or
* (F2) `x_a·x_b = 1` and `x_c·x_R = 1`.

*Proof.* By Lemma 1 choose `P_i` over `F̄` with `P₁ + P₂ + P₃ = −R'`,
`R' ∈ {R, −R}`, a rational point.  By Lemma 2, `σP_i = ε_i P_i` with
`ε_i = −1` for `i ∈ N` and `+1` otherwise (either sign choice of `P_i`).
Applying `σ` gives `Σ ε_i P_i = −R'`; subtracting, `2 Σ_{i∈N} P_i = O`.

* `|N| = 1`: `P_i ∈ E[2]`, but both points of `E[2]` are rational (`O` is
  excluded, `T` has `x = 0`, which lifts), contradicting `i ∈ N`.
* `|N| = 3`: `P₁ + P₂ + P₃ ∈ E[2]`, hence `R' ∈ E[2]`, hence `x_R = 0`;
  excluded.
* `|N| = 2`: `P_a + P_b = τ ∈ E[2]`.  If `τ = O`, `P_b = −P_a` (so
  `x_a = x_b`) and `P_c = −R'` (so `x_c = x_R`).  If `τ = T`,
  `P_b = T − P_a`; translation by `T` is `x ↦ 1/x` on this curve
  (`b = 1`), checked below, so `x_b = 1/x_a`, and `P_c = −R' − T` gives
  `x_c = 1/x_R`. ∎

**Corollary.**  If `x_R ∉ V` and `1/x_R ∉ V` (and `x_R ≠ 0`), every root of
`S₄(·,·,·,x_R)` in `V³` has three liftable abscissae, hence is a rational
decomposition `±P_a ± P_b ± P_c = ±R`; conversely every such decomposition
is a root (Lemma 1).  So the group-law enumeration and the baseline return
the same root set.  When `x_R ∈ V` or `1/x_R ∈ V` the baseline may return
extra roots of the two families above, all with two *non-liftable*
abscissae, and the group-law formulation must not be used: the kernels call
the baseline for such a target.  For a uniformly random point of
`F_{2ⁿ}` the probability is `2·2^{ℓ−n}`.

(`x(P + T) = 1/x` is not assumed: §3.2 verifies it on every liftable `x ≠ 0` in every field of degree 5 to 12 and on samples at the PS-0 field degrees.)

### 3.2 What was checked

All commands are in [`run_ps1.sh`](run_ps1.sh); every record is in [`raw/`](raw/).

* **Translation by `T`.**  `x(P + T) = 1/x` and `T + T = O` hold for every liftable `x ≠ 0` and both signs, in every
  field of degree 5…12 (exhaustive) and, at 200,000 sampled points each, at `n = 24, 27, 30, 33, 36`
  (`raw/torsion.jsonl`: 0 mismatches).
* **Exhaustive agreement, `n ≤ 9`.**  For each of the 35 pairs `(n, l)`, `n ∈ {5,…,9}`, `l ∈ {1,…,n}` (`l = n` makes `V` the
  whole field), every rational target `x_R ≠ 0` (4168 targets): the root set of the generic summation-polynomial
  enumerator (independent of the baseline's gcd machinery, §PS-0) against the exact group-law reference.  **2868
  targets agree exactly; 1300 differ, and every difference is an extra root of the family of the theorem (F1/F2) at a
  target with `x_R ∈ V` or `1/x_R ∈ V`; 1334 targets are degenerate in that sense, 34 of them with no extra root.
  Differences at non-degenerate targets: 0; unexplained differences: 0** (`raw/agree.jsonl`; 15,628,190 roots in all).
* **The bit-sliced pipeline against the exact reference**, per pair, same 35 cases, every non-degenerate rational
  target (2834 targets, 1,259,432 pairs): the 64-lane program (the code that was gate-counted) and the exact
  reference agree on whether the pair has a third in `V` for every pair that is not flagged exceptional;
  **0 mismatches**.  10,814 pairs (0.86%, small fields) have a zero denominator and are flagged for the exact path; those are
  not compared (`raw/bscheck.jsonl`).
* **Sampled agreement at `n = 11, 13`.**  `l = 3, 4, 5` at `n = 11` and `l = 4, 5, 6` at `n = 13`, `10^6` random pairs of base
  points each with a random rational target: the roots of the pair from the summation polynomial against the four
  group-law thirds: **0 disagreements in 6·10^6 pairs** (the numbers of those pairs that had a root: 2,051 / 16,265 / 46,834 / 4,706 / 13,141 / 28,234; `raw/sample.jsonl`).
* **The PS-0 corpus.**  The exact reference reproduces the full root set of every target of `baseline_v1` at
  `l = 8…12`: 400 targets, 0 differences, none degenerate (`raw/corpus_check.jsonl`).

## 4. Gates per pair

Gates are counted as **fused ternary-logic issues** (an XOR-AND is one gate; `fused` below) and, for the program, also unfused
(`unfused`: an XOR-AND counts two).  The field is the one of the rung (`n = 3l`, PS-0's modulus); a product is the schoolbook `n²`
XOR-ANDs plus the modular reduction; a squaring is the linear map `a ↦ a²` (`popcount − 1` XORs per output plane); an
inversion is the Itoh–Tsujii chain.

**One task of the group-law test** = one evaluation of the straight-line program: stage A (`d_U, d_V`, zero fix-up), the
Montgomery share of the inversion across a batch of 64 tasks and two chains, and stage C (four candidates and the
membership masks).  In the bit-sliced kernel one task serves all lanes of a vector at once, so the plan's `G` is this
count and the gates per pair are `G / lanes`.

| `l` | `n` | product | squaring | inversion | stage A | stage C | **task, fused** | task, unfused | baseline gcd test, derived | baseline / group law |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 8 | 24 | 668 | 41 | 5,475 | 99 | 3,195 | **7,410** | 12,766 | 70,695 | 9.5 |
| 9 | 27 | 833 | 46 | 6,013 | 111 | 3,919 | **9,137** | 15,961 | 94,729 | 10.4 |
| 10 | 30 | 958 | 15 | 7,493 | 123 | 4,339 | **10,354** | 18,869 | 118,099 | 11.4 |
| 11 | 33 | 1,153 | 21 | 6,958 | 135 | 5,187 | **12,349** | 22,653 | 151,371 | 12.3 |
| 12 | 36 | 1,366 | 22 | 9,813 | 147 | 6,087 | **14,608** | 27,011 | 190,275 | 13.0 |

The squaring cost varies with the modulus (a trinomial with a small middle term is cheap).  **The baseline's number is derived, not a
program**: per pair, its counted operations (PS-0 `opcount/counts.jsonl`: `mul`, `batch_mul` at the cost of a product, `sqr`, stand-alone
`clmul` at `n²`, `reduce` at the reduction part of a product, and `(clmul + sqr_wide − reduce)` extra `2n − 1`-plane XORs to sum
unreduced products; `sqr_wide` is free wiring).  It assumes the generic case of the baseline's straight-line code (no vanishing
leading coefficient).

Against the plan's thresholds (`PLAN.md` §5): the baseline formulation is `7.1·10^4` at `n = 24` and
`1.5·10^5` at `n = 33`, above the `9·10^4` that a 10× over 310 ns needs at `n = 33` and, at `n = 36`, above `1.8·10^5` (5×
unreachable by bit-slicing it); the group-law test is `1.2·10^4` at `n = 33`, below both.  **The formulation for PS-2 is the group-law
test (F3).**  This is a statement about gates in the program model, not about wall time.

The program's code (`examples/pdp_speedup/{bitslice,bspairs,gatecount}.rs`) runs unchanged on `u64` (64 lanes), on the counting lane and on a
program-recording lane; PS-2 exports and hashes the recorded program.

## 5. Identity certificates: blocked, with the exact gap

The protocol requires an `identity.certificate/v1` record for every identity relied on.  The schema
(`docs/identity-certificates/README.md`) fixes `modulus = 2^61 − 1` and evaluates both sides over `F_p`; every
identity here (the coefficient formulas, the row-point identity, `x(P + T) = 1/x`) holds in characteristic 2 over
`F_{2^n}` and has no meaning modulo `2^61 − 1`.  **Original statement:** every identity relied on gets a certificate.
**Issue:** the format cannot represent them, so none can be issued.  **Proposed change:** a `v2` schema whose evaluation field is
`F_{2^n}` (modulus polynomial in the record, points drawn by the same SHA-256/splitmix64 rule).  **Substituted evidence**
(not a certificate): exhaustive agreement over every field of degree 5–9 (§3.2), which is stronger than Schwartz–Zippel
sampling at those sizes, plus sampled agreement at degrees 11, 13 and the PS-0 corpus at degrees 24–36.  PS-1 is therefore
not marked complete on this requirement.

## 6. Not claimed

No wall-time statement; the group-law test is a different algorithm from the baseline's, decides the rational question,
and applies to targets outside `x_R ∈ V ∪ V⁻¹ ∪ {0}` (PS-2 and later must call the baseline for those, and count it).
Nothing transfers to ECC2K-130.  The incremental-update hypothesis of the plan (Gray-code updates of the quartic) is not
tested here.  Nothing was measured on Arm64, Metal or an FPGA.
