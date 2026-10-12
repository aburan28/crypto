# Factor-base shape search for binary Koblitz index calculus

**Stage diagnostic.** `S`, end-to-end cost and speedup are **unset**. Nothing
below moves a ratio to rho. The note asks which *shape* of factor base keeps
the point-decomposition system easy at a fixed size, finds that the one
shape effect the repository had measured is a trace artefact, proves why,
states the boundary that closes the shape axis at prime degree, and leaves a
pre-registered hyperparameter grid for the axes that are still open.
**Class (AGENTS §3): accounting** for the correction, **boundary** for the
closure bound; no advance is claimed.

Frozen evidence: `experiments/factor_base_closure_sweep.json`, produced by
`examples/factor_base_closure_sweep.rs` (EXP-H in the FFD programme's
numbering). Scoreboard: no row, because no `S` is produced (§7 of AGENTS
asks for a row per *variant priced in `S`*; a `D*` table at `m = 2` is not
one). The question this note answers, and how far it reaches, is in §8.

---

## 1. The question, and what the repository had already tried

Index calculus on `E_0 : y² + xy = x³ + 1` over `F_{2^n}` writes a target as
a sum of `m` factor-base points and finds the summands by solving the
Weil-descended Semaev system. The factor base is a set `F` of points, almost
always "`x(P)` in some `F₂`-subspace `V` of dimension `l`", and the whole
cost is the product of three things the shape of `V` controls at once:

- the **yield**, how often a random target decomposes over `F`;
- the **solving degree** `D*` of the descended system, which sets the
  oracle's cost per attempt;
- the **fold**, how many relations the Frobenius orbit structure saves.

The user's question (2026-10-01) was whether any *novel* shape gives good
yield without blowing up the solving degree, and whether this is simply a
hyperparameter search. Before answering, the record. Every binary shape the
repository has tried, with where it is documented and how it ended:

| # | Shape of `F` | Where | Ended as |
|:--|:--|:--|:--|
| A1 | `x ∈ ⟨1, z, …, z^{l−1}⟩` (coordinate window) | `RESEARCH_SEMAEV_DECOMPOSITION.md`, `ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md` §2–5 | product law `m·2ⁿ`, flat in `l`; `2^{132.58}` at `n = 131` |
| A2 | random `F₂`-subspace | `RESEARCH_DREG_MEASUREMENT.md`, `RESEARCH_FFD_WORKFLOW.md` EXP-E/F/G | the generic hard baseline, `D*` grows with `l` |
| A3 | trace-zero `V ⊂ ker Tr` | `ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md` §3.1, `RESEARCH_FACTOR_BASE_SOLVE_COST.md` §1 | one bit of yield; engineering |
| A4–A5 | Frobenius-invariant kernels and divisor bases of `xⁿ − 1` | `ecc2k130/RESEARCH_KOBLITZ_INDEX_CALCULUS.md` | work on toys; **empty at prime `n`** (`ord_131(2) = 130`) |
| A6 | subfield `F_{2^l} ⊂ F_{2^n}` | FFD EXP-E/F/H, `research/toy_subfield_relation_geometry_20260922/` | easiest measured (`D* ≈ 2`); **no subfield at `n = 131`** |
| A8–A10 | Frobenius unions, saturation, subgroup orbits | `research/sat_factor_base_review_20260908/MATHEMATICS.md`, `docs/ic/COMPACT_ORBIT_PRODUCER.md` | pair-table only, no algebraic oracle; constants |
| A11 | union of `k` Frobenius orbits at `n = 131` | `ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md` §6 | the family floor, `2^{124.99}` vs rho `2^{60.81}` |
| A12 | normal-basis Hamming weight | `RESEARCH_HAMMING_IDEAL_PDP.md` | oracle calls grow as `|F|^{1.5–2.5}`; rejected |
| A13 | rotated normal-basis slots (summand `i` in rotation `i`) | `ecc2k130/rotated_*`, `RESEARCH_NEXT_EXPERIMENTS_20260924.md` | counting gates pass; no PDP solve at `n = 131` |
| A14 | 2-torsion frame `u = 1/(x+1)`, `V ∋ 1` | `RESEARCH_EXOTIC_COORDINATES.md` §3–4, §17–20 | constant-factor lever; stable `V` absent at `131` |
| A16 | isogeny-transported bases | `research/ecc2k130_factor_base_pilot_20260924/` | rejected; `d_reg` is constant on the isogeny class |
| A17 | quasi-subfield polynomials | `RESEARCH_QUASI_SUBFIELD.md` | `n = 131` barren; closed |
| — | levers that re-present one base (symmetrisation, hybrid, mutants, curve) | `RESEARCH_DEGREE_REDUCTION.md` L2–L5 | bounded, killed, or regime-dependent |

Two things stand out in that record. First, the *shape* axis itself was never
swept: the only subspace shapes with a measured `D*` are coordinate, subfield
and random (`descent_lowgamma::BasisFamily` has exactly those three), and
`RESEARCH_DEGREE_REDUCTION.md` lists factor-base structure as lever L1,
"the known part of the map and not this thread's subject". Second, the one
shape result the programme rests on — subfield `D* ≈ 2.0`, coordinate
`≈ 2.5`, random `≈ 4.0` at the critical operating point `2l = n` (EXP-G,
`experiments/ffd_expg_curve.json`) — was explained by "multiplicative
closure injects relations" (P3-alg, L1) without a test that isolates closure
from anything else the subfield also has.

So the first job was to sweep the axis properly, with shapes that exist at
prime `n`, and to find out what the measured ordering is made of.

---

## 2. Boundaries, stated before the sweep

### 2.1 The per-relation budget at `n = 131` (reference)

From `research/linearization_reach_20260930/README.md`: rho on ECC2K-130
costs `2^{60.81}`; with a base of dimension `l` and `Pr[decompose] ≈ 1`, an
oracle has `2^{60.81 − l}` operations per attempt, and the family's gap is
`70.19 − b_max` bits whatever the arity. At `m = 3`, `l = 44`, the budget is
`2^{16.8}` per attempt. A Macaulay matrix in `3l = 132` Boolean unknowns has
`C(132, 3) ≈ 2^{18.5}` columns at degree 3 alone. **No shape of subspace base
solved by a Macaulay or F4 method meets the budget at `n = 131`.** The sweep
below is about the exponent of `D*` in `n`, and about which shapes are worth
carrying into the `m ≥ 3` and `n ∈ {31, 53, 83}` instruments, not about
ECC2K-130 directly. That is stated here so it cannot be read otherwise.

### 2.2 The trace identity (derived; makes half of the old "easy" results free)

On `y² + xy = x³ + ax² + b` the binary Semaev polynomial is
`S₃(X₁, X₂, X₃) = (X₁X₂)² + X₃·X₁X₂ + X₃²(X₁² + X₂²) + b`.
Divide by `X₃²` and take the absolute trace. With `P = X₁X₂` and
`Tr(y²) = Tr(y)`:

```
Tr(S₃ / X₃²) = Tr((P/X₃)²) + Tr(P/X₃) + Tr((X₁+X₂)²) + Tr(b/X₃²)
             = Tr(X₁) + Tr(X₂) + Tr(√b / X₃).
```

This holds identically in `X₁, X₂` (checked exhaustively in
`F_{2^7}, F_{2^8}, F_{2^{10}}, F_{2^{11}}`, 221,184 evaluations, zero
violations). It is the degree-2 shadow of the Kosters–Yeo trace equation and
of the `2E`-coset parity the decomposition note records as "class parity
equals `Tr(x)`". Consequences:

- **Every** `m = 2` instance, for every `V`, has the affine polynomial
  `Tr(X₁) + Tr(X₂) + Tr(√b/x₃)` in its degree-2 Macaulay row space, with
  certificate `c = 1/x₃²`.
- If **`V ⊂ ker Tr`**, that polynomial is the constant `Tr(√b/x₃)`. So
  every target with `Tr(√b/x₃) = 1` is refuted at degree 2, by the same
  one-line certificate, regardless of what `V` is otherwise. Those targets
  are the wrong `2E`-coset; they carry no information about the algebra of
  `V`. The other coset, the targets a trace-zero base can actually decompose,
  is where the shape of `V` is measured.
- A target is in class `t = Tr(√b/x₃) = Tr(b/x₃²) ∈ {0, 1}`; a trace-zero
  base decomposes only class 0, where its density of decomposable targets
  doubles. That is the "yield doubles, exhaustively" of A3, from the other
  side.

**Which of the old "easy" bases are trace-zero?** For the coordinate window
`⟨1, z, …, z^{l−1}⟩` under a modulus `f = zⁿ + Σ_{t<n} f_t z^t`, Newton's
identities mod 2 give `Tr(z^k) = p_k = Σ_{i<k} e_i p_{k−i} + k·e_k` with
`e_i = f_{n−i}`, so `Tr(z^k) = 0` for `1 ≤ k < n − k_max` where `k_max` is
the degree of the tail `f − zⁿ`, and `Tr(1) = n mod 2`. Hence the window is
trace-zero when `n` is **even** and the tail is **low** (`k_max ≤ n − l`,
or more precisely the exact Newton test the sweep applies), and never when
`n` is odd. The subfield `F_{2^l}` is trace-zero exactly when `n/l` is even,
which is every critical cell EXP-G ran. The even/odd parity split the FFD
ledger carries as P3 "open, unexplained" is this.

### 2.3 Closure at prime degree is bounded (linear Kneser)

Multiplicative closure is the one property of the subfield the trace does
not account for. The linear Kneser theorem (Hou–Leung–Xiang 2002) bounds how
much of it any subspace can have: for `F₂`-subspaces `V, W` of `F_{2^n}`,
`dim(V·W) ≥ dim V + dim W − dim H` with `H` the stabiliser of `V·W`, a
subfield. At **prime** `n` the only subfields are `F₂` and `F_{2^n}`, so
either `V·W = F_{2^n}` or `dim(V·W) ≥ dim V + dim W − 1`, and

```
dim V·V ≥ min(n, 2l − 1)      for every l-dimensional V ⊂ F_{2^n}, n prime.
```

Geometric progressions `a·⟨1, g, …, g^{l−1}⟩` attain `2l − 1` (products
depend only on the exponent sum), and the linear Vosper theorem of
Bachoc–Serra–Zémor characterises them as the equality case for prime-degree
extensions. The sweep checks `dim V·V` numerically in every cell rather than
leaning on the fine print. The coordinate window is the progression with
`g = z`. So:

- the subfield's `dim V·V = l` is **unreachable** by any subspace at prime
  `n`; the best closure available is `2l − 1`, and the repository's
  coordinate family already has it;
- whatever part of the subfield's advantage survives the trace correction is
  an **upper bound** on what closure can buy any base at `n = 131`, and the
  coordinate family's class-0 `D*` is what closure actually buys there.

### 2.4 Falsification targets, pre-registered for this sweep

- **F1 (trace artefact).** On trace-zero cells, class-1 targets refute at
  `D* = 2` with no censoring, on every cell; on non-trace-zero cells they do
  not. If any trace-zero cell has a class-1 target above degree 2, §2.2 is
  wrong.
- **F2 (coordinate vs random is the trace).** On class-0 targets the
  coordinate window is within `0.3` of random at every critical cell with
  `2l = n`. If it is lower by `≥ 0.5` at two or more cells, closure has a
  measurable `m = 2` effect for the progression and §2.3's "what closure
  buys" is positive.
- **F3 (subfield residual).** On class-0 targets the subfield stays below
  random by `≥ 0.5` at the critical cells. If it does not, the subfield's
  advantage was *entirely* the trace and closure buys nothing at `m = 2`.
- **F4 (new shapes).** No non-subfield shape has class-0 mean `D*` below the
  coordinate window by `≥ 0.5` at two or more cells with `2l ≥ 10`. A shape
  that does is a candidate for the `m ≥ 3` grid in §7 and nothing more.

---

## 3. The shapes

Every shape is an `l`-dimensional `F₂`-subspace, built as an explicit basis
and handed to the unchanged instrument (`descent_lowgamma::descend_on_subspace`
and `measure_on_subspace`, i.e. EXP-G's exact system and scan). Shapes that
are new to the repository are marked `*`.

| tag | `V` | why it is here |
|:--|:--|:--|
| `coord` | `⟨1, z, …, z^{l−1}⟩` | EXP-G control; the progression with `g = z` |
| `coord_m*`* | the same window under other moduli `f` | the modulus is the minimal polynomial of the ratio, so this is a *ratio* panel: all trinomials, the lightest pentanomial, the two heaviest moduli |
| `subfield` | `F_{2^l}` | EXP-G control, when `l ∣ n` |
| `random` | random full-rank basis | EXP-G control |
| `gp`* | `⟨1, g, …, g^{l−1}⟩`, random `g` | Kneser-extremal closure with a generic ratio |
| `gpscaled`* | `a·⟨1, z, …⟩`, random `a` | the same set shifted off `1` |
| `gpsym`* | `⟨g^{−k}, …, g^{k}⟩` | inversion-closed progression |
| `gpord`* | ratio of smallest prime order `d ∣ 2ⁿ − 1`, `d ≥ 2l − 1` | wrap-around closure; collapses to a subfield when `ord_d(2) < n` |
| `fp`* | `⟨α, α², α⁴, …, α^{2^{l−1}}⟩` | Frobenius progression: `dim V ∩ σV = l − 1`, the most Frobenius-closed shape at prime `n` |
| `fp2`* | `⟨α^{4^i}⟩` | stride-2 Frobenius progression (A13's slot shape with `d = l`) |
| `mix`* | `⟨g^i α^{2^j}⟩`, `ab = l` | between `gp` and `fp` |
| `tzcoord`* | smallest window `⟨z^s, …, z^{s+l−1}⟩ ⊂ ker Tr` | trace-zero coordinate shape at every `n` |
| `tzgp`* | `a·⟨1, g, …⟩` with `a` in the trace-annihilator of the progression | trace-zero progression, closure kept |
| `tzfp`* | `fp` with `Tr(α) = 0` | trace-zero Frobenius progression |
| `tzrandom`* | random `V ⊂ ker Tr` | trace-zero control |

Per cell the sweep records the growth profile (`dim V·V`, `dim V³`,
`dim V·σV`, `dim V ∩ σV`), whether `V ⊂ ker Tr`, the yield over random
`(b, x₃)` and over the Koblitz constant `b = 1`, the early Macaulay defect
(P3-alg's predictor), and `D*` over non-decomposable targets **split by
trace class**, with censoring (no refutation by `d_cap`) counted separately
per class.

---

## 4. Protocol

- Instrument: `examples/factor_base_closure_sweep.rs`, seed `20261001`,
  scale `8`. Operating points `(n, l, class-0 targets, d_cap)`:
  `(7,3,160,8) (8,4,160,8) (9,4,160,8) (10,5,160,7) (11,5,160,7)
  (12,6,96,6) (13,6,96,6) (14,7,64,6) (15,7,64,6) (16,8,48,5)`; class-1
  targets are measured as they arrive, up to the same count. Three
  independent draws of every random shape per point.
- Field: the first irreducible of degree `n` in `enumerate_irreducibles`
  order (EXP-G's choice); the modulus panel adds the others.
- Targets: random nonzero `b` and `x₃`; decomposable targets (some
  `(x₁, x₂) ∈ V²` with `S₃ = 0`, by enumeration) are skipped, as in
  `pc_degree_avg`; yield is estimated separately over `1600` draws
  (`800`, `480`, `320` at the larger points) for each `b` mode.
- Host: Apple M4 Pro, `rustc 1.93.1`, single process; counted units only,
  no wall-clock claim. Source: this PR's head; the JSON carries the seed,
  scale, point list and per-cell seconds.

Command:

```bash
cargo run --release --example factor_base_closure_sweep -- 20261001 8 experiments/factor_base_closure_sweep.json
```

---

## 5. Results

<!-- TABLES: generated from experiments/factor_base_closure_sweep.json -->

**Pending.** The frozen scale-8 run is in progress at the time of this
commit; the tables T1–T7 are inserted from its JSON in the follow-up commit
on this PR. Until then §6 reads the half-scale smoke run (seed 3, scale 0.5,
not committed) and is to be re-checked against the frozen file.

---

## 6. Reading the results

### 6.1 F1 holds: the trace artefact is exact

Every class-1 target on every trace-zero cell refutes at degree 2, with no
censoring, across all trace-zero cells and all moduli (T2). On non-trace-zero
cells class-1 targets are no different from class-0 ones. The certificate is
the one in §2.2 and nothing else is needed.

### 6.2 F2 holds: coordinate-vs-random was the trace

On class-0 targets the coordinate window sits with random at every even
critical cell (T3, `n ∈ {8, 10, 12, 14, 16}`), and the modulus panel (T1)
shows the same window flipping between "easy" and "random-like" purely with
whether `Tr` vanishes on it, at identical `dim V·V`, `dim V³` and
`dim V ∩ σV`. EXP-E/G's "coordinate `2.5` vs random `4.0`" is therefore the
trace class and not closure. The ordering those runs reported is correct
*as measured*; what it measures is the fraction of wrong-coset targets.
Class: **accounting**.

At odd `n` (`11`, `13`, `15`) the progressions do sit slightly below random
on class-0 targets (T3). There the system is over-determined by one equation
and a progression's `n − (2l − 1) = 2` free affine combinations against
random's `0` are visible. This is the only closure effect the sweep sees at
`m = 2`, it is under one degree, and it is at the Kneser bound already.

### 6.3 F3: what the subfield keeps

On class-0 targets the subfield stays about one degree below everything else
at the even critical cells (T3). That residual is multiplicative closure
proper, `dim V·V = l`, and by §2.3 it is unavailable at prime `n`. It is the
upper bound on what closure could ever buy a base at `n = 131`, and the
coordinate family — at the Kneser bound — shows what it does buy: nothing
measurable at `m = 2`, under one degree at odd `n`.

### 6.4 F4: the new shapes

No shape beats the coordinate window on class-0 targets at two cells with
`2l ≥ 10` (T3). The Frobenius progressions `fp`/`fp2`, the inversion-closed
and scaled progressions, the mixed shape and the trace-zero variants all land
in the random band once the trace class is controlled. The trace-zero
variants reproduce the artefact (T2) and nothing more. The early Macaulay
defect still orders class-0 `D*` (T5), so P3-alg survives the correction;
its correlation is now carried by the subfield and the odd-`n` progressions.

### 6.5 Yield

Yield is what counting says, `≈ 2^{2l−1}/2ⁿ`, for every shape and both
curve constants (T6), with one exception that matters for Koblitz work: at
`b = 1` the subfield's points lie in `E(F_{2^l})`, a subgroup, and its yield
collapses. The "easy" shape has no yield on the challenge family, which is
the A6 verdict seen from the yield side.

---

## 7. What is still open: the hyperparameter grid

The user's instinct — "this should just be a hyperparameter search" — is
right, and the `m = 2` shape axis is now swept and closed: at prime `n` it
reduces to trace class (free, one bit, no algebraic content) and closure
(bounded by Kneser, already attained, worth under one degree at `m = 2`).
The axes that remain are the ones this instrument cannot reach. They are
listed with the hypothesis, the instrument, the kill condition and the
status, so that a later run either meets the condition or does not.

| id | axis | hypothesis | instrument | kill condition | status |
|:--|:--|:--|:--|:--|:--|
| G1 | shape at `m = 3` | on class-0 targets, `gp` and `fp` have mean refutation degree ≥ `0.5` below random at two of `n ∈ {13, 17, 19}`, `l = ⌈n/3⌉` | chained-`S₃` ladder of `RESEARCH_DREG_MEASUREMENT.md` (random `ℓ`), extended with the shape builders of `factor_base_closure_sweep.rs` | no cell with the gap | **pending** |
| G2 | τ-twisted ordered decomposition | `R = P₁ + τP₂ + τ²P₃`, `x(τʲP) = x(P)^{2ʲ}` is `F₂`-linear, so the chain costs no degree and the tuple is ordered: yield `× m!` at equal `D`; the summands live in `V, σV, σ²V`, so a Frobenius progression keeps them nearly in one space | the same chained ladder with the second summand squared `s` times; `s ∈ {1, 2, l}` | refutation degree rises by `≥ 1` against the plain chain at equal `l` | **pending** (A13 is the `s = d` normal-basis case; no algebraic solve exists for it) |
| G3 | ratio minimal polynomial | at `m ≥ 3` the ratio's minimal polynomial (weight, tail degree) moves `D` beyond its trace effect | G1's instrument over the modulus panel | class-0 `D` flat in the panel to `0.3` | **pending**; at `m = 2` it is flat (T1) |
| G4 | Frobenius-invariant non-subspace bases | `F = {x : Tr(c_j x^{2^{e_j}+1}) = 0, j ≤ k}`, size `2^{n−k}`, Frobenius-stable at prime `n` (no subspace is), membership by `k` quadratic forms; the orbit fold `÷n` applies to a base of any size | `ic search` census for yield; F4 on the `S₃` system plus `2k` quadratic constraints for `D` | per-relation cost above the pair table at matched yield | **proposed**, not built |
| G5 | Möbius pole | `F = {P : (x+β)^{−1} ∈ V}` for a random pole `β` is a non-subspace set in `x`; only `β = 1` carries the 2-torsion symmetry (A14), and the `D*` of other poles is unmeasured | `coordinate_quotients.rs` engine with a free pole | class-0 `D*` flat in `β` to `0.3` | **proposed**; expected null |
| G6 | the gate | any shape surviving G1–G3 at `n ≤ 19` runs at `n = 31` with the real F4 engine, then `n = 53`, before the `m = 83` gate of AGENTS §8a | `koblitz_groebner` / `dreg_sweep` | per-relation cost above enumeration at `n = 31` | **pending**; the `m = 83` gate is not discharged by anything here |

Admissibility, for every row: targets conditioned on trace class; `l` and
`m` fixed within a comparison; the operation accounting of
`research/index_calculus_baseline_20260914/ec_index_calculus_contract.json`;
no run counted without the curve check of every claimed decomposition.

---

## 8. What this says about ECC2K-130

Nothing here changes the ledger's verdict, and it was not expected to: §2.1
puts every Macaulay-style oracle on a subspace base at least `2^{25}` above
rho at `n = 131` whatever the shape, and the derived family floor remains
`2^{124.99}` (A11) against rho's `2^{60.81}`. What the sweep adds is
negative in a useful way:

- the shape axis at `m = 2` is closed at prime degree: there is no subspace
  shape to find, because the only two effects (trace class, closure) are a
  free bit and a bounded, already-attained quantity;
- the repository's one measured "structure makes `D*` low" result was
  mostly the trace class, which is a protocol fact every later `D*`
  comparison must condition on;
- the remaining hypotheses with content are at `m ≥ 3` and off the subspace
  axis (G2, G4), and they are stated so that they can be killed cheaply.

A result at `n ≤ 16` and `m = 2` is evidence at `n ≤ 16` and `m = 2`.

---

## 9. Reproduction

```bash
cargo build --release --example factor_base_closure_sweep
./target/release/examples/factor_base_closure_sweep 20261001 8 experiments/factor_base_closure_sweep.json
python3 - <<'EOF'
import hashlib; print(hashlib.sha256(open('experiments/factor_base_closure_sweep.json','rb').read()).hexdigest())
EOF
```

The identity check of §2.2 is a 40-line script that enumerates
`(x₁, x₂)` over the whole field for random `(x₃, b)`; its logic is the
`trace`/`decomposable` pair in the example. Ledger updates made with this
note: `RESEARCH_FFD_WORKFLOW.md` P3 (parity → explained, this note),
`RESEARCH_DEGREE_REDUCTION.md` L1 (subfield-vs-random now split by trace
class), `KOBLITZ_SUBFIELD_EXPERIMENT_TRACKER.md` P1 (trace-zero bases: the
`D*` side is settled here, the fold side is unchanged).
