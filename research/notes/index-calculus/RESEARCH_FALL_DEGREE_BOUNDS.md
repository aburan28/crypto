# Bounds on the first fall degree and the last fall degree of index-calculus systems

**Status:** literature ledger plus a toy measurement, 2026-10-08.
**Question answered:** what can be *proved* (and what is only measured or
heuristic) about the first fall degree `d_ff` and the last fall degree
`d_F` of the Weil-descended summation-polynomial systems that index
calculus on `E/F_{2^n}` (and, by the same arguments, `E/F_{p^k}`) has to
solve.
**Raw data:** `research/fall_degree_bounds_20261008/` (Sage scripts,
JSONL per cell). Experiments that test the claims below are preregistered
in its `PREREGISTRATION.md` and scored in its `RESULTS.md`.
**Builds on:** `RESEARCH_FFD_MEASUREMENT.md`, `RESEARCH_DREG_MEASUREMENT.md`,
`RESEARCH_FFD_PROOF_COMPLEXITY.md`,
`research/dreg_degree7_ell5_20260930/`.

---

## 0. Definitions used here (one convention, stated once)

Work over `k = F_q`, ring `R = k[x_1..x_N]`, system `F ⊂ R`; for
Boolean systems the field equations `x_i^2 − x_i` are *in* `F`.

- `V_{F,c}` (Huang–Kosters–Yeo, HKY): the smallest `k`-space with
  `F ∩ R_{≤c} ⊆ V_{F,c}` and `g ∈ V_{F,c}, deg(hg) ≤ c ⇒ hg ∈ V_{F,c}`.
  This is the row space of the degree-`c` Macaulay matrix *after*
  iterating degree falls (MutantXL closure); Caminata–Gorla prove
  `V_{F,c} = rowsp(M_c)` for every degree-compatible order.
- **Last fall degree** `d_F = min{c : f ∈ V_{F,max(c,deg f)} ∀ f ∈ (F)}`.
  Equivalently (HKY Prop. 2.6 iii, Caminata–Gorla Thm 2.8) the *largest*
  `c` with `V_{F,c} ∩ R_{≤c−1} ≠ V_{F,c−1}`: the last degree at which a
  fall is still needed.
- **First fall degree**, operational (HKY/HKYY): the *smallest* `c` with
  `V_{F,c} ∩ R_{≤c−1} ≠ V_{F,c−1}`. The syzygy-based definition
  (Hodges–Petit–Schlatter, Caminata–Gorla: first `d` with a non-trivial
  syzygy of `F^top` in degree `d`) is *not* the same number; see §1.
- **Solving degree** `sd_σ(F)`: least `d` such that `rowsp(M_d)` contains
  a `σ`-Gröbner basis. **Degree of regularity** `d_reg`: Bardet–Faugère–
  Salvy, first `d` with `(F^top)_d = R_d`.
- In the Boolean ring `B = F_2[x]/(x_i^2 − x_i)` with the rule "multiply
  `g` by `m` only if `deg m + deg g ≤ c`", the closure equals the
  reduction of `V_{F∪FE,c}`; so `d_F` computed in `B` *is* HKY's `d_F`
  for `F ∪ {field equations}`. (Proof: reduction mod field equations is
  a degree-non-increasing map that stays inside `V_{F∪FE,c}`, and the two
  closures are each contained in the other by induction on the
  construction.) This is the convention of the Sage script.

## 1. General bounds (THEOREM level, any system)

| statement | status | source |
|---|---|---|
| `sd_σ(F) = max{ d_F , max.GB.deg_σ(F) }` for every degree-compatible `σ` | THEOREM | Caminata–Gorla, arXiv:2112.05579, Thm 3.1 |
| `d_F = max{ d_g : g in a reduced GB, deg g < d_g }` where `d_g = min{c : g ∈ V_{F,c}}` | THEOREM | Caminata–Gorla Thm 2.9 |
| `d_F ≤` the largest degree reached by any F4/F5/XL run that terminates, in particular `d_F ≤ d_reg` when the latter is finite | THEOREM | HKY Prop 2.6 ii, HKYY §1 |
| `sd(F) ≤ reg(F^h)` (Castelnuovo–Mumford regularity of the homogenisation) for `F` containing the field equations | THEOREM | Caminata–Gorla 2021, cited as Thm 5.1 in arXiv:2112.05579 |
| `sd(F) ≤ 2·d_reg(F) − 2` when `F` contains the field equations and `d_reg(F) ≥ max(q, deg f_i)`; hence `d_F ≤ 2·d_reg − 2` | THEOREM (the only proved link from the Bardet–Faugère–Salvy degree of regularity to the solving degree) | Semaev–Tenti, J. Algebra 2021 (Thm 2.1), via Caminata–Gorla Thm 5.2 |
| `d_ff(operational) ≤ d_F` whenever a fall exists | THEOREM (one line from the two definitions) | HKY/HKYY §2.2 |
| `d_ff(syzygy) ≶ d_F`, `d_ff(syzygy) ≶ sd`, `d_ff(syzygy) ≶ d_reg`, with arbitrarily large gaps both ways | THEOREM (explicit examples) | Caminata–Gorla §4, Ex. 4.1–4.4 |
| zero-dimensional radical `F` with ≤ `e` solutions is solved by linear algebra in `V_{F,max(d_F,e)}` plus univariate factoring | THEOREM | HKY Prop 2.8 |
| Boolean system (field equations in `F`): `d_F ≤ N + max deg F` (in general `N(q−1) + max deg F`) | THEOREM (trivial: every reduced cofactor has degree ≤ `N(q−1)`) | this note |
| For an **unsatisfiable** Boolean system, `d_F` = the Polynomial-Calculus refutation degree of `F` with Boolean axioms | THEOREM (the PC derivation rules *are* the `V_{F,c}` closure rules) | Clegg–Edmonds–Impagliazzo 1996 + Def. above |

Consequences. (a) `d_F` is the quantity a MutantXL-type solver actually
pays for; plain XL pays `sd ≥ d_F`. (b) Any *lower* bound on `d_F` for an
unsatisfiable instance is a PC-degree lower bound, and vice versa; this is
the bridge `RESEARCH_FFD_PROOF_COMPLEXITY.md` builds on. (c) The
first fall degree, in either definition, bounds nothing above: it is only
a lower bound on `d_F` (operational) or not comparable at all (syzygy).

## 2. Index-calculus systems: what is proved

Setting: `S_{m+1}(x_1,…,x_m, x(R)) = 0`, `x_i ∈ V`, `V ⊂ F_{2^n}` an
`F_2`-subspace of dimension `n'` (`m n' ≈ n`), Weil-descended to `n`
equations in `N = m n'` Boolean variables.

### 2.1 First fall degree — upper bounds

| `m` | bound | status | source |
|---|---|---|---|
| 2 | `d_ff = 2` for ordinary curves, in general | THEOREM for the *relation* (an explicit `F_2`-linear combination of the descended equations equals a degree-1 polynomial, coming from the trace morphism `E(F) → F_2`, `P ↦ Tr((x(P)+a_2)/a_1^2)`); "usually a fall" because the quadratic parts cancel in practice | Kosters–Yeo, arXiv:1503.08001, Cor 4.11, Rem 4.12 |
| ≥ 3, `n' ≥ m` | `d_ff ≤ m² − m + 1` | THEOREM | Kousidis–Wiemers, J. Math. Cryptol. 2019 (arXiv:1906.05594), Thm 3.2 |
| any | `d_ff ≤ m² + 1` | THEOREM (older, weaker) | Petit–Quisquater ASIACRYPT 2012, Prop 1 |
| any descended `f` | `d_ff ≤ 1 + (max over monomials of the sum of the base-`q` digit sums of its exponents)` | THEOREM (the general Weil-descent weight bound; the two rows above are its specialisations) | Hodges–Petit–Schlatter, FFA 2014 |

| 2, subfield base `V = F_{2^{n/2}}` | `dim(span_{F_2}(F) ∩ R_{≤1}) ≥ dim span_{F_2}(F) − n'`, all in the bits of `X_1 + X_2`; measured `n' − 1 + Tr(b/x_3²)` on 39/40 E1 draws, so `d_ff = 2` and degree 2 already pins all but about one bit of `X_1 + X_2` | THEOREM for the inequality (§5, Prop. B); exact count TOY-EVIDENCE |
| 2, any `V` | the Kosters–Yeo relation restricted to `V` is `L = Tr(X_1) + Tr(X_2) + Tr(b/x_3²)` for `y² + xy = x³ + b`; so `d_ff = 2` iff `Tr|_V ≢ 0` or another fall exists, and when `Tr|_V ≡ 0` and `Tr(b/x_3²) = 1` the target is refuted at degree 2 | THEOREM given Kosters–Yeo Cor 4.11 (§5, Prop. A); 112/112 E1 draws, 17/17 refutations |

None of these depends on `n`. Kousidis–Wiemers conjecture `m² − m + 1`
is sharp for `m ≥ 3` and measured `d_ff = 7 (m=3), 13 (m=4)` at
`n ≤ 21`; their **first fall degree equals their observed F4 top degree**
at those sizes, which is the whole empirical basis of the "first fall
degree assumption".

### 2.2 First fall degree — lower bounds

Nothing beyond the trivial `d_ff ≥ 2` is proved for `m ≥ 3`. For `m = 2`
the value `2` is attained, so the bound is tight and *useless* as a cost
proxy: Kosters–Yeo (§5) show that if "`d_reg ≈ d_ff`" held for these
systems, the split system `{S_3(a_1,a_2,X_1), S_3(a_3,X_1,X_2), …}` would
decide "`S_m(a) = 0`" in polynomial time, a problem they prove NP-complete.
So any lower bound that matters must be on `d_F`, not `d_ff`.

### 2.3 Last fall degree — upper bounds

| bound | status | source / remark |
|---|---|---|
| `d_F ≤ N + 2 = m n' + 2` for the descended `S_3` system (`N + deg` in general) | THEOREM (trivial) | §1 |
| `d_{F'} ≤ max( τ(max(d_F, deg F, (m+1)s), q, m), m·τ(2s, q, 1), q )`, `τ(r,c,t) = ⌊2t(c−1) log_c(r/2t) + 1⌋`, for the descent of a zero-dimensional radical `F ⊂ F_{q^n}[X_1..X_m]` with ≤ `s` solutions and an injective coordinate projection | THEOREM | Huang–Kosters–Yang–Yeo, arXiv:1505.02532, Thm 1.1 (HFE case `m = 1`: HKY CRYPTO 2015, Thm 4.5) |
| …applied to `F = {S_3, L_V(X_1), L_V(X_2)}` (`L_V` = vanishing polynomial of `V`, degree `2^{n'}`, so `deg F = 2^{n'}`, `s ≤ 2^{n'+1}`): `d_{F'} ≤ 4n' + 6 ≈ 2n + 6` | THEOREM modulo the injective-projection hypothesis, which fails for `S_3` (two `X_2` per `X_1`); HKYY §4/§6 say the condition "might" be removable | computed in `papers/`; numeric table in §3 |
| `sd(Weil(F)) ≤ n·reg(F^h) − n + 1` | THEOREM, but `reg(F^h) ≥ 2^{n'}` here, so exponential | Caminata–Ceria–Gorla, arXiv:2112.10506, Thm 4.7 / Cor 4.9 |
| any `n`-independent bound | **none proved** for the subspace factor base; HKYY §6 explain why: the big-field system is not zero-dimensional until the subspace constraint (degree `2^{n'}`) is added | — |
| quasi-subfield variant: replace `V` by the roots of `X^{q^{n'}} − λ(X)`, `deg λ = q^{n'·α}` small | not a fall-degree bound at all: Thm 3.2 prices Rojas' resultant method, `m!·q^{n−n'm+n'}·Õ(m^{5.188}(3d)^{4.876 m²} + m q^{2n'})`, under a zero-dimensionality assumption; `RESEARCH_QUASI_SUBFIELD.md` shows the needed `λ` do not exist in the useful range | Huang–Kosters–Petit–Yeo–Yun, JMC 2020, Thm 3.2 |

So the only *proved* upper bounds on the last fall degree of the actual
index-calculus system are linear in `n`, i.e. they certify nothing better
than exponential cost. The subexponential claims (Petit–Quisquater
`2^{O(n^{2/3} log n)}`, Semaev 2015, Kousidis–Wiemers' sharpened
`2^{c log n (n^{2/3} − n^{1/3} + 1)}`) all rest on the **unproved**
`d_F ≈ d_ff` (resp. `d_reg = m² − m + 1 + o(1)`).

### 2.4 Last fall degree — lower bounds

| statement | status |
|---|---|
| `d_F ≥ d_ff(operational)`; for `m = 2` this is `≥ 2`, for `m ≥ 3` nothing better than `≥ 2` is proved | THEOREM / trivial |
| `d_F ≥ d_g` for every reduced-GB element with a fall; for unsatisfiable targets `d_F` is exactly the degree at which `1` enters the mutant closure (= PC refutation degree) | THEOREM |
| Conditional: if ECDLP over `F_{2^n}` has no `2^{o(n)}` algorithm, then for `m = n^{1/3}`, `n' = n^{2/3}`, the last fall degree of the descended `S_{m+1}` system is `ω(m²)`, because `d_F = O(m²)` with HKY Prop 2.8 gives `2^{O(n^{2/3} log n)}` | MODEL-BOUND (conditional, not a theorem about `d_F`) |
| Unconditional `d_F = Ω(n)` (or even `ω(1)`) for the generic subspace factor base | **OPEN**. The route in `RESEARCH_FFD_PROOF_COMPLEXITY.md` (PC-degree lower bound after quotienting out the one trace fall) is the candidate proof obligation |

## 3. Toy measurement: the HKY last fall degree itself (new to this repo)

Earlier repo measurements record the first fall (`ffd_harness`) and the
*plain-Macaulay* refutation degree `D*` (`pc_degree_harness`,
`dreg_ladder`). Neither iterates falls, so neither is `d_F`. The Sage
script `lfd.sage` computes, on the descended `S_3` system restricted to
a random (or subfield) `V` of dimension `n' = n/2`:

- `ffd`: first `c` with a new element of degree `< c` in `rowsp(M_c)`;
- `xl`: first `c ≥ max.GB.deg` with `rowsp(M_c) = (F)_{≤c}` (plain XL,
  no mutants; an upper bound on `sd`);
- `d_last`: HKY `d_F` via the mutant closure and a deglex Gröbner basis.

Semi-regular reference for `n` quadrics in `n` Boolean unknowns
(Bardet–Faugère–Salvy series `(1+z)^n/(1+z^2)^n`): `D_reg = 3, 3, 4, 4, 4,
5, 5` for `n = 6, 8, 10, 12, 14, 16, 18`.

Measured cells (one value per draw; `research/fall_degree_bounds_20261008/TABLE.md` is regenerated by `table.py`):

| n | n' | family | target | draws | d_ff | plain-XL full-ideal degree | d_F (HKY last fall) | semi-regular D_reg |
|--:|--:|---|---|--:|---|---|---|--:|
| 8 | 4 | random | sat | 7 | 2 2 2 2 2 2 2 | 5 5 5 5 5 5 5 | 3 3 3 3 3 3 3 | 3 |
| 8 | 4 | random | unsat | 1 | 2 | 5 | 3 | 3 |
| 8 | 4 | subfield | sat | 6 | 2 2 2 2 2 2 | 5 5 5 5 5 5 | 3 3 3 3 3 3 | 3 |
| 8 | 4 | subfield | unsat | 2 | 2 2 | 5 5 | 2 2 | 3 |
| 10 | 5 | random | sat | 5 | 2 2 2 3 2 | 6 6 6 6 6 | 3 3 3 3 3 | 4 |
| 10 | 5 | random | unsat | 3 | 2 2 2 | 6 6 6 | 3 3 2 | 4 |
| 10 | 5 | subfield | sat | 5 | 2 2 2 2 2 | 6 6 6 6 6 | 3 3 3 3 3 | 4 |
| 10 | 5 | subfield | unsat | 3 | 2 2 2 | 6 6 6 | 2 2 3 | 4 |
| 12 | 6 | random | sat | 5 | 2 2 2 2 2 | 7 7 7 7 7 | 3 3 3 3 3 | 4 |
| 12 | 6 | random | unsat | 3 | 2 2 2 | 7 7 7 | 3 3 3 | 4 |
| 12 | 6 | subfield | sat | 5 | 2 2 2 2 2 | 7 7 7 7 7 | 3 3 3 3 3 | 4 |
| 12 | 6 | subfield | unsat | 3 | 2 2 2 | 7 7 7 | 3 2 3 | 4 |
| 14 | 7 | random | sat | 1 | 2 | 8 | 3 | 4 |
| 14 | 7 | random | unsat | 1 | 2 | 8 | 3 | 4 |
| 14 | 7 | subfield | sat | 2 | 2 2 | 8 8 | 3 3 | 4 |
| 16 | 8 | random | sat | 5 | 2 2 2 2 2 | – – – – – | 3 3 3 3 3 | 5 |
| 16 | 8 | random | unsat | 1 | 2 | – | 3 | 5 |
| 18 | 9 | random | sat | 3 | 2 2 2 | – – – | 4 4 4 | 5 |
| 18 | 9 | random | unsat | 3 | 2 2 2 | – – – | 4 4 4 | 5 |


Reading.

- `d_ff = 2` on almost every draw. The exceptions are random bases with
  `Tr|_V ≡ 0`, exactly as Prop. A predicts.
- `d_F` is **3** on every satisfiable cell for `n = 8 … 16` and **4** on
  every draw at `n = 18, 20, 22`. Unsatisfiable subfield targets with
  `Tr(b/x_3²) = 1` read 2 (Prop. A).
- The plain-XL degree at which the Macaulay space reaches the whole ideal
  is exactly `n' + 1` on every draw where it was measured. Without
  mutants, the tower pays the full subspace degree.
- Kousidis–Wiemers' Magma F4 runs (their Table 1, random subspace,
  `n' = n/2`) report a highest step degree of 4 at `n = 34–36` and 5 at
  `n = 37–48`. F4's elements all lie in `V_{F,D}` for its top step degree
  `D`, so by HKY Prop 2.6 ii these are **upper bounds** on `d_F`.
- Together: for `m = 2`, `d_F` steps 3 → 4 near `n = 18` and is at most
  5 through `n = 48`. A generic system of the same shape would need 8.
  Growth is real but slow. The data cannot separate `Θ(log n)` from a
  small linear slope.
- E2 (random-quadratic controls with a planted linear relation) sits
  **one degree above** Semaev at every `n = 10 … 18`. So the low `d_F` is
  structural beyond the trace fall, not generic.

TOY-EVIDENCE only.

## 4. Falsification criteria and next steps

- **Kill the "`d_F` is bounded" reading** (for `m = 2`): one more step
  `5 → 6` between `n = 48` and `n ≈ 64` with the same slope as `4 → 5`
  would make `d_F ≈ n/12 + 1` the better fit than `log_2 n`.
- **Kill the "`d_F = d_ff + O(1)`" assumption for `m = 3`**: the
  Kousidis–Wiemers `m = 3` cells (`d_ff = d_reg = 7` for `n ≤ 21`) need
  `n ≥ 40` to see a step; `dreg_ladder` on the chained system already
  shows the refutation degree stepping with `ℓ` at fixed `n`.
- Conservative: extend `lfd.sage` to `n = 16 … 24` (`d_F ≤ 5` keeps the
  closure within `C(24,≤5) ≈ 55k` columns), split by parity of `n`, and
  fit `d_F` against `log n` and `n`.
- Representation-changing: measure `d_F` of the *symmetrised* system.
  Adding the trace relation as a generator is a no-op: it already lies in
  `V_{F,2}`, so `V_{F,c}` is unchanged for `c ≥ 2` (HKY Prop 2.6 viii).
  The informative comparison is E2's planted-linear control, which
  isolates the trace fall from the rest of the Semaev structure.
- Speculative: prove a PC-degree lower bound for the descended system
  after removing the one known fall, via the Alekhnovich–Razborov
  immunity criterion applied to the multiplication tensor of a dense
  basis.

## 5. Two statements proved here (`m = 2`, characteristic 2)

Curve `E_{a_2}: y² + xy = x³ + a_2 x² + b` over `F = F_{2^n}`,
`S_3(x_1, x_2, x_3) = (x_1x_2 + x_1x_3 + x_2x_3)² + x_1x_2x_3 + b`. This
polynomial does not depend on `a_2`. `F` is the Weil descent of
`S_3(X_1, X_2, x_3)` with `X_i ∈ V`.

**Proposition A (trace relation, restricted to `V`).** Let `x_3 ≠ 0`.
Then `Tr(X_1) + Tr(X_2) + Tr(b/x_3²) ∈ span_{F_2}(F)` (after reducing
modulo the field equations). Restricted to `V` with basis `v_j` and bits
`y_{ij}`, this is
`L = Σ_j Tr(v_j)(y_{1j} + y_{2j}) + Tr(b/x_3²)`.

*Proof.*

1. Choose `a_2 ∈ {0, ω}` with `Tr(ω) = 1` so that
   `Tr(x_3 + a_2 + b/x_3²) = 0`. Then `x_3` is the x-coordinate of a
   rational point `P` of `E_{a_2}`.
2. Kosters–Yeo Cor. 4.11 (with `a_1 = 1`, `a_3 = 0`) gives an explicit
   `F_2`-combination of the descended coordinates of `S_3(X_1, X_2, x(P))`
   equal to `Tr(x(P) + a_2) + Σ_j Tr(α_j)(X_{1j} + X_{2j})`. Here `α_j`
   is the basis of `F/F_2`, and the combination's coefficients are
   `Tr(α_j/b'²)` with `b' = x(P)`.
3. By the choice of `a_2`, `Tr(x_3 + a_2) = Tr(b/x_3²)`.
4. `Σ_j Tr(α_j) X_{ij} = Tr(X_i)` by linearity of the trace.
5. Restricting `X_i = Σ_j y_{ij} v_j` is `F_2`-linear and commutes with
   taking spans. ∎

*Consequences.*

- If `Tr|_V ≢ 0`, then `L` has degree 1 and is a combination of
  degree-2 rows, so `d_ff = 2`.
- If `Tr|_V ≡ 0`, then `L` is the constant `Tr(b/x_3²)`. When that
  constant is 1, `1 ∈ V_{F,2}`, so the target is not decomposable and
  `d_F ≤ 2`. A uniformly random target has this constant equal to 1 for
  about half of all draws.
- The EXP-J degree-3 syzygy is `(L+1)·L = 0` (E1b).

**Proposition B (subfield base).** Let `n = 2n'` and
`V = F_{2^{n'}} ⊂ F`. Then
`dim(span_{F_2}(F) ∩ R_{≤1}) ≥ dim span_{F_2}(F) − n'`, and every such
element is affine-linear in the bits of `s_1 = X_1 + X_2`. In particular,
if the `n` descended equations are `F_2`-independent, there are at least
`n'` of them.

*Proof.*

1. With `s_1 = X_1 + X_2` and `s_2 = X_1X_2`, both in `V`:
   `S_3 = s_2² + s_1² x_3² + s_2 x_3 + b`.
2. Over `F_2`, the bits of `s_1²` are linear in the bits of `s_1`, which
   are linear in the `y_{ij}`.
3. The bits of `s_2` and `s_2²` are `F_2`-linear combinations of the `n'`
   bits of `s_2 ∈ V`, each a quadratic form in the `y_{ij}`.
4. So the quadratic parts of the `n` descended equations lie in an
   `F_2`-space of dimension at most `n'`.
5. The map from `span_{F_2}(F)` to quadratic parts therefore has rank
   at most `n'`. Its kernel, of dimension at least
   `dim span_{F_2}(F) − n'`, consists of the combinations with no
   quadratic part. Each lies in the span of the bits of `s_1²` and the
   constant, so each is affine-linear in the bits of `s_1`. ∎

The `n` equations are *not* always independent, and Prop. A says when.

- For the subfield, `Tr|_V ≡ 0`. When `Tr(b/x_3²) = 0`, the Kosters–Yeo
  combination is identically zero, a dependency, so the bound gives
  `n' − 1`.
- When the constant is 1, the bound gives `n'`, including the constant
  1 itself.

E1 matches this exactly on 39 of 40 subfield draws:

| `Tr(b/x_3²)` | relations of degree ≤ 1 | draws |
|---|---|--:|
| 1 | `n'` (the constant 1 and `n' − 1` linear) | 16 |
| 0 | `n' − 1` | 23 |
| 0, with `x_3 ∈ F_{2^{n'}}` | 2 (`n = 16`, seed 112005) | 1 |

The outlier has every coefficient of `S_3` in the subfield. The descent
then has only `n'` independent coordinates, and the target is refuted at
degree 2.

This is the algebraic reason the subfield case is easy. At degree 2 the
system already pins all but about one bit of `X_1 + X_2`, as in the
classical subfield decomposition.

- The inequality is a THEOREM.
- The exact count `n' − 1 + Tr(b/x_3²)` for `x_3 ∉ F_{2^{n'}}` is
  TOY-EVIDENCE. Its proof needs the dependency of the `n` equations to be
  exactly the one from Prop. A.

## References

- M.-D. Huang, M. Kosters, S. L. Yeo, *Last fall degree, HFE, and Weil
  descent attacks on ECDLP*, CRYPTO 2015 (ePrint 2015/573).
- M.-D. Huang, M. Kosters, Y. Yang, S. L. Yeo, *On the last fall degree
  of zero-dimensional Weil descent systems*, arXiv:1505.02532.
- M. Kosters, S. L. Yeo, *Notes on summation polynomials*, arXiv:1503.08001.
- S. Kousidis, A. Wiemers, *On the first fall degree of summation
  polynomials*, J. Math. Cryptol. 13 (2019), arXiv:1906.05594.
- C. Petit, J.-J. Quisquater, *On polynomial systems arising from a Weil
  descent*, ASIACRYPT 2012.
- T. Hodges, C. Petit, J. Schlatter, *First fall degree and Weil
  descent*, Finite Fields Appl. 30 (2014).
- A. Caminata, E. Gorla, *Solving degree, last fall degree, and related
  invariants*, J. Symb. Comput. (2023), arXiv:2112.05579.
- A. Caminata, M. Ceria, E. Gorla, *The complexity of solving Weil
  restriction systems*, arXiv:2112.10506.
- M.-D. Huang, M. Kosters, C. Petit, S. L. Yeo, Y. Yun, *Quasi-subfield
  polynomials and the elliptic curve discrete logarithm problem*, J.
  Math. Cryptol. 14 (2020).
- M.-D. Huang, *On the last fall degree of Weil descent polynomial
  systems*, arXiv:2103.07282; *Last fall degree of semi-local polynomial
  systems*, arXiv:2311.02804 (bounds `d_F` by the number of closed-point
  solutions; not yet applied to summation systems).
- M. Clegg, J. Edmonds, R. Impagliazzo, *Using the Groebner basis
  algorithm to find proofs of unsatisfiability*, STOC 1996.
