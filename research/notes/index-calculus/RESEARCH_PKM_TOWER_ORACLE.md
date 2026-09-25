# The PKM tower oracle: an algebraic decomposition oracle for prime-field curves

**Status:** design and pre-registration (§§0–9) followed by a design-validation
pilot (§10, 2026-09-24).
- §§0–9 were fixed before any run, per `AGENTS.md` §1 and §4, in commit
  `c2157282`, and are unchanged since.
- **The pilot contradicts the pre-registered expectation.** F4's solving degree
  does not grow like `N`. In all three tower families it stays at 4–5 for `m = 2`
  through `N = 18`, and at 5–6 for `m = 3` through `N = 12`. At `m = 4` it is 6
  and 7 at the only two sizes run, `N = 8` and 12.
- Under §5.2 this reads *inconclusive*, not H1. The shape-matched null behaves
  the same way, and the plateau is too short to tell a bounded degree from a
  slowly growing one.
- §10 says what the pilot does and does not show.
- **Round 2 (§11, 2026-09-25)** does four things:
  - builds the sparse, tower-aware engine `f4_fp_tower`;
  - cross-checks it against `f4_fp` on the pilot's systems;
  - pre-registers larger sizes before running them, disclosing one system
    measured at a new size first (§11.3);
  - runs them.
- **The solving degree is not bounded.**
  - At `m = 2` it rises to 6 at `N = 20` and stays 6 at `N = 22`, in the
    Kummer and the isogeny family, at every prime and in the null. This refutes
    the pilot's Conjecture 3.
  - At `m = 3` it reaches 7 at `N = 15`.
  - Whether it grows linearly in `N` is not decided. `m = 4` ran out of
    memory at `N = 16`, where a diagnostic gives only `D ≥ 7`.
- **Round 3 (§12, 2026-09-25) is pre-registered.** It makes a memory-only change
  to the engine: the basis keeps the echelon's rows. It checks that the change
  reproduces round 2 exactly, then runs `m = 4` at `N = 16`. §12.3 discloses a
  profiling run of that system.

**Thread:** the prime regime of the index-calculus framework (`docs/ic/FRAMEWORK.md`).
**Siblings:** `RESEARCH_IC_BOUNDARY_LEDGER.md` (the table family and its law),
`RESEARCH_INDEX_CALCULUS_FACTOR_BASE.md` (structured `GF(p)` bases),
`research/ecdlp_autolab/yokoyama_lower_bound.md`,
`research/notes/cm-isogeny/RESEARCH_PKM_CRITERION.md`.

---

## 0. Summary

The prime regime of `ic bench` has two decomposition oracles, `subtract` and
`mitm`. Both are table methods. The boundary ledger derived their family law,
`S_family = 0.75·r^{1/6}` on a prime-order curve (ledger §11.3), so no base size
and no fold can make them cross rho. The scoreboard says what is left open: "an
algebraic decomposition oracle builds no such table and its asymptotics are the
question this ladder cannot reach." For binary and Koblitz curves the framework
has that oracle, `descent-algebraic`, with eleven engines behind it. For prime
fields it has none.

Petit, Kosters and Messeng (PKC 2016) published the one algebraic mechanism for
prime-field curves. Its factor base is `F = {P ∈ E(F_p) : x(P) ∈ V}`, where
`V = {x : L(x) = 0}` and `L = L_t ∘ … ∘ L_1` is a composition of low-degree maps
whose preimage tree splits completely over `F_p`. The decomposition system keeps
one auxiliary variable per level of the composition, so every membership
equation has low degree (Amadori–Pintore–Sala 2017, §2.2, system (6)). This is
the prime-field stand-in for the `F_2`-subspace of binary Weil descent.

This note designs that oracle for the framework:

- a factor-base builder, `prime-tower`, with three tower families: Kummer,
  Dickson/torus and 2-isogeny;
- a decomposition oracle, `pkm-tower`;
- an `F_p` counterpart of the `SystemSolver` plug point, with a reference engine
  and a measuring instrument;
- a staged plan whose first stage is a cheap test on the solver axis that decides
  whether the rest is worth building.

**The prediction, stated before measuring.** In §3 the design derives that the
tower presentation supplies the *field-equation* half of Weil descent (a finite
quotient ring with a square-free monomial basis) but not the *conjugate-equation*
half, which needs a Galois group that `F_p` does not have. A PKM system is
therefore a *one-generator* system over a `2^N`-point quotient, `N = m·t ≈ log₂ r`.
For such systems the XL (Macaulay) solving degree is provably at least
`d + δ`, where `δ` is the degree of the inverse of the summation polynomial on the
grid (§3.3), and `δ` is generically about `N`.

The pre-registered prediction is that F4's solving degree also grows linearly in
`N` (H0, §5.3). If it does, the tower family is closed on the solver axis with a
named mechanism, and Stage B is not built. If the solving degree stays bounded
(H1), prime-field index calculus becomes subexponential under a measured
heuristic, the exponent-moving outcome. Stage B then prices it end to end against
rho.

---

## 1. Why this oracle is the missing one

### 1.1 What the suite has for prime-field curves

| component | what it is | why it does not answer the question |
|:--|:--|:--|
| `ic bench` prime regime (`src/bin/ic/bench.rs::run_prime`) | `prime-abscissa` base with `subtract` or `mitm` | Table oracles only. Their family law is `Θ(r^{1/6})` above rho (ledger §11.3) |
| `ic prime` (`cryptanalysis` repo, PR #64) | automorphism-orbit bases, 2-point x-lookup | A table lookup, not an algebraic solve |
| `ec_index_calculus.rs` | Semaev `S₃`, small-x base, `m = 2` | `S₄` is computed but never used in a decomposition |
| `pkm_criterion.rs` | PKM *resistance audit*: Solinas weight, smoothness, arithmetic-progression "divisor sets" | Its own header says "Not a full PKM attack" |
| `f4_fp.rs`, `groebner_f4.rs` | `F_p` F4 (degree-bounded, `p < 2³²`) and Buchberger | Engines sized for 3–6 unknowns; not wired to any prime-field oracle |

### 1.2 What was tested elsewhere, and why it is not PKM

Every prime-field "structured base" measured so far used a **single univariate
membership equation** per coordinate:

- the multiplicative-subgroup base with constraint `X^d − 1`
  (`RESEARCH_INDEX_CALCULUS_FACTOR_BASE.md` §1.1, 2.4× slower at 19 bits);
- the autoresearcher's toy suite `crypto_autoresearcher.index_calculus`. Its
  `subgroup` base uses membership `x^d − g^d`; its draft msolve engine (PR #1417
  there) solves that system and refuses `p ≥ 2^16` because the packaged msolve
  0.6.5 build is wrong above it;
- the autoresearcher's `EXP-ALPF-005` "Kummer / rational-map" cell. It
  considered the parameterised version and set it aside ("S4 becomes degree 24 in
  t → not tractable"), then measured a product of two quadratics instead: again
  one univariate equation.

These are all the *naive* presentation `⟨S_{m+1}, F_V(x_1), …, F_V(x_m)⟩`.
Yokoyama, Yasuda, Takahashi and Kogure (JMC 2020) bound exactly that presentation
from below and show it cannot beat generic methods. The autoresearcher's own
literature review records the gap
(`experiments/EXP-ECDLP-SOURCE-TAG-JOIN-001/literature-refresh-2026-07-17.md`):
its instruments "do not compile the recursive equations of `L` … The verified
source-tag negative therefore does not test the central algebraic mechanism of
the rational-map proposal." Its PKM arm, `ARM-PKM-PRIME` (TASK-20260725-661,
GOAL-RELN-001), has been blocked since 2026-07-25 on exactly the missing piece:
"a literature-derived fixture binding: … factor-base rule, decomposition map,
and polynomial system." §2 supplies that binding.

### 1.3 Two scope corrections this design depends on

- **The "PKM divisor-set construction".** `RESEARCH_INDEX_CALCULUS_FACTOR_BASE.md`
  §1.1 describes it as a point-coordinate interval base, and
  `RESEARCH_PKM_CRITERION.md` audits arithmetic-progression divisor sets. Neither
  is PKM's construction. PKM's `V` is the zero set of a *composed* map, with one
  auxiliary variable per level. No interval, progression or single-polynomial
  base tests it.
- **Proposition 1 of that note (sparsity invariance)** fixes the set `V` and
  re-encodes the generator of `(F_V)` *in the same ring*. The tower introduces new
  variables, so the PKM ideal lives in a larger ring, and only its elimination
  ideal equals the naive one. Proposition 1 does not apply to it. As far as the
  repository's summary of Yokoyama et al. goes, neither does their lower bound,
  which is stated for the naive presentation. (The paper's scope should be
  checked against its full text; `yokoyama_lower_bound.md` is a paraphrase.)

The gap is therefore exactly one mechanism: the auxiliary-variable tower. §3
explains why that mechanism is plausible, and why it is probably not enough.

---

## 2. The construction

Notation:
- `E : y² = x³ + a x + b` over `F_p`, with a prime-order subgroup of order `r`
  (the framework's instances have cofactor 1);
- `m` summands, `t` tower levels, and `N = m·t` tower variables;
- `V ⊂ F_p` with `|V| = 2^t`, and `F = {P ∈ E(F_p) : x(P) ∈ V}`;
- about half of `V` are abscissae of `E`, so the base has about `2^{t−1}` columns
  and `2^t` signed points.

### 2.1 Tower factor bases

A **tower** is a list of degree-2 maps `L_1, …, L_t` and a constant `c` such that
the full preimage `V = (L_t ∘ … ∘ L_1)^{−1}(c)` has `2^t` distinct elements, all in
`F_p`. Membership is written with auxiliary variables `y_0 = x, y_1, …, y_{t−1}`:

```text
    T_j :  y_{j+1} − L_{j+1}(y_j) = 0      (j = 0 … t−2, cleared of denominators)
    T_t :  L_t(y_{t−1}) − c       = 0
```

Each `T_j` has degree 2. Three families give a completely split preimage tree.

| kind | `V` | step map | condition on `p` | construction cost |
|:--|:--|:--|:--|:--|
| `kummer` | `g·μ_{2^t}`, a coset of the `2^t`-th roots of unity | `y ↦ y²`, with `y_0 = x/g` and `y_{t−1}² = 1` | `2^t \| p − 1` | one exponentiation for `ζ`, then `2^t` multiplications |
| `dickson` | `{w + w^{−1} : w ∈ u·μ_{2^t}}`, with `u ∈ F_p^*` or `u` in the norm-1 torus `T ⊂ F_{p²}^*` | `y ↦ y² − 2` (Dickson `D_2`), final `y_{t−1}² − 2 = c`, `c = u^{2^t} + u^{−2^t}` | `2^t \| p − 1` (split case) or `2^t \| p + 1` (torus case) | `F_{p²}` arithmetic for `u` |
| `isogeny` | `x(R' + ⟨Q⟩)` on an auxiliary curve `E'/F_p` with a rational point `Q` of order `2^t` | the x-map of a 2-isogeny (Vélu): `y_{j+1}(y_j − ξ_j) = y_j² − ξ_j y_j + τ_j` | none: `E'` is searched for | an auxiliary-curve search, charged (below) |

Conditions on the choice of coset or point, so that `|V| = 2^t` and the tower is
unramified:
- `kummer`: any `g ∈ F_p^*` works, drawn from a public seed.
- `dickson`: require `u² ∉ μ_{2^t}`, which makes the coset disjoint from its
  inverse; also `y_j ≠ 0` on `V`, i.e. `c` avoids the critical value `−2` along
  the chain.
- `isogeny`: require `2R' ∉ ⟨Q⟩`, so the coset and its negative are disjoint and
  `x` is injective on it.

**The isogeny tower in full.** Let `Q_0 = Q` on `E'_0 = E'`. For
`j = 0, …, t−1`:
- the kernel point is `K_j = 2^{t−1−j}·Q_j`, with `ξ_j = x(K_j)`;
- `τ_j = 3ξ_j² + a'_j`;
- the step is the Vélu map `x ↦ x + τ_j/(x − ξ_j)`;
- the codomain has `a'_{j+1} = a'_j − 5τ_j` and `b'_{j+1} = b'_j − 7ξ_jτ_j`;
- the next point is `Q_{j+1} = φ_j(Q_j)`.

The last step's value is `c = x(Φ(R'))`, where `Φ = φ_{t−1} ∘ … ∘ φ_0`. Then
`L(x) = c` has exactly the `2^t` roots `x(R' + ⟨Q⟩)`, all rational.

This is the family that needs no condition on `p`. An auxiliary curve with a
cyclic rational subgroup of order `2^t` is found by trying random curves until
`2^t | #E'(F_p)` (about `2^t` point counts, polynomial per count with SEA,
`O(p^{1/4})` with the repository's baby-step giant-step at toy sizes). This is
PKM's "pre-computation specific for every curve". **It is charged to the factor
base's phase, never free**: at toy sizes it can dominate, and the report shows it
as its own line.

Mixed-radix towers (steps of degree 3 or 5 where `p ± 1` has those factors) are a
parameter (`degrees=2,2,3`). They are the case the Nikolaev criterion flags on
smooth `p ± 1` primes. Stage A covers them only at the end.

### 2.2 The decomposition system

A target `R` with `x_R = x(R)` gives, for `m` summands with tower variables
`y_{i,j}` (`i = 1…m`, `j = 0…t−1`) and `x_i = g·y_{i,0}` (Kummer) or `x_i = y_{i,0}`:

- **`full`** (`m ≤ 3`): `S_{m+1}(x_1, …, x_m, x_R) = 0` together with all `m·t`
  tower equations.
  - `S₃` is the closed form in `ec_index_calculus.rs`.
  - `S₄ = Res_U(S₃(x_1,x_2,U), S₃(x_3,x_R,U))`, by the two-quadratic formula
    `Res = (AC′ − A′C)² − (AB′ − A′B)(BC′ − B′C)`. No general resultant code is
    needed.
- **`chain`** (Semaev 2015; any `m ≥ 3`):
  `S₃(x_1,x_2,u_1) = 0`, `S₃(u_k, x_{k+2}, u_{k+1}) = 0`, `S₃(u_{m−2}, x_m, x_R) = 0`,
  with free `F_p` unknowns `u_k` (abscissae of partial sums), plus the towers.
- **`tree`** (Karabina's topology; `m = 4`): `S₃(x_1,x_2,u_1)`, `S₃(x_3,x_4,u_2)`,
  `S₃(u_1,u_2,x_R)`. The autoresearcher's `KN-TECH-b18366` names the chain/tree
  comparison as an open measurement.

**Presentation** is a parameter:
- `raw` hands PKM's system (6) to the engine as written;
- `reduced` first replaces `S` by its normal form modulo the towers (§2.4).

The two generate the same ideal. Engines may behave differently on them, and
that difference is data.

**Symmetry** is a parameter:
- `none` keeps the `m!` orderings;
- `disjoint` gives each summand its own coset `V_i` (Galbraith–Gebregiyorgis,
  Amadori–Pintore–Sala). That breaks the permutation symmetry at the price of an
  `m`-times larger base.

### 2.3 Sizing

A random target decomposes with probability about `|V|^m/(m!·r)`. The balanced
choice of `|V|` against the linear algebra is derived in §3.7. For the solver-axis
probe (Stage A) the natural choice is `|V|^m ≈ m!·r`, about one planted solution
per target, which gives `N = m·t ≈ log₂ r + log₂ m!`.

### 2.4 The quotient ring: what the tower is, algebraically

**Lemma 1.** Order the variables of each block `y_{i,0} ≻ y_{i,1} ≻ … ≻ y_{i,t−1}`
under grevlex. Then:
1. every tower equation has leading monomial `y_{i,j}²`;
2. the `m·t` tower equations form a Gröbner basis, since their leading monomials
   are pairwise coprime;
3. the standard monomials are the `2^N` square-free monomials;
4. the quotient `R = F_p[y]/(towers)` is isomorphic to `F_p^{V^m}` by the Chinese
   Remainder Theorem, because the towers are radical with `2^N` rational zeros.

For the isogeny step, part 1 holds because `y_{j+1} ≺ y_j` makes
`y_j² ≻ y_j·y_{j+1}` in grevlex.

So the tower does for `F_p` what the field equations `v² = v` do for `F_2`: it
turns a degree-`2^t` univariate constraint into `t` quadrics over a Boolean-like
basis. This is the half of Weil descent PKM transplants.

For the **Kummer** tower the correspondence is exact:
- `R` is the group algebra `F_p[(Z/2^t)^m]`;
- the square-free monomial `∏_j y_{i,j}^{b_j}` is the character `k ↦ ζ^{k·B}`,
  where `B = Σ_j b_j 2^j`;
- `x_i^e = g^e·∏_{j ∈ bits(e)} y_{i,j}`, a single monomial of degree `popcount(e)`.

The binary descent has the same property: `x^{2^j}` is `F_2`-linear, so `x^e` has
degree `popcount(e)` in the Boolean variables. The difference is §3.2.

---

## 3. What the tower can and cannot buy

### 3.1 What it buys

In the naive presentation the membership equation has degree `|V| = 2^t`, and the
regularity Yokoyama et al. use grows like `m·|V|`, exponential in `t`. In the
tower presentation every equation has degree 2 or `deg S`, and the Macaulay bound
on the regularity is `N + d`, linear in `t`. If the solving degree tracked a small
constant instead, as the first-fall-degree assumption hopes for binary Weil
descent, the per-target cost would be polynomial in `log p`. Choosing `m` near
`√(log p)` would then make prime-field index calculus subexponential. That is the
exponent-moving outcome, and it is why the mechanism deserves a real test.

### 3.2 What it does not buy: the conjugate equations

Binary Weil descent turns **one** equation over `F_{2^n}` into **n** equations
over `F_2`. The `n` descended equations are the coordinates of the equation's
Galois conjugates, and they are valid because the unknowns are `F_2`-rational.
With `n ≈ m·ℓ` generators in `m·ℓ` Boolean unknowns, the system is square, and
the semi-regular degree is a small fraction of `N`.

`F_p` has no subfield and no Galois group to conjugate by
(`RESEARCH_INDEX_CALCULUS_FACTOR_BASE.md` §2.3). The tower supplies the finite
quotient of §2.4, but the summation polynomial stays **one** equation:
- the `full` presentation is `1 + N` equations in `N` unknowns, of which `N` are
  the "field equations";
- the `chain` presentation is `m − 1` summation equations and `m − 2` free `F_p`
  unknowns, so still one net equation.

That single `F_p` equation cuts `2^N ≈ p` grid points to about one, so it carries
`log₂ p` bits, where the binary system spreads the same information over `n`
equations. The autoresearcher already has a matching hypothesis: `KN-OPEN-d6ad3f`
calls it `H_count`, the mechanism needing "a descent block of many generators",
and notes that "two generators return deficit 0 at `p = 2` as well". A PKM system
is the one-generator case in odd characteristic.

### 3.3 A bound for one generator

Let `S̄ = NF(S)` be the normal form of the summation polynomial modulo the
towers, so `S̄ ∈ R`. Let `d = deg S` in the presentation handed to the engine, and
let `Z ⊂ V^m` be its zero set on the grid, with `σ = |Z|` (typically 0, 1 or a few
orderings of one decomposition). Let `I_(D)` be the degree-`D` Macaulay space of
the generators: every combination `Σ h_k f_k` with `deg(h_k f_k) ≤ D`.

**Proposition 2.** For every `D`:
1. `NF(I_(D)) = S̄·R_{≤D−d}`, the products reduced in `R`, where `R_{≤e}` is the
   span of square-free monomials of degree at most `e`. Hence
   `dim NF(I_(D)) ≤ Σ_{k ≤ D−d} C(N, k)`.
2. A linear form `ℓ` vanishing on `Z` lies in `I_(D)` only if some function on
   `V^m` that agrees with `ℓ/S̄` off `Z` has degree at most `D − d` in the
   square-free basis. Write `δ_Z(ℓ/S̄)` for the least such degree. Suppose the
   reduced Gröbner basis contains a linear form `ℓ`. This always holds for
   `σ = 1`; for `σ = 0` read `ℓ = 1`, the refutation. Then the Macaulay (XL)
   solving degree is at least `d + δ_Z(ℓ/S̄)`.

*Proof sketch.*
- Tower multiples reduce to 0, and `NF(h·S) = NF(NF(h)·S̄)` with
  `deg NF(h) ≤ deg h`. That gives part 1.
- `R ≅ F_p^{V^m}` turns products into pointwise products. So `ℓ = NF(h·S)` forces
  `NF(h) = ℓ/S̄` at every point off `Z`, and leaves it free on `Z`. That gives
  part 2. ∎

`δ_Z` is computable exactly on toy grids. For the Kummer tower it is a
`(Z/2^t)^m` discrete Fourier transform followed by a popcount weight, `O(2^N·N)`.
A function on `2^N` points with `σ` free values generically has degree `N − O(1)`:
lowering it to `e` imposes `Σ_{k>e} C(N,k)` linear conditions on `σ` unknowns.
**So unless `1/S̄` is unusually structured on the grid, the XL solving degree of a
PKM system is about `N + d`.** For the naive presentation the same argument gives
`m·|V|`, the regularity Yokoyama et al. use. The tower moves the exponent of the
*degree* from `2^t` to `t`; it does not make the degree small.

### 3.4 Kummer and Dickson towers are degenerate at infinity

Take the top-degree parts of the `raw` system. The towers contribute `y_{i,j}²`.
`S_{m+1}` has degree `2^{m−1}` in each `x_i`. So a monomial of the largest possible
total degree, `m·2^{m−1}`, has every exponent equal to `2^{m−1} ≥ 2`: the top form
is `c(x_R)·∏ x_i^{2^{m−1}}`, which lies in `(y_{i,0}²)` whenever `c(x_R) ≠ 0`. For
`S₃` the top form is `x₁²x₂²`, with `c = 1`. Stage 0 checks `c` for `m ≥ 3`.

The top-part ideal is therefore just `(y²)`, and its Hilbert series is
`(1+t)^N`, which stays positive through degree `N`. The Bardet–Faugère–Salvy
degree of regularity is `N + 1` whatever `S` is: the top parts carry no
information.

The `reduced` presentation does not escape this. `x_i^e` with `e ≤ 2^{m−1}` touches
only levels `j ≤ m − 1`, so the top part of `S̄` lives in the lowest `m` levels.
The higher levels are again free at infinity, and the regularity is still about
`N`.

Two consequences follow:
- **The semi-regular prediction is the wrong yardstick here.** The generic
  `F_p` series for these degrees is `(1+t)^N·(1 − t^d)`, whose first non-positive
  coefficient is at `⌈(N + d)/2⌉`. The report prints it, but it is not a bound.
- **Isogeny towers are not degenerate in this way.** In their top-part algebra,
  `y_j² ≡ y_j·y_{j+1}`, so `top(S)` survives. The three families may therefore
  scale differently, and Stage A measures all three.

The framework's `semi_regular_degree` is the Boolean formula; the `F_p` engines
need the `F_p` series (§4.4).

### 3.5 What F4 can still do

Proposition 2 bounds the Macaulay (XL) computation. F4 can do better. When a row
reduces to a lower degree (a "fall"), F4 keeps the result and multiplies it by
monomials up to the degree bound. Those products lie outside `S̄·R_{≤D−d}`.

In the Kummer algebra falls are everywhere, because `y_j·y_j = y_{j+1}` is a
carry. Whether those carries cascade into the linear forms below degree `N` is not
something a counting argument decides. The invariant that bounds F4 is the last
fall degree (Huang–Kosters–Yeo). **This is the one question Stage A exists to
answer, and the reason the oracle is worth building at all.**

### 3.6 What each outcome costs

If F4's solving degree grows as `β·N`, the matrix at that degree has width
`Σ_{k ≤ βN} C(N, k) = 2^{H(β)·N + o(N)}`, where `H` is the binary entropy function.
Since `N ≈ log₂ p`, that width is `p^{H(β)+o(1)}`, and eliminating the matrix costs
a power of its width. The per-target cost is then a power of `p`, not of `log p`.

- **`β ≥ 1/2`.** The width alone is at least `2^{N−1}`, and `2^N ≥ |V|^m`. Every
  target then costs at least about `√r`, which is rho's whole budget, before a
  single relation is kept.
- **Smaller `β`.** The per-target exponent races `m` in §3.7's balance, which
  becomes about `r^{2(1+γ)/(m+1)}` for a per-target cost of `r^γ`.

The decision rule of §5.2 is set so that only a bounded degree (`β ≈ 0`) reaches
Stage B. If the solving degree is bounded, or grows as `o(N/log N)`, the cost per
target is polynomial in `log p`. That is the only regime in which the rest of this
design pays.

### 3.7 The linear-algebra floor: only `m ≥ 4` can matter

With per-target cost `T`, a base of `|V|` points and sparse linear algebra, the
end-to-end cost is roughly `m!·r·T / |V|^{m−1} + m·|V|²`. The optimum is
`|V| ≈ (m!·r·T)^{1/(m+1)}`, giving a total of about `r^{2/(m+1)}` for `T = r^{o(1)}`.
- For `m ≤ 2` this is above `r^{1/2}`, and for `m = 3` it is `r^{1/2}` times
  factors that only grow. No oracle, however fast, beats rho asymptotically at
  `m ≤ 3`; `m = 3` at best ties rho's exponent.
- For `m = 4` it is `r^{2/5}`, and for `m = 5`, `r^{1/3}`.

A large-prime variation shifts this balance; it is a Stage B lever, not a Stage A
question.

So the decisive Stage A cells are `m ≥ 4`, in the `chain` or `tree` presentation.
Cells at `m = 2, 3` measure the mechanism, not the verdict.

---

## 4. Framework integration

### 4.1 New types (in `ic_framework/stages.rs`)

```rust
/// A polynomial system over F_p, the currency between a prime-field
/// algebraic oracle and an F_p engine.  Polynomials use f4_fp::Poly.
pub struct FpSystem {
    pub p: u64,
    pub n_vars: usize,
    pub equations: Vec<f4_fp::Poly>,
    pub var_names: Vec<String>,
    /// The tower equations' indices, so an engine may treat them as field
    /// equations (e.g. reduce into the square-free basis) instead of generators.
    pub tower_equations: Vec<usize>,
}

pub enum FpVerdict { Solved(Vec<Vec<u64>>), Unsatisfiable, BudgetExceeded }

pub trait FpSystemSolver: Send + Sync {
    fn name(&self) -> &str;
    fn describe(&self) -> String;
    fn parameters(&self) -> &[(&str, &str)] { &[] }
    fn accepts(&self, shape: &SystemShape) -> bool { let _ = shape; true }
    fn finds_every_solution(&self) -> bool { true }
    fn solve(&self, system: &FpSystem, params: &Params, budget: Option<Duration>)
        -> (FpVerdict, SolverCost);
}
```

`SolverCost`, `SolverTotals` and `SystemShape` are reused unchanged. The runner
already prices `solver_totals()` into `S` (§4.6). The contract of `SystemSolver`
applies word for word:
- count something;
- respect the budget;
- never guess `Unsatisfiable`;
- report `solving_degree`, not `degree_reached`;
- decline what you cannot do.

A parallel trait is preferred to a generic one. `SolverVerdict::Solved` carries
bitmasks, and every existing engine would otherwise have to change.

### 4.2 `prime-tower` (a `FactorBaseBuilder<PrimeCurve>` in `plugins.rs`)

| parameter | meaning |
|:--|:--|
| `kind` | `kummer`, `dickson` or `isogeny` |
| `levels` | `t`; `\|V\| = 2^t` (or the product of `degrees`) |
| `degrees` | optional mixed radix, e.g. `2,2,3` |
| `seed` | public seed for the coset or auxiliary point; never target-dependent |
| `symmetry` | `none` or `disjoint` (then `m` cosets) |
| `aux_budget` | `isogeny` only: maximum auxiliary curves tried |

**Charges.** Every field operation is charged to `FactorBase::cost` as a counter,
priced by the calibration:
- construction of the tower, including the auxiliary curve search and its point
  counts;
- the enumeration of `V`;
- the Legendre symbol and square root per element, exactly as
  `prime_factor_base` charges them.

**Refusals.** It refuses with a reason, rather than falling back, when the
instance's `p` lacks the structure the kind needs (`kummer` needs `2^t | p − 1`),
and suggests `isogeny` or a filtered instance. A fallback would change the
experiment.

**Output.** The column map is `ColumnFold::Abscissa`, one column per abscissa with
`−P` at coefficient `r − 1`, the same as `prime-abscissa`. The table oracles
therefore run on the same base unchanged, and `mitm` on `prime-tower` is the
cell's reference row. `FactorBase` gains
`pub tower: Option<TowerSpec>`, next to `subspace_basis`, which `pkm-tower` reads.

### 4.3 `pkm-tower` (a `DecompositionOracle<PrimeCurve>`)

Modelled line for line on `DescentAlgebraicOracle`, with an `FpSolverHarness`
mirroring `SolverHarness`.

| parameter | meaning |
|:--|:--|
| `m` | summands, 2–5 |
| `mode` | `full` (`m ≤ 3`), `chain` or `tree` |
| `presentation` | `raw` or `reduced` |
| solver | any `FpSystemSolver`, via `--solver` |

**`prepare`**
- read `fb.tower`, or refuse (the base is not a tower base);
- build one system to learn its shape;
- let the engine decline the shape once, with a reason a sweep can report.

**`decompose(R)`**
1. Build the system for `x_R`.
2. Solve it within the budget.
3. For each solution, map `y_{i,0} → x_i`. A value of `x_i` that is not an
   abscissa of `E` (it names a point on the twist) counts as a lift failure.
4. Call `lift_abscissae`, which charges group additions.
5. Return the first sum equal to `R`.

**Counting rules.**
- A system solved with no liftable solution counts as `unliftable_systems`.
- `BudgetExceeded` is counted and never folded into "does not decompose".
- `last_system` and `solver_totals` behave as in `descent-algebraic`.

**Cross-check.** Every target of a full run is also decided by `tower-exhaustive`
(§4.4) and the two answers must agree, as `AGENTS.md` §6 requires of a new oracle.

### 4.4 `F_p` engines (a new `ic_framework/fp_solvers.rs`, with `fp_solver_registry()`)

| name | what it is | unit | role |
|:--|:--|:--|:--|
| `tower-exhaustive` | For `m = 2`, each `v ∈ V` gives `S₃(v, X, x_R)`, a quadratic; its roots are tested for membership in `V` by hash. `chain` enumerates `V^{m−1}` and roots the rest | square roots | **the reference**: complete, and the cost to beat on the same system |
| `f4-fp` | adapter to `f4_fp::solve`: grevlex, a degree bound defaulting to `N + d + 2`, a deadline | `F_p` multiply-adds (elimination only) | first real engine |
| `macaulay-fp` | Macaulay matrices in `R` (square-free columns, rows `u·S̄` reduced), `D = d, d+1, …` | `F_p` multiply-adds | **the instrument**: Hilbert function, first-fall degree, the degree at which linear forms appear, and `δ_Z` on toy grids. Tests Proposition 2 |
| `f4-fp-tower` | *Stage B only.* F4 whose monomials are `u64` square-free masks, multiplied with the tower's rewriting (carry-add for `kummer`, the Chebyshev product rule for `dickson`, a general rewrite for `isogeny`), with sparse `F_p` elimination and Gebauer–Möller; built on the architecture of `pq_f4_f2.rs` | `F_p` multiply-adds | the engine that could exploit falls at scale |
| `msolve` | optional external engine | wall time (`measured`) | only after a known-answer battery passes at the run's `p`; msolve 0.6.5 fails it above `2^16` |

Changes needed in `f4_fp.rs`, which are Stage 0 work:
1. A multiply-add counter in `eliminate`, returned in `F4Report`.
2. Replace `roots_univariate`'s `O(p)` scan with `gcd(f, x^p − x)` and
   equal-degree splitting. Keep the scan as a cross-check for `p < 2^16`.
3. `semi_regular_degree_fp(n_vars, degrees)`: the first non-positive coefficient
   of `Π(1 − t^{d_i})/(1 − t)^{n_vars}`.

### 4.5 Instances

`find_prime_order_curve` draws `p` freely. Two things are added:
- `find_prime_order_curve_where(bits, seed, pred)` with predicates
  `v2(p−1) ≥ t` and `v2(p+1) ≥ t`;
- the sweep key `"p_filter": "v2m1>=6"`.

**A filtered instance is a restricted prime family, and every report on one says
so.** Conclusions about generic primes come only from `kind = isogeny`, which runs
on unfiltered instances. Rho and the reference oracles always run on the same
curve, so every ratio is matched.

### 4.6 Accounting

- Everything is inside `S` (`AGENTS.md` §2): tower construction, the auxiliary
  search, the base, target generation, building each system, every solver call
  (refutations included, and they are the common case), lifting, the linear
  algebra and verification.
- The new unit `F_p multiply-adds` gets a calibration factor `ns_per_fp_muladd`,
  measured per process.
- The runner's pricing arm (`ic_framework/mod.rs`, next to `word XORs`) prices it
  by count once `docs/ic/calibration.json` carries a pinned ratio, which is added
  the way §12 of the ledger note pinned the others. Until then, rows are priced
  `measured` by wall time and say so.
- A count that covers only the elimination carries the qualifier
  `(elimination only)`, as the F4-family `F_2` engines do.

---

## 5. Boundaries and the falsification target (pre-registered)

### 5.1 Boundaries

- **Floor:** `S ≥ √(π/2A)`, which is `0.886` for `A = 2` on a generic curve.
- **Reference:** counted Pollard rho on the same instance (`rho_reference`, the
  `vs rho` column).
- **Cell references, on the same `prime-tower` base:**
  - `mitm` is the strongest table method;
  - `tower-exhaustive` is the strongest enumeration of the same system.

  An algebraic engine earns a cell only by beating both, per `FRAMEWORK.md` §4.
- **Family floor for this oracle (§3.7):** `Θ(r^{2/(m+1)})` for any per-target
  cost `r^{o(1)}`. For `m ≤ 3` it is at or above rho.

### 5.2 Stage A: the solver axis, and the only measurement that decides

**Cells.**
- `kind ∈ {kummer, dickson, isogeny}`;
- `m ∈ {2, 3, 4, 5}`: `full` for `m ≤ 3`; `chain` for `m ≥ 3`; `tree` for `m = 4`;
- `presentation ∈ {raw, reduced}`;
- `N = m·t` over at least four consecutive values of `t` per `m`, as far as the
  engine goes within a 600 s budget and an 8 GiB cap per call;
- `p` sized so that `|V|^m ≈ m!·r` (§2.3).

**Targets.** 8 planted-decomposable and 8 random targets per cell, from public
seeds.

**Statistics.**
- `D_F4`: `solving_degree` from `f4-fp`, or from `f4-fp-tower` if it exists by
  then. A run whose `pairs_above_bound > 0` at the end reports a lower bound, never
  a value.
- `D_ref`: the refutation degree on non-decomposable targets.
- `D_XL` and `δ_Z` from `macaulay-fp`.
- The largest matrix, as a cost proxy.

**Fit.** Per (`kind`, `m`, `mode`, `presentation`), fit `D = α + β·N` by least
squares. The 95% confidence interval for `β` comes from bootstrapping over targets
within each `N`.

**Hypotheses and decision rule.**

| outcome | condition on `β` | decision |
|:--|:--|:--|
| **H0**, the one-generator law | lower 95% bound above `0.25` | close the tower family on the solver axis for that kind and presentation, with the mechanism of §3.2 recorded. Stage B is not built for it |
| **H1**, PKM's hope | upper 95% bound below `0.10`, and `D_F4` below the shape-matched null (§5.3) at the largest `N` | go to Stage B |
| inconclusive | anything else | extend `N` (`f4-fp-tower`, more memory) before deciding |

Pre-registered expectation: **H0 for all three kinds**, with `β` near `1` for
`kummer` and `dickson` (§3.4) and possibly smaller for `isogeny`.

**Theorem check.** Proposition 2 predicts `D_XL ≥ d + δ_Z` on every system. A
violation is a bug in the instrument or the proof, and is reported as one before
anything else is read.

### 5.3 Controls

Every control runs at the same (`kind`, `N`, `m`) as the cell it controls.

| control | what changes | what it rules out |
|:--|:--|:--|
| **naive** | the same `V`, membership as one univariate `F_V(x_i)`, no auxiliaries | Measures what the tower buys over the Yokoyama regime (§3.1) |
| **shape-matched null** | `S̄` replaced by a random element of `R` with the same monomial support and a planted zero | Whether any low degree comes from the summation polynomial's structure or only from the shape |
| **random set** | `V` a random `2^t`-subset of `F_p`, membership interpolated | Whether the tower's algebraic structure matters, as opposed to any set of that size |
| **generator-count ladder** (instrument power) | `g − 1` extra random quadrics vanishing at the planted solution, `g = 1, 2, 4, …, N/2` | This is `H_count` of `KN-OPEN-d6ad3f` tested in odd characteristic. The solving degree must fall as `g` grows; **if the instrument cannot see that fall, it has no power and any negative from it is void** |
| **binary twin** | `descent-algebraic` on a random binary curve with matched `m` and `N = m·ℓ` | The characteristic-2 comparison, from the existing `ic descent --solver` tables |

### 5.4 Inadmissible

The following void a row or a comparison:
- changing the tower, the coset or the base per target;
- counting `BudgetExceeded` as "does not decompose";
- leaving refutations out of a cost;
- charging the auxiliary-curve search as free;
- choosing seeds after seeing results;
- quoting a per-call solver cost as a speed (`AGENTS.md` §2: Stage A is a stage
  diagnostic end to end);
- comparing a restricted-`p` row to a generic-`p` row without saying so.

---

## 6. Staged plan

| stage | builds | produces | proceeds when |
|:--|:--|:--|:--|
| **0: construction** | `prime_tower.rs` (the three kinds, enumeration, normal form), `pkm_system.rs` (full, chain, tree, raw, reduced, `S₄` by resultant), the `FpSystemSolver` trait, `tower-exhaustive`, the `f4-fp` adapter with its counter and root finder, instance filters | tests only (§7) | every correctness test passes |
| **A: solver axis** | `macaulay-fp`, the Stage A cells and controls, as `ic descent --field prime` (paired engines, fingerprints, as the binary `ic descent --solver`) | `docs/ic/runs/pkm-tower-degrees-YYYY-MM-DD.json` and a round in this note | H0 closes a family; H1 opens Stage B |
| **B: end to end** *(only after H1)* | `prime-tower` and `pkm-tower` in `ic bench`, `f4-fp-tower`, and a pinned `ns_per_fp_muladd` | `S` against rho and against `mitm` on the same base, and exponent fits over at least four sizes at `m ≥ 4` | an `S` row on the scoreboard, classified by `AGENTS.md` §3 |

Stage B also needs `AGENTS.md` §8: the frozen regression suite does not apply to
prime-field systems. An equivalent matched suite is frozen under the parent
accounting contract, with the reason recorded.

---

## 7. Files and tests

**New files**
- `src/cryptanalysis/prime_tower.rs`: `TowerSpec`, the three constructions,
  `enumerate_v`, `tower_equations`, `normal_form`.
- `src/cryptanalysis/pkm_system.rs`: system builders and the lift map.
- `src/cryptanalysis/ic_framework/fp_solvers.rs`: engines and registry.

**Edits**
- `ic_framework/stages.rs`: the `F_p` types.
- `ic_framework/plugins.rs`: `PrimeTowerBase`, `PkmTowerOracle`, `FpSolverHarness`.
- `ic_boundary.rs`: `FactorBase::tower`, the instance filter, `ns_per_fp_muladd`.
- `ic_framework/mod.rs`: the pricing arm.
- `f4_fp.rs`: the counter, root finder and `F_p` semi-regular series.
- `src/bin/ic/bench.rs` and `src/bin/ic/descent.rs`: CLI wiring.
- `docs/ic/sweeps/pkm-tower-probe.json`.

**Tests (all at toy `p`)**

1. **Towers**: each kind yields `2^t` distinct elements of `F_p`; every tower
   equation vanishes on `V`; the leading monomials are `y²` in grevlex; the
   standard monomial count is `2^t`; `isogeny` codomains satisfy Vélu, checked by
   mapping points.
2. **Systems**: a planted decomposition is a zero of every presentation; `raw` and
   `reduced` have the same zero set on `V^m`; `S₄` by resultant vanishes at planted
   triples; `chain` and `full` agree on `m = 3`.
3. **Engines**: `f4-fp` and `tower-exhaustive` return the same solution sets on
   every target of a sweep; an unsatisfiable system is never `Solved`; the budget
   is respected.
4. **Oracle**: `pkm-tower` agrees with `tower-exhaustive` on every target of a full
   run; lift failures and unliftable systems are counted, not dropped; a recovered
   logarithm is verified by `[k]G = Q`.
5. **Instrument**: `macaulay-fp` reproduces Lemma 1's Hilbert function
   `Σ C(N,k)` on the towers alone; Proposition 2's inequality holds on every probe
   system; `δ_Z` from the Fourier route equals `δ_Z` from linear algebra on
   `N ≤ 10`.

---

## 8. What this unblocks elsewhere (pointers only)

In the autoresearcher repository (`aburan28/crypto-autoresearcher`):
- **`ARM-PKM-PRIME`** (GOAL-RELN-001): §2 is the literature-derived fixture binding
  its materialisation gate asks for.
- **`KN-OPEN-020`**: Stage A is a measurement on the solver axis, where
  `IDEA-20260902-701458` says the universal no-go must be decided.
- **`KN-OPEN-002`**: the growth of the solving degree over prime fields, in the
  presentation the literature proposes.
- **`KN-OPEN-d6ad3f`**: the generator-count ladder of §5.3 separates characteristic
  from generator count in odd characteristic, the confound its digit twin could
  not break.
- **`IDEA-20260921-ffd330`**: its power-residue cosets are the Kummer family on the
  *yield* axis. `KN-FIND-007` (mean yield depends on `|V|` only) is why this design
  tests the *solver* axis instead.

A message carries a pointer, not a permission. None of these records changes
because this note exists.

---

## 9. What this note does not claim

- It does not claim PKM's construction as published uses exactly these three tower
  families. The description of PKM's system follows Amadori–Pintore–Sala §2.2; the
  PKM paper itself is not in this container. The three families are the natural
  instances of "a composition of low-degree maps whose preimage tree splits
  completely", and Stage 0 should check the choice against the paper.
- It does not claim H0. §3 derives a bound for XL (Proposition 2) and a heuristic
  for F4. Stage A measures the heuristic.
- Nothing here is a statement about a deployed curve. Every `N` in Stage A is
  toy-sized, and conclusions are scoped to the cells measured.

---

## 10. Pilot, 2026-09-24: the pre-registered expectation fails at toy sizes

This section was written after §§0–9 were committed and after every run it
reports. It is a pilot, not Stage A. It runs a small part of Stage A's cells
with the repository's dense `f4_fp`, to test the prediction the staged plan
rests on before anything else is built. It re-grades nothing:
- §10.3 applies the rule of §5.2 as written;
- §10.8 proposes an amendment for Stage A and marks it post hoc.

The systems come from `examples/pkm_tower_pilot.rs`. The raw rows, their exact
flags, and the two scripts that tabulate and cross-check them are in
`research/pkm_tower_pilot_20260924/`.

### 10.1 What was run

| item | pilot | §5.2 asks for |
|:--|:--|:--|
| kinds | `kummer`, `dickson` (split case only), `isogeny` | the same, `dickson` in both cases |
| `m`, mode | 2, `full`: `S₃` and the `2t` tower equations. 3, `chain`: `S₃(x₁,x₂,u)` and `S₃(u,x₃,x_R)`, with a free `u`. 4, `chain`: three `S₃` links with two free unknowns | 2–5; `full`, `chain`, `tree` |
| presentation | `raw` only | `raw` and `reduced` |
| `N` | `m = 2`: 2–18 (20 attempted); `m = 3`: 6–12 (15 attempted); `m = 4`: 8, 12 | at least four consecutive `t` per `m` |
| prime | `p = 786433 = 3·2^18 + 1` in every cell, so `2^t \| p − 1` up to `t = 18`. Two checks at `p = 3221225473 = 3·2^30 + 1` | sized per `N`, so that `\|V\|^m ≈ m!·r` |
| curve | a random curve over `F_p` per cell; the solver axis needs no prime order | an instance of prime order |
| targets | 2 planted and 2 random per cell up to `t = 7`, random only above | 8 and 8 |
| engine | `f4_fp::f4`: dense elimination, grevlex, normal strategy, degree bound `n + d + 6`. No finished row had a pair above the bound | `f4-fp`, then `f4-fp-tower` |
| budget | 300 s per system (900–3600 s for the refutation-only cells at the largest `N`), and a 13 GB address-space cap | 600 s, 8 GiB |
| controls | `null`: `S₃`'s monomial support, random coefficients, and a planted zero on planted targets. `naive`: `x^{2^t} = g^{2^t}` per unknown, Kummer only. Generator-count `ladder` at `N = 8, 12` | also `random set` and the binary twin |
| not run | `macaulay-fp` and `δ_Z`, so Proposition 2 is untested here; the random-set control; the binary twin; `tree`; `m = 5` | — |

**Seeds.** Every cell draws its tower, curve and targets from `0x504B4D54`
combined with (kind, `m`, `t`, `g`). A cell is therefore the same system in every
run that contains it. Overlapping runs, made with three builds of the example,
repeated 215 measurements. Every repeat agreed on every deterministic field,
including those re-measured with the committed source.

**Planted targets stop early.** Once a planted system is solved, F4 keeps
processing pairs that reduce to zero, to certify the basis, and it climbs to
degree 7–11 doing so. At `N = 16` that tail, not the solve, reached about 10 GB.
The larger cells are therefore refutations, on random targets: above `t = 7`
at `m = 2`, and above `t = 3` or 4 at `m = 3`. `f4_fp` gained an opt-in stop
(`F4Options::stopping_below`). It halts once the staircase is at most a given
size, which bounds the number of solutions but leaves the basis uncertified.
It was used only for the re-runs named below.

**An accounting correction.** `f4_fp`'s `solving_degree` records the degree of
the *last* productive step. On tower systems the normal strategy climbs to
degree 5, learns low-degree elements, and descends again. The field therefore
under-reports: the Kummer `N = 12` refutation has its last productive step at
degree 3, after a step at degree 5. `F4Report` now also carries:
- `solving_degree_max`, the highest productive step degree, which is the
  framework's definition (`FRAMEWORK.md` §3);
- `max_cols_to_solution`, the widest matrix up to that step;
- `steps_to_solution`.

Every figure below uses them. The one smoke run made before the change is kept
with the data, and the analysis skips it.

**Correctness.** `verify.py` recounts every finished tower row by exhaustive
search over `V^m`, without F4.
- For `m = 3` it counts triples with `S₄ = Res_u(S₃(x₁,x₂,u), S₃(x₃,x_R,u)) = 0`,
  the chain system with `u` eliminated.
- For `m = 4` it counts quadruples with `S₅ = 0`, both free unknowns eliminated.
  `test_verify.py` checks that counter against the other order of elimination,
  and against planted decompositions.

All 176 rows agree: F4 refuted exactly the systems with no solution on the
grid, and found at least one solution in every other. No planted solution was
lost, across 310 distinct systems.

### 10.2 Results

**`m = 2`, `full`, `raw`.** At every `N`, every target gave the same `D`, so each
`D` below is exact, not a median. The other columns are defined as follows:
- *width* is the widest F4 matrix up to the solution;
- *steps* is the number of F4 steps to the solution (median);
- *time* is the median wall-clock F4 time on a shared 4-core container, a
  practicality note only.

| `N` | Kummer `D` | width | steps | time | isogeny `D` | width | steps | null `D` | naive `D` |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 2 | 4 | 9 | 2.5 | < 1 ms | 4 | 9 | 2.5 | 4 | 4 |
| 4 | 4 | 37 | 4.5 | < 1 ms | 5 | 48 | 4.5 | 4 | 6 |
| 6 | 4 | 137 | 6.5 | 1.5 ms | 5 | 146 | 8 | 4 | 10 |
| 8 | 4 | 282 | 8.5 | 7.9 ms | 5 | 327 | 10 | 4 | 18 |
| 10 | 4 | 615 | 11 | 91 ms | 5 | 632 | 13 | 4 | 34 |
| 12 | 5 | 4,427 | 14 | 1.1 s | 5 | 3,172 | 17 | 5 | 66 |
| 14 | 5 | 8,348 | 16 | 7.0 s | 5 | 7,225 | 20 | 5 | 130 |
| 16 | 5 | 9,427 | 18 | 21 s | 5 | 14,138 | 23 | 5 | — |
| 18 | 5 | 15,172 | 21 | 73 s | 5‡ | 17,754 | 25 | 5 | — |
| 20† | out of memory | | | | | | | | |

† At `p = 3221225473`. F4 ran 23 steps, none above degree 5, then exhausted
the 13 GB cap before refuting. No `D` is recorded.
‡ The second of the two systems is decomposable. Its first run spent 19
minutes certifying the basis and exhausted the memory cap. Re-run with the
staircase stop, F4 had pinned it at degree 5 to two points, the two orderings
of one decomposition.

- **`dickson` equals `kummer`** in `D`, width and steps at every `N` from 2 to 18.
  Their towers have the same leading terms and supports, and so do their `S₃`
  systems; only the coefficients differ. The two runs are therefore the same F4
  computation on different numbers, and `dickson` is not an independent
  replication. `isogeny` is: its step maps change from level to level, and it is
  not degenerate at infinity (§3.4).
- **`null` matches the tower** in `D`, width and steps. This holds for Kummer at
  every `N` it ran, and for Dickson and isogeny up to `N = 10`. One null system at
  `N = 16` timed out at 300 s, still certifying its basis at degree 7. Re-run with
  the staircase stop, F4 had pinned it down to at most one point at degree 5.
- **`naive` gives `D = 2^t + 2` exactly.**
- **At `p = 3221225473`,** Kummer gives the same `D`, widths and step counts at
  every `N` from 2 to 16. At that prime `|V|² ≪ p`, so the flat degree does not
  come from `|V|²` approaching `p` at the smaller prime.

**`m = 3`, `chain`, `raw`.**

| `N` | Kummer `D` | width | Dickson `D` | isogeny `D` | width |
|--:|--:|--:|--:|--:|--:|
| 6 | 5 | 368 | 5 | 6 | 399 |
| 9 | 6 | 4,407 | 6 | 6 | 3,543 |
| 12 | 6 | 12,482 | 6 | 6 | 13,249 |
| 15 | out of memory§ | | — | — | |

§ Two attempts. The traced one ran 19 steps, none above degree 6, with
matrices up to 27,379 columns, then exhausted the 13 GB cap before an answer.

**`m = 4`, `chain`, `raw`.** Random and planted targets agree. The planted runs
used the staircase stop at 64 points. It fired at 24–43, since every ordering
of a decomposition (at least `4! = 24`) solves the chain.

| `N` | Kummer `D` | width | isogeny `D` | width |
|--:|--:|--:|--:|--:|
| 8 | 6 | 1,916 | 6 | 1,897 |
| 12 | 7 | 21,768 | 7 | 21,841 |

At `N = 12` the solving degree rises by one per summand: 5, 6 and 7 for
`m = 2, 3, 4`. At `m = 4` there are only two sizes, and one refutation at
`N = 12` takes about two minutes.

At `N = 12`, two of the four Kummer systems are planted ones that hit the 300 s
budget in their certification tails. The two random ones finished at `D = 6`,
as did every Dickson and isogeny system.

**Generator-count ladder (instrument power).** On planted targets, `g − 1`
random quadrics through the planted point are added to the system:

| kind, `N` | `g = 1` | 2 | 4 | 8 | 12 |
|:--|--:|--:|--:|--:|--:|
| Kummer, 12 | 5 | 4 | 4 | 4 | 4 |
| isogeny, 12 | 5 | 5 | 5 | 5 | 4 |

At `N = 8`, Kummer and Dickson stay at 4 for `g = 1, 4`, and isogeny falls from 5
to 4. `D` cannot fall below 4, the degree of the raw `S₃`. Wherever it started
above 4, it fell as `g` grew, so the instrument registers a fall. It fell at
`g = 2` for Kummer at `N = 12`, at `g = 4` for isogeny at `N = 8`, and only at
`g = 12` for isogeny at `N = 12`.

### 10.3 The rule of §5.2, as written

| cell | β̂ (least squares) | pre-registered 95% interval | leave-one-`N`-out slopes | upper-half slope | §5.2 reading |
|:--|:--|:--|:--|:--|:--|
| Kummer, `m = 2` | 0.084 (`N` = 2…18) | [0.084, 0.084], degenerate | 0.079–0.095 | 0.100 (`N` = 10…18) | inconclusive |
| Dickson, `m = 2` | 0.085 (2…18) | degenerate | 0.079–0.095 | 0.100 (10…18) | inconclusive |
| isogeny, `m = 2` | 0.042 (2…18) | degenerate | 0.000–0.044 | 0.000 (10…18) | inconclusive |
| Kummer and Dickson, `m = 3` | 0.190 (6…12) | degenerate | 0.000–0.333 | — (three values) | inconclusive |
| isogeny, `m = 3` | 0.000 (6…12) | degenerate | 0.000–0.000 | — | inconclusive (no null at `m = 3`) |
| Kummer and isogeny, `m = 4` | 0.250 (8, 12) | degenerate | — (two values) | — | inconclusive (no null at `m = 4`) |
| naive (Kummer, `m = 2`) | 5.743 (2…12) | degenerate | 3.6–7.2 | 12.0 | H0 |

- **The pre-registered interval is degenerate.** Every target at a given `N`
  gave the same `D`. Resampling targets within `N` therefore returns β̂ every
  time, and the "95% interval" is a single point. The rule cannot be read at face
  value, so the jackknife and upper-half slopes are reported beside it.
- **H0 needs the lower bound above 0.25.** No tower cell meets it. The naive
  control does, as it should.
- **H1 needs two things:** the upper bound below 0.10, and `D_F4` below the
  shape-matched null at the largest `N`. The first holds only on the degenerate
  interval. The second fails in every cell where both ran: the null gives the
  same `D`, the same widths and the same step counts.
- **The reading is inconclusive for every kind.** §5.2 prescribes what to do
  next: extend `N` before deciding.
- **The expectation itself fails.** It was H0 for all three kinds, with β near 1
  for Kummer and Dickson. From `D = 4` at `N = 10`, β = 1 would give `D ≈ 12` at
  `N = 18`; the measured `D` is 5. Proposition 2 is a theorem about XL, and §3.3's
  estimate of about `N + d` is a statement about XL. Neither is tested here, and
  neither says anything about F4 with falls. What is refuted over the measured
  range is the prediction built on them (§0, §5.2): that F4 grows the same way.
  §3.5 had named the way out, falls.

### 10.4 What the controls say

- **naive against tower.** On the same `V` and the same targets, the solving
  degree drops from `2^t + 2` to 4 or 5: 66 against 5 at `N = 12`. That is the
  effect §3.1 hoped the tower would have, and it is much larger than §3.3
  predicted.
- **null equals tower.** The low degree is a property of the tower algebra and
  of `S₃`'s monomial support, not of the summation polynomial. That is exactly
  what the null was built to detect. For the oracle it cuts both ways: it needs
  nothing curve-specific, and it exploits nothing curve-specific either.
- **isogeny is as flat as Kummer.** Its `D` is 5 from `N = 4` on; Kummer's is 4,
  then 5 from `N = 12`. The effect therefore does not come from the Kummer group
  algebra, nor from its carries `y_j·y_j = y_{j+1}` (§2.4, §3.5). What the three
  families share is that every tower equation is a quadric with leading term
  `y_j²`.
- **ladder.** The instrument registers a fall once enough generators are
  added. The power check was designed to validate a negative. The pilot's
  result is not a negative, so the check carries less weight here.

### 10.5 What F4 does: an observation, not a proof

Two diagnostics were run on Kummer `m = 2` random targets.
- **The step trace** (`F4_DEBUG=1`). The first twelve steps, at degrees 3 and 4,
  have the same pairs, rows and pivots for `t = 6, 7, 8`. For `t = 5`, the first
  eight match. After those steps, each further level of the tower adds two or
  three steps, at degrees 4 and 5. The steps to the solution are 11, 14, 16, 18
  and 21 for `t = 5 … 9`. The steps at degree 5 number 1, 2 and 4 for
  `t = 6, 7, 8`.
- **The degree-capped closure** (`--cap 4 --dump 3`: F4 run to completion with
  every pair above degree 4 dropped). The result is the same for `t = 6, 7, 8`:
  - two quadrics, on levels 0–1 and 0–2 of both blocks;
  - two cubics on levels 0–2;
  - 88 cubics on levels 0–5, and none on any higher level, whatever `t` is.

  At degree 4, then, F4 knows a fixed window of the tower's lowest levels. What
  it needs from the levels above that window, it learns in the steps that
  follow, the first of which is at degree 5.

These observations suggest a mechanism.

> **Conjecture 3.** For `m = 2`, in all three families, F4 with the normal
> strategy solves the `raw` tower system at degree at most 5, for every `t`.
>
> *Suggested mechanism (level-window propagation).* Every level of a tower
> has the same shape. So the degree-5 closure of what F4 knows on levels
> `j … j + k` contains the corresponding equations on levels
> `j + 1 … j + k + 1`, and F4 climbs the tower one window at a time.

The conjecture fits the data, and nothing here proves it. For the isogeny
tower, whose maps change from level to level, only the shape repeats. A proof
would bound the last fall degree (Huang–Kosters–Yeo) of `m = 2` tower systems.
§10.9 ranks that as the third next step.

### 10.6 What it costs (a stage diagnostic, `AGENTS.md` §2)

- **Width.** While `D` stays fixed, the width is at most the number of
  monomials of degree at most `D`, `C(n + D, D)`, so it grows at most
  polynomially in `N`. At `N = 18` it is 15,172 of the possible 33,649. The
  slopes of `log₂(width)` per unit `N` fall from 0.70 (Kummer) and 0.73
  (isogeny) over `N = 2…18` to 0.52 and 0.59 over the upper half, `N = 10…18`,
  as a polynomial of fixed degree predicts. A slope over so short a range is not
  an exponent.
- **Time.** F4 time grows much faster than the width, because the elimination is
  dense: 1.1 s, 7.0 s, 21 s and 73 s at `N = 12, 14, 16, 18`.
- **Against the reference.** On the same system, the cell reference
  `tower-exhaustive` costs `2^t` square roots, which is 512 at `N = 18`. F4 has no
  operation counter yet (a Stage 0 item, §4.4), so the two are not in a common
  unit. For scale, F4's matrices at `N = 18` reach 15,172 columns and 25,791
  rows. This is no speed claim; nothing here is priced end to end.
- **§3.7 is unchanged.** At `m ≤ 3` no oracle, however fast, beats rho
  asymptotically. The measured cells test the mechanism, not the verdict.
- **No scoreboard row.** The pilot prices nothing, so it has no `S` to draw
  (`AGENTS.md` §7).

### 10.7 What the pilot shows, and what it does not

At the measured sizes (`N ≤ 18` at `m = 2`, `N ≤ 12` at `m = 3, 4`, `raw`, dense
`f4_fp`, two primes), the pilot shows four things:
1. The tower presentation removes the degree of the Yokoyama regime:
   `2^t + 2` becomes 4–7.
2. F4's solving degree on tower systems does not follow `N + d`. At `m = 2` it
   is 4–5 through `N = 18`, in all three families. At `m = 3` it is 5–6 through
   `N = 12`.
3. The effect depends neither on the summation polynomial (the null) nor on the
   Kummer group algebra (isogeny).
4. The degree grows with the number of summands. At `N = 12` it is 5, 6 and 7
   for `m = 2, 3, 4`.

It does not show the following:
1. **That `D` is bounded.** For Kummer and Dickson at `m = 2`, `D` moved once,
   from 4 to 5 at `N = 12`. For isogeny it has not moved since `N = 4`, but every
   `N` measured is small. One increase every 10 in `N` (β = 0.1) fits the data as
   well as a constant, and §3.6 says what β = 0.1 would mean: width
   `p^{H(0.1)} ≈ p^{0.47}` per target, a negative for the oracle. Only `N` of
   about 30–40 separates the two. In 13 GB the dense `f4_fp` reaches `N = 18`
   at `m = 2` and `N = 12` at `m = 3`. At `N = 20` and `N = 15` it ran out of
   memory before an answer, having stepped no higher than degree 5 and 6
   respectively.
2. **How `D` grows with `N` at `m ≥ 4`.** These are the only cells §3.7 allows
   to matter. `m = 4` has two sizes, with `D = 6` and 7, which is no slope at
   all. `m = 5` and `tree` were not run.
3. Any speed, any `S`, or anything about a deployed curve.

### 10.8 A proposed amendment to §5.2 (post hoc, for Stage A)

This amendment was written after seeing the pilot's data, and it is not applied
to it.
- **A1(a): a plateau criterion replaces the bootstrap.** The solving degree is
  an integer and is constant within each `N`, so a within-`N` resampling
  interval collapses. Report instead:
  - the per-`N` values;
  - the leave-one-`N`-out slopes;
  - the length `L` of the final plateau, the `N` range since `D` last rose.

  A slope below 0.10 needs `L > 10`. H0 needs an increase at least every 4 in
  `N` over the upper half.
- **A1(b): H1 splits in two.**
  - **H1a (viability):** `D` is bounded, judged by A1(a). This gates Stage B.
  - **H1b (curve-specificity):** `D` is below the null. It is reported and gates
    nothing.

  §5.2 conflated the two. A low degree that the null shares means the mechanism
  is generic in the polynomial, and that does not weaken the oracle.

Under A1, Kummer and Dickson still read inconclusive: at `m = 2` their final
plateau covers `N = 12…18`, so `L = 6`. Isogeny's plateau covers
`N = 4…18`, `L = 14`, which would pass A1(a). A1 was written after these data,
though, so that is a prediction for Stage A to test, not a result.

### 10.9 Next steps, ranked

1. **Extend `N` to 30–40 at `m = 2` and `m = 3`.** Nothing else separates a
   bounded degree from one increase every 10 in `N` (A1(a)). This needs the
   sparse, tower-aware engine of §4.4 (`f4-fp-tower`), or an external engine that
   has passed the known-answer battery. §10.5 suggests a design: apply the tower
   as a rewriting rule, and keep in the matrix only the window of levels F4 is
   working on.
2. **Run `m = 4` at four or more sizes, and `m = 5`, in `chain` and `tree`.**
   These are the decisive cells. The pilot's two `m = 4` sizes already take two
   minutes per refutation at `N = 12`, so they too need the sparse engine. The
   growth of `D` with `m` (5, 6, 7 at `N = 12`) decides how the per-target cost
   scales with `m` in §3.7's balance.
3. **Prove or refute Conjecture 3.** It is concrete enough to attack directly:
   show that the degree-5 closure of a window of levels contains the window one
   level up.
4. **Replicate independently** with another Gröbner engine (msolve after the
   battery, or Magma), and compare against whatever degrees PKM report
   (§10.10).
5. **Then build Stage 0,** as planned. Nothing in the pilot changes §4.

### 10.10 Relation to the literature, and what is at stake

- **PKM's paper could not be retrieved in this container.** The authors' host
  failed TLS, and the publisher's copy is paywalled. Whether PKM saw the same low
  degrees is therefore unchecked. Amadori–Pintore–Sala (ePrint 2017/609, §3.3)
  quote PKM's per-system times at `m = 3`, from 0.02 s at 11 bits to 5163 s at
  22 bits, and report no degrees. The flat degree here may reproduce PKM's own
  observation rather than add a new one; this note claims neither.
- **The stake.** If `D` is bounded for some `m ≥ 4`, §3.7's floor `r^{2/(m+1)}`
  becomes the method's exponent, up to factors polynomial in `log r`; at `m = 4`
  that is `r^{2/5}`. This is the exponent-moving outcome §0 describes. It rests on
  two things the pilot did not measure (§10.7), and is recorded as the stake, not
  as a finding.
- **Which primes the Kummer and Dickson families reach.** A Kummer tower of
  length `t` needs `2^t | p − 1`. A Dickson torus tower needs `2^t | p + 1`.
  Isogeny towers need neither. The values below, given for scale, are properties
  of the primes, not measurements:

  | prime of | `v₂(p − 1)` | `v₂(p + 1)` |
  |:--|--:|--:|
  | P-224 | 96 | 1 |
  | P-256 | 1 | 96 |
  | P-521 | 1 | 521 |
  | Ed448 | 1 | 224 |
  | P-192, SM2 | 1 | 64 |
  | P-384 | 1 | 32 |
  | secp256k1 | 1 | 4 |
  | brainpoolP256r1, BN254 | 1 | 3 |
  | Curve25519 | 2 | 1 |
  | BLS12-381 | 1 | 2 |

  Nothing in §10 runs at these sizes. The Dickson torus case is unmeasured,
  because the pilot's prime is split. A long tower is a precondition for these
  two families, not a weakness.

---

## 11. Round 2, 2026-09-25: a sparse tower engine, and `N` past the pilot

§§11.1–11.6 were written and committed before any round-2 cell ran; §11.3
discloses the one measurement at a new size made before them. The results are
added below them as §11.7 onward, and §§11.1–11.6 stay as committed. The data,
commands and build hashes are in `research/pkm_tower_round2_20260925/`.

### 11.1 The engine: `f4_fp_tower`

§10.9 ranked a sparse, tower-aware F4 first. The dense `f4_fp` had run out of
13 GB at `N = 20` (`m = 2`) and at `N = 15` (`m = 3`).
`src/cryptanalysis/f4_fp_tower.rs` is that engine: F4 on the quotient
`R = F_p[y, u]/(towers)` of §2.4.
- **The tower is a rewriting rule.** A monomial is square-free in the tower
  variables. It is stored as a 64-bit mask plus up to four exponents for the
  chain unknowns `u`, packed into one `u128` whose integer order is grevlex. A
  product is rewritten by `y_j² = a·y_j·y_{j+1} + b·y_{j+1} + c·y_j + d` as it
  is formed, so the tower equations never enter a matrix.
- **Pairs.** Critical pairs go through Gebauer–Möller (`UPDATE`). A *tower
  pair* `(g, y_j)`, one for each tower variable of `LM(g)`, stands for the
  S-polynomial of `g` against `y_j² − …`. It enters the matrix as the single row
  `NF(y_j·g)`. The selection is the normal strategy, as in `f4_fp`.
- **Elimination, in the manner of Faugère–Lachartre.**
  - Symbolic preprocessing gives reducer rows with distinct leads, and at every
    lcm column the shortest S-row joins them.
  - The columns left without a reducer are then exactly the monomials without a
    divisor.
  - The engine reduces the reducer block to dense rows over those columns
    (`B' = A⁻¹B`, right to left), then reduces every S-row by it in one pass
    (`D − C·B'`), then echelonizes the residues in parallel chunks. Every pivot
    it returns is a new leading monomial.
- **Arithmetic.** Below `2³¹` the dense kernel reduces lazily: it subtracts a
  multiple of `p` only when an accumulator crosses `2⁶³`. That loop vectorizes,
  and on x86-64 it is dispatched at run time to AVX2 or AVX-512. Above `2³¹`
  every step is reduced (Barrett).
- **Rows are formed on demand.** A row is kept as the product that makes it (a
  multiplier and an element) and formed again, a block at a time, when it is
  eliminated. The memory is then the dense block `B'`, not the matrix.
- **What it measures.** The same statistics as `f4_fp`: `solving_degree_max`,
  the width to solution and the steps to solution. It also counts `F_p`
  multiply-adds, the engine's own unit. Its input is the *reduced*
  presentation of §2.2: a summation polynomial enters as its normal form in `R`.

The engine has eight unit tests (`cargo test --release --lib f4_fp_tower`):
- the grevlex keys;
- products against a generic reduction;
- refutation and solution counts against brute force on random Kummer blocks;
- a free unknown eliminated like a chain link;
- the staircase stop keeping every solution;
- the degree bound;
- the elimination against a plain sequential echelon on random matrices, at
  three primes on both sides of `2³¹`;
- every kernel against exact arithmetic.

For scale only (wall clock is a practicality note, `AGENTS.md` §2): the
pilot's Kummer `m = 2` random system at `N = 18` takes 6.2 s and 139 MB of
memory, where `f4_fp` took 73 s (§10.2).

### 11.2 Cross-check against `f4_fp` on the pilot's systems

`AGENTS.md` §6 does not accept a new oracle that has not been cross-checked
against the one it replaces, on every input of a full run. The pilot's cells
were re-run with `--engine tower`, with the same seeds and the pilot's flags.
These are runs `XV1`–`XV5`. Some cells were widened to every kind and control,
but no size exceeds the pilot's largest for its `m`.

The five runs cover:
- `m = 2`: all three kinds, tower and null, `N = 2–18`, with the planted
  targets up to `N = 14`, plus Kummer at `p = 3221225473`;
- `m = 3`: all three kinds, `N = 6–12`;
- `m = 4`: Kummer and isogeny, `N = 8` and 12, on random targets and on
  planted ones with the staircase stop at 64;
- the generator-count ladder at `N = 8` and 12.

Together they took about nine minutes. Every system finished, and none lost a
planted solution.

| part | systems | finished by both engines | same `D` | same verdict |
|:--|--:|--:|--:|--:|
| `m = 2`, tower (both primes) | 126 | 126 | 110, and 16 at `N = 2` | 126 |
| `m = 2`, null | 96 | 72 | 60, and 12 at `N = 2` | 72 |
| ladder | 40 | 28 | 28 | 28 |
| `m = 3` | 30 | 28 | 28 | 28 |
| `m = 4` | 16 | 16 | 16 | 16 |

- **The verdict agrees on every system both engines finished.** Whether a
  system was refuted never differs between them. `analyze.py` counts 269
  systems. It skips one more, the Kummer null at `N = 16` that `f4_fp` finished
  only in its staircase-stop re-run (R2 of the pilot), and that system agrees as
  well.
- **The solving degree agrees on every one of those systems from `N = 4` up**,
  in every kind, control and `m`, and at both primes.
- **It differs only at `N = 2`, on 28 systems**, where the tower engine gives
  `D = 3` and `f4_fp` gives 4. With one level, `S̄` has degree 2 in `R`, whereas
  the raw `S₃` has degree 4, which is `f4_fp`'s floor (§10.2). From `N = 4` on,
  the input degree no longer decides `D`. There, Kummer's and Dickson's `S̄`
  have degree 2 in `R` and isogeny's has degree 4, against 4 for the raw `S₃`,
  and the two engines still agree.
- **38 systems have no `f4_fp` counterpart to compare with**:
  - 24 are Dickson and isogeny nulls above `N = 10`, and 12 are ladder rungs at
    `N = 8`. The pilot never ran either.
  - Two are Kummer `m = 3` random targets at `N = 12`. The pilot's cell also
    planted at that size, so it drew different targets. Both give `D = 6`, as
    the pilot's did.
- **`verify.py` confirms all 172 tower rows** against exhaustive search.

The cross-check therefore reproduces every solving degree of §10.2 above `N = 2`
with the new engine. That includes the two systems `f4_fp` could finish only
with the staircase stop.

### 11.3 What was measured before this section (a disclosure)

Tuning the engine used the pilot's Kummer `m = 2` random target 0 at
`N = 14, 16, 18`, which are systems the pilot had measured. It also used, once,
the same cell at `N = 20` at `p₀ = 786433`. The pilot never finished a system at
that size: its `N = 20` attempt was at `p = 3221225473` and ran out of memory.
That system:
- **gave `D = 6`**. F4 ran 22 steps at degrees 3–5. The degree-5 pairs thinned
  out, to 200 and then 20 S-rows in the last two such steps, without a fall to
  the unit. Step 23, at degree 6, had 18,305 critical and 12,030 tower pairs
  and 33,446 columns. It gave 2,042 new elements, the lowest of degree 3, and
  step 24, at degree 4, found `1`;
- was re-run with a later build, the one that forms rows on demand, and its
  step trace agrees with the first line for line.

It is the only measurement at a new size before this section. It contradicts
Conjecture 3 (`D ≤ 5` at `m = 2` for every `t`), if the engine is right at
`N = 20`, which only this engine reaches.

It shaped the cells below. At `N = 20`, `|V|² = 2^20 = 1.33·p₀`. That is the
first size on the pilot's line where the grid is larger than the field; at
`N = 18` it is `0.33·p₀`. PKM's own regime is `|V|^m ≈ m!·r` (§2.3), so the rise
may come from `N` or from `|V|²/p`, and the round separates the two. It did not
shape the decision rule, which is A1 as written in §10.8, before it.

### 11.4 Cells

The new prime is `p₁ = 2013265921 = 15·2^27 + 1`. Its `v₂(p₁ − 1) = 27` allows
Kummer towers up to `t = 27`, and `p₁ < 2³¹` keeps the lazy kernel. The isogeny
family stays at `p₀`, because the pilot's construction of an isogeny tower
counts points naively, in `O(p)`.

| cell | kind | `m` | prime | `N` | targets per `N` | control | note |
|:--|:--|--:|:--|:--|:--|:--|:--|
| K1 | Kummer | 2 | `p₁` | 12, 14, …, 22 | 2 random | tower | `\|V\|²/p₁ ≤ 0.002`: every system refutes with high probability |
| K1n | Kummer | 2 | `p₁` | 12, 14, …, 20 | 2 random | null | the shape-matched null of §5.3 |
| K0 | Kummer | 2 | `p₀` | 20 | 2 random | tower | target 0 is §11.3's system; staircase stop at 8 |
| I0 | isogeny | 2 | `p₀` | 20 | 2 random | tower | staircase stop at 8 |
| M4 | Kummer | 4 | `p₁` | 8, 12, 16 | 2 random | tower | `chain` |
| M3 | Kummer | 3 | `p₁` | 9, 12, 15, 18 | 2 random | tower | `chain` |

- **Budgets.** 7,200 s per system, a 14 GB address-space cap and
  `--max-nnz 2000000000`. A system that times out, stops for size or runs out of
  memory reports no `D`, at most a lower bound, and never "does not decompose".
- **Fixed as in the pilot.** The degree bound is `n + d + 6`, and the seed is
  `0x504B4D54` with the pilot's cell rule, so a cell is the same system in every
  run.
- **Build.** The engine and the example are as committed with this section.
- **Order.** The cells run in the table's order, one system at a time.
- **Check.** Every finished tower row goes through `verify.py`.
- **Confirmation runs.** After the cells, every `m = 2` cell with a finished
  system at `D ≥ 6` is re-run once, whole, with the degree bound at 5
  (`--cap 5`) and otherwise the same flags. §11.6 says what those re-runs must
  show.

### 11.5 Predictions, before the runs

1. **K1 at `N = 12`–18:** `D = 5` on every system, as at `p₀` and at
   `3221225473` (§10.2).
2. **K1 at `N = 20` and 22:** `D = 6`. That is, the rise of §11.3 is a function
   of `N`, not of `|V|²/p`. The reading behind it is that, once `t = 10`, the
   level window of §10.5 no longer closes at degree 5. The confidence is low:
   the prediction rests on one system.
3. **K0:** `D = 6` on every system that refutes, and target 0 reproduces §11.3.
4. **K1n:** the null gives the tower's `D` at every `N` (§10.4).
5. **I0:** `D ≥ 5`, and nothing more specific.
6. **M4 and M3:** `D = 6` and 7 at `N = 8` and 12 for `m = 4`, and `D = 6` at
   `N = 9` and 12 for `m = 3` (§10.2). At the new sizes (`N = 16` for `m = 4`;
   15 and 18 for `m = 3`) there is no numerical prediction beyond `D` not
   falling.

### 11.6 Decision rule, fixed now

- **A1 of §10.8, adopted as written, per cell** (kind, `m`, prime). It reports
  the per-`N` values, the leave-one-`N`-out slopes and the length `L` of the
  final plateau.
  - H1a (bounded `D`) needs `L > 10`, and it gates Stage B.
  - H0 needs an increase at least every 4 in `N` over the upper half of the
    cell's range.
  - Anything else is inconclusive.
  - H1b (`D` below the null) is reported and gates nothing.
  - §5.2's bootstrap is not computed, because it is degenerate (§10.3).
- **What this round can reach, stated in advance.**
  - Kummer's plateau at `m = 2` starts at `N = 12` (§10.2), so H1a needs `D`
    constant through `N ≥ 24`, past this round's largest size. For Kummer at
    `m = 2`, this round can therefore return H0 or inconclusive, not H1a.
  - The same holds at `m = 3`, whose plateau starts at `N = 9` so that H1a needs
    `N ≥ 21`, and at `m = 4`.
- **Conjecture 3** is refuted by any finished `m = 2` tower row with `D ≥ 6`
  that meets two conditions. `verify.py` must confirm its verdict. Its
  confirmation run (§11.4) must confirm the degree-6 step: with the degree bound
  at 5, the same system must end with pairs above the bound, without refuting
  and without reaching a staircase stop.
- **`N` or `|V|²/p`.**
  - "`N`" holds if both K1 systems at `N = 20` have `D = 6`.
  - "`|V|²/p`" holds if all four K1 systems at `N = 20` and 22 have `D = 5`
    while K0's refuted systems have `D = 6`.
  - Anything else is reported as mixed.
- **Instrument first.** A verdict that `verify.py` contradicts voids the
  round's reading until it is explained, as in §10.1.
- **Inadmissible**, beyond §5.4:
  - changing the engine, the cells, the budgets or the seeds after the first
    round-2 run;
  - dropping a target;
  - re-running a system in the hope of a different answer.

### 11.7 Results

**Build.**
- After the pre-registration commit (`a2a851ed`), `origin/main` was merged into
  the branch (`d8b1df8f`).
- The merge touched no file the example builds from except the module list in
  `src/cryptanalysis/mod.rs`, so the engine and the example are the sources of
  §11.1. Their hashes are in the round's README.
- The merged build reproduces XV2 and XV5 exactly: 60 systems, every
  deterministic field.
- Every round-2 run used that build.

**K1: Kummer, `m = 2`, `p₁`.** Every system refuted.

| `N` | `D`, targets 0 and 1 | width | steps | F4 time, each | nonzeros, widest step |
|--:|:--|--:|--:|--:|--:|
| 12 | 5, 5 | 1,231 | 13 | 0.05 s | 0.3 M |
| 14 | 5, 5 | 2,304 | 15 | 0.2 s | 1.2 M |
| 16 | 5, 5 | 4,789 | 17 | 1.1 s | 4.6 M |
| 18 | 5, 5 | 8,644 | 20 | 6 s | 18 M |
| 20 | **6, 6** | 33,446 | 24 | 57 s | 121 M |
| 22 | **6, 6** | 51,909 | 28 | 10 min | 473 M |

- **The shape of the system decides everything measured here.**
  - At each `N` the two targets give the same `D`, width and step count.
  - So do the systems at `p₀` and at `3221225473`, wherever they overlap.
  - Their counts of nonzeros and multiply-adds agree exactly at the two large
    primes. At `p₀` they differ by a few parts per million, from cancellations
    modulo a 20-bit prime.
- **Memory.** The largest system, at `N = 22`, peaked at 2.3 GB. Its matrices
  hold 473 million nonzeros, which are never in memory together.

**K1n: the null at `p₁`.** It has the same `D`, width and step count as K1 at
every `N` from 12 to 20, including `D = 6` at `N = 20`.

**K0 and I0: `m = 2` at `p₀`, `N = 20`.** Every system refuted, so the
staircase stop never fired.
- **K0** gives `D = 6` on both targets, at 33,446 columns in 24 steps. Target 0
  reproduces §11.3 field for field.
- **I0, the isogeny tower, gives `D = 6` on both targets**, at 35,158 columns in
  28 steps. The isogeny `D` had been 5 at every `N` from 4 to 18. It rises at
  the same level as Kummer's, `t = 10`.

**M4: Kummer, `m = 4`, `p₁`.**
- `N = 8` and 12 give `D = 6` and 7 on both targets, as in the pilot at `p₀`:
  1,049 and 9,452 columns, in 21 and 29 steps.
- A trace of the `N = 12` system (D1 in the README) shows how `m = 4` refutes.
  One early step reaches degree 7 (step 16). The refutation then comes through a
  long tail of degree-6 steps, in which the basis grows to 1,886 elements and
  the pending pairs to 25,873.
- **`N = 16` ran out of memory.** Its first system hit the 14 GB cap after 29
  minutes ("memory allocation of 131072 bytes failed"). The process ended
  there, so the second target never ran, and there is no row and no `D`.
- **A diagnostic re-run of that system (D2) gives a lower bound.** It is labelled
  as such, runs on the same build with a 25-minute budget so that the engine
  stops itself before memory runs out, and records its trace.
  - In 48 steps it reached degree 7 (at step 12) and no higher, over up to
    76,837 columns. It had not refuted when the budget ran out, so
    `D ≥ 7` at `N = 16`, against 7 at `N = 12`.
  - Every step from step 12 on is at degree 7. The last steps grow fast: step
    48 had 4,117 critical pairs and 229 million nonzeros, with a basis of 6,036
    elements and 123,032 pending pairs. That tail is what filled the memory.

**M3: Kummer, `m = 3`, `p₁`.**
- `N = 9` and 12 give `D = 6` on both targets, as in the pilot.
- **`N = 15` gives `D = 7` on both targets**, at 31,397 columns in 24 steps,
  taking about 6 minutes each.
- **`N = 18` stopped for size on both systems.** Each ran about 44 minutes and
  23 steps, reaching degree 7 at 104,600 columns. The `--max-nnz` cap of
  `2·10⁹` then stopped step 24. `D ≥ 7` is a lower bound, not a value.

**Confirmation runs.**
- Every `m = 2` cell with a system at `D ≥ 6` was re-run whole, with the degree
  bound at 5 (`confirm_round2.sh`: five cells, ten systems).
- Every system ended with pairs above the bound, without refuting and without a
  staircase stop. That is 30,335 pairs for each Kummer system and the null, and
  14,828 for isogeny. Degree 5 does not suffice under this engine's strategy.
- **At degree 5, the level window of §10.5 has stopped growing.** The capped
  Kummer runs at `N = 20` and `N = 22` end with the same 22 steps, the same
  2,922 basis elements and the same 30,335 pending pairs. Only their column
  counts differ (15,740 and 16,572). What F4 learns at degree 5 no longer
  depends on `t`.

**Exhaustive check.**
- `verify.py` confirms the verdict of all 198 finished tower rows, 26 of them
  new in round 2. It skips the 8 capped tower rows, which claim none.
- No planted solution was lost.
- The 218 repeated measurements agree on every deterministic field. They include
  K0 against §11.3's profiling rows, and D1 against M4.

**The exit codes in `progress.txt` are all 0 and carry no information.** The
committed scripts wrote them with `$(date …)` in the same line, which resets
`$?` before it is read. The logs record how each run ended, and
`confirm_round2.sh` reads the status first.

### 11.8 The rule of §11.6, applied

**The predictions of §11.5.**

| | prediction | outcome |
|:--|:--|:--|
| 1 | K1, `N = 12–18`: `D = 5` | held on all 8 systems |
| 2 | K1, `N = 20, 22`: `D = 6` | held on all 4 |
| 3 | K0: `D = 6`; target 0 reproduces §11.3 | held |
| 4 | K1n: the null's `D` is the tower's | held at every `N` |
| 5 | I0: `D ≥ 5` | held: `D = 6` |
| 6 | M4: 6 and 7 at `N = 8, 12`; M3: 6 at `N = 9, 12`; `D` not falling at new sizes | held where measured. M3 gives 7 at `N = 15`, and both `N = 18` systems stopped for size at `D ≥ 7`. M4 gave only `D ≥ 7` at `N = 16` (D2) |

**`N`, not `|V|²/p`.** Both K1 systems at `N = 20` have `D = 6`, so "`N`"
holds. At `p₁` the grid is `2^20/p₁ ≈ 5·10⁻⁴` of the field, three orders of
magnitude below §11.3's system, and the run has the same shape step for step.

**Conjecture 3 is refuted**, in the sense §11.6 fixed.
- Eight finished `m = 2` tower systems have `D = 6`: K1 at `N = 20` and 22, K0
  and I0.
- `verify.py` confirms each of their verdicts.
- Their confirmation runs confirm each degree-6 step.
- The refutation concerns the reduced presentation that the engine measures.
  The raw presentation agrees with it on every pilot system from `N = 4` up
  (§11.2), but no raw system at `N = 20` has finished.

**A1, per cell.** The same `D(N)` reads differently depending on where a cell's
range starts, so each reading comes with its range.

| cell: kind, `m`, prime (runs) | `N` | `D` | final plateau `L` | upper-half slope | reading |
|:--|:--|:--|--:|:--|:--|
| Kummer, 2, `p₁` (K1) | 12–22 | 5 5 5 5 6 6 | 2 | 0.25 (`N = 18–22`) | H0, at the threshold |
| null, Kummer, 2, `p₁` (K1n) | 12–20 | 5 5 5 5 6 | 0 | 0.25 (16–20) | H0, at the threshold |
| Kummer, 2, `p₀` (XV1, K0) | 2–20 | 3 4 4 4 4 5 5 5 5 6 | 0 | 0.10 (12–20) | inconclusive |
| isogeny, 2, `p₀` (XV1, I0) | 2–20 | 3 5 5 5 5 5 5 5 5 6 | 0 | 0.10 (12–20) | inconclusive |
| Kummer, 3, `p₁` (M3) | 9–15; 18 stopped at `≥ 7` | 6 6 7 | 0 | — | inconclusive (three sizes) |
| Kummer, 4, `p₁` (M4) | 8–12; 16 only `≥ 7` (D2) | 6 7 | 0 | — | inconclusive (two sizes) |

- **Both H0 readings are at the threshold, and both come from where the range
  begins.**
  - At `p₁` the cells start at `N = 12`, so each upper half is three points with
    one rise, and its slope is exactly the threshold, 1/4.
  - The same `D(N)`, sampled from `N = 2` at `p₀`, has an upper-half slope of
    0.10 and reads inconclusive.
  - The rule's reading is recorded as it stands, and it closes nothing. §3.7
    already rules `m = 2` out: no oracle beats rho there.
- **H1a holds in no round-2 cell.** No final plateau is longer than 2 (K1's).
  The isogeny plateau of `N = 4–18`, which §10.8 noted would pass A1, ends at
  `N = 20`.
- **One line still reads H1a: the isogeny null at `p₀`.** The cross-check ran it
  only to `N = 18`, so it keeps that plateau. The null gates nothing, and at
  `N = 20` round 2 ran the isogeny tower, not its null.
- **`m = 4`, the regime that decides the oracle (§3.7), has two sizes**, and
  the rule cannot read it.

### 11.9 What round 2 shows, and what it does not

It shows five things.
1. **A sparse engine that matches the dense one.** `f4_fp_tower` gives
   `f4_fp`'s verdict on every pilot system, and its solving degree from `N = 4`
   up. It reaches `N = 22` at `m = 2` in 2.3 GB, where the dense `f4_fp` ran
   out of 13 GB at `N = 20`.
2. **At `m = 2` the solving degree is not bounded by 5.** It rises to 6 at
   `N = 20`, in both families, at both primes and in the null. That refutes
   Conjecture 3 (§11.8).
3. **The rise depends on `N` alone** among the things varied here. It does not
   depend on `|V|²/p`, the target, the curve or the summation polynomial: within
   a family, the whole run has the same shape.
4. **At `m = 3` the degree rises too**, to 7 at `N = 15`. At `N = 18` both
   systems stopped for size at degree 7, a lower bound.
5. **The rises fall at regular levels.** This pattern was noticed after the
   data, and is stated as a pattern, not a law.
   - At `m = 2`, the Kummer degree rises at `t = 6` and at `t = 10`, which is
     every four levels. The isogeny degree, 5 from `t = 2`, also rises at
     `t = 10`.
   - At `m = 3` the rises fall at `t = 3` and `t = 5`.
   - At `m = 4`, at `t = 3`.

It does not show three things.
1. **Whether `D` grows linearly in `N`.**
   - Two rises at `m = 2` fit one rise every 8 in `N`, `β = 1/8`. By §3.6 that
     is a width of `p^{H(1/8)} ≈ p^{0.54}` per target, a negative for the oracle
     if it holds at every `m`.
   - They fit slower growth equally well. Rises at geometrically spaced `N`
     (12, 20, then about 33) would make `D` logarithmic in `N`, which §3.6
     counts as viable.
   - Telling the two apart needs `N ≈ 28–33` at `m = 2`.
2. **Anything at `m = 4` beyond `N = 12`, the regime that decides the oracle,**
   except the lower bound `D ≥ 7` at `N = 16`.
3. **Any speed, `S` or scoreboard row.** Nothing is priced end to end, so there
   is no row to draw (`AGENTS.md` §7), as in §10.6.

### 11.10 Next steps, ranked

1. **Reach `m = 4` at `N ≥ 16`.** The traces of §11.7 show where the memory
   goes: into a long tail of steps at the top degree. That is degree 6 at
   `N = 12` and degree 7 at `N = 16`, where in 25 minutes the basis (every
   element of which is kept, retired ones included) passes 6,000 elements and
   the pending pairs 123,000. Two levers:
   - **Memory only, the same algorithm:** retire basis elements that no pending
     pair references, and store elements over a shared monomial table rather
     than as `u128` keys.
   - **A different algorithm:** signature-based criteria (F5 or GVW) to skip
     most of the pairs that reduce to zero. Its solving degree must be checked
     against this engine's first.
2. **Separate linear from logarithmic growth at `m = 2`.** Test the next rise
   the "every four levels" pattern predicts, at `t = 14` (`N = 28`). This needs
   the second lever.
3. **Prove or refute the pattern.** It is a `t`-dependent version of
   Conjecture 3: one more degree for every four more levels.
4. **Replicate `D = 6` at `N = 20`** with an independent Gröbner engine, as
   §10.9's fourth step proposed.

## 12. Round 3, 2026-09-25: a compact basis, and `m = 4` at `N = 16`

§§12.1–12.6 were written and committed before any round-3 cell ran; §12.3
discloses what was measured at a new size before them. The results are added
below them as §12.7 onward, and §§12.1–12.6 stay as committed. The data,
commands and build hashes are in `research/pkm_tower_round3_20260925/`.

### 12.1 Where the memory went, and the change

§11.10 ranked first reaching `m = 4` at `N ≥ 16`, and its first lever was memory
only: retire unused basis elements and store the rest compactly. Before choosing,
a profiling run (§12.3) measured the basis and the dense block on round 2's
`N = 16` system.
- **After step 48 the basis held 11,304 elements with 274 million terms.** A term
  was a `u128` monomial and a `u32` coefficient, 20 bytes, so the basis took
  about 5.5 GB. Step 48 alone added 4,926 elements with 159 million terms. A new
  element's tail is dense over the step's columns without a divisor.
- **Retiring unused elements would free almost nothing.** Only 6 elements were
  inactive with no pending pair using them. 169 million terms sat in elements
  that are inactive but still used by pending pairs. An element made inactive
  keeps its pending pairs, and its tower pairs are never pruned (§11.1).
- **Step 49 is larger again.** It is at degree 7, with 31,902 S-rows, 29,454
  pivot rows and 76,463 columns, 47,009 of them without a divisor. Its `B'` has
  1.1·10⁹ entries (4.4 GB), and the echelon of its residues can grow as large.
  With the basis at 5.5 GB that passes the 14 GB cap, where round 2's run died.

**The change: the basis keeps the echelon's rows.**
- A new element is kept as the pivot row the elimination leaves: dense over the
  step's columns without a divisor, from its lead on, zeros included, at 4 bytes
  an entry. The elements a step adds share that step's column list, and products
  skip the zeros.
- The elimination no longer converts its pivots into sparse rows and then into
  polynomials, copies at 8 and 20 bytes a term that were held alongside the
  dense rows.
- Step 48's new elements are 203 million entries in this form, 0.8 GB, where
  they took 3.2 GB. For the whole basis after step 48 the same ratio gives about
  1.4 GB. That figure is an estimate, not a measurement.

**Nothing else changes.**
- The pairs and their criteria, the selection, symbolic preprocessing, the
  elimination and its arithmetic, and the order in which elements enter the
  basis are all as in §11.1.
- Every multiply-add is the same one, so the unit of §11.1 does not move.
- The trace gains the kept elements and their entries, `B'`'s entries and the
  process's memory. The example prints each step as it ends, so that a run that
  dies still leaves its trace.
- A ninth unit test checks a stored element against its polynomial: its terms
  and its products, with zeros inside it and after it and with entries before
  its lead that belong to others.

By `AGENTS.md` §3 this is engineering: it moves memory, not the algorithm.

### 12.2 The identity check, before any cell

A memory change is admissible only if it changes nothing else. Before any round-3
cell:
- **Replay.** `replay_round2.sh` re-runs, with the new build and the same flags,
  every system that round 2 and its cross-check finished: XV1–XV5, K1, K1n, K0,
  I0, M4 at `N = 8` and 12, M3 at `N = 9`–15, the five confirmation cells and D1.
  It leaves out the sizes round 2 did not finish.
- **Every field.** `compare_builds.py` requires every field of every row to equal
  round 2's, except the wall clock (`ms`, `wall_s`). The multiply-add counts,
  nonzeros and residue sizes are among those fields. Every step of every trace
  that both builds printed must also equal round 2's, on every field round 2
  printed except its times.
- **Paired runs.** `pair_memory.py` then runs four systems that round 2 finished,
  each with the round-2 build and then the new one, and records each run's peak
  resident memory and wall time. The systems are target 0 of K0 (`N = 20`), D1
  (`m = 4`, `N = 12`), M3 at `N = 15` and K1 at `N = 22`. The two rows of each
  pair must agree too.
- **Any difference voids the change.** No cell runs until the difference is
  explained and fixed, and then the whole check runs again.

### 12.3 What was measured before this section (a disclosure)

The profile of §12.1 comes from one run of round 2's D2 system: Kummer, `m = 4`,
`p₁`, `N = 16`, target 0. It ran on the round-2 engine with temporary counters
(`profile/profiling.patch`, in no build since), from 07:46:45 to 08:12:55 UTC
today. Its log is `profile/profile-kummer-m4-p1-N16.log`.
- It repeated round 2's 48 steps, the steps D2's trace shows, in 292 s, and
  reported step 49's matrix before eliminating it. Those are the numbers of
  §12.1.
- It computed step 49's `B'` (227 s). It had been in step 49 for 21 minutes,
  about 17 of them on the S-rows, when it was stopped by hand (exit 143). That
  was before the step ended, so that the step's outcome would not be known before
  this section. Nothing it printed says whether step 49 finds new elements, or
  `1`.
- **It shaped the prediction of §12.5.** Step 49 has 4.9 times as many S-rows as
  step 48 (31,902 against 6,524). At `N = 12`, the step that found `1` had 3.1
  times as many as the step before it (7,440 against 2,374; D1).

### 12.4 Cells

| cell | kind | `m` | prime | `N` | targets | control | note |
|:--|:--|--:|:--|:--|:--|:--|:--|
| M4b | Kummer | 4 | `p₁` | 16 | 2 random, round 2's M4 targets | tower | the systems round 2 could not finish |

- **Budgets.** 14,400 s per system (round 2 allowed 7,200), the 14 GB
  address-space cap and `--max-nnz 2000000000`. A system that times out, stops
  for size or runs out of memory reports at most a lower bound. It never reports
  a value, and never "does not decompose".
- **Fixed as before.** The degree bound is `n + d + 6`, and the seed and cell rule
  are the pilot's, so M4b's systems are round 2's M4 systems at `N = 16`. The live
  trace (`--trace`) is on.
- **Build.** The engine and the example are as committed with this section, and
  run only once §12.2's check has passed.
- **Check.** Every finished row goes through `verify.py`, and `analyze.py` reads
  round 3 together with the pilot and round 2. M4b's target 0 is D2's system, so
  its first 48 steps must repeat D2's trace, which `compare_builds.py` checks.
- **Confirmation.** If a system finishes with `D = 7`, the cell is re-run once,
  whole, with the degree bound at 6 (`--cap 6`). Each such system must then end
  with pairs above the bound, without refuting and without a staircase stop.
  That confirms that degree 6 does not suffice, as §11.6 required at `m = 2`.

### 12.5 Predictions, before the runs

1. **Identity.** The replay and the paired runs reproduce round 2 on every
   compared field.
2. **Memory.** The new build's peak is below round 2's on every paired system.
   Both M4b systems get through step 49 inside the 14 GB cap.
3. **M4b: `D = 7` on both targets, and both refute at degree 7.** The confidence
   is low. The prediction rests on step 49's shape (§12.3): a jump in S-rows like
   the one that ended `N = 12`. The alternative is that degree 7 does not close
   and F4 moves on to degree 8. This machine may not reach that far, and the
   round would then give only a lower bound.
4. **Both targets give the same `D`, width and step count**, as in every cell of
   round 2.

### 12.6 Decision rule, fixed now

- **A1 of §10.8, per cell, as in §11.6.** The Kummer `m = 4` line at `p₁` then has
  `N = 8, 12, 16`.
  - If `D` is 6, 7, 7, the final plateau is `L = 4` and the upper half
    (`N = 12`–16) has no increase: inconclusive.
  - If `D` is 6, 7, 8, the upper half has one increase in 4, the threshold: H0,
    at the threshold.
  - H1a needs `L > 10`, which no round-3 outcome can give. The round can return
    H0 or inconclusive at `m = 4`, not H1a.
  - A lower bound is not a value. The line then keeps two sizes, and the rule
    cannot read it.
- **What each outcome would mean for §3.7**, stated now and read after the runs.
  - **`D = 7` at `N = 16`:** one level without a rise at `m = 4`. That fits slow
    growth and closes nothing.
  - **`D = 8`:** a rise at every level of `m = 4` so far (`t = 2, 3, 4`). If that
    continued, it would be one rise every 4 in `N`, `β = 1/4`, and a width of
    `p^{H(1/4)} ≈ p^{0.81}` per target by §3.6. The balance of §3.7 at `m = 4`
    cannot absorb that.
- **Instrument first.** A verdict that `verify.py` contradicts voids the round's
  reading until it is explained.
- **Inadmissible**, beyond §5.4 and §11.6:
  - changing the engine, the cell, the budgets or the seeds after the first
    round-3 run;
  - dropping a target;
  - re-running a system in the hope of a different answer.

  A run that the machine ends, not the engine, may be re-run once, and is
  reported as such.

---

## References

- C. Petit, M. Kosters, A. Messeng, *Algebraic approaches for the elliptic curve
  discrete logarithm problem over prime fields*, PKC 2016, LNCS 9615, 3–18.
- A. Amadori, F. Pintore, M. Sala, *On the discrete logarithm problem for
  prime-field elliptic curves*, Finite Fields Appl. 51 (2018); ePrint 2017/609,
  §2.2 (system (6)).
- K. Yokoyama, M. Yasuda, Y. Takahashi, J. Kogure, *Complexity bounds on Semaev's
  naive index calculus method for ECDLP*, J. Math. Cryptol. 14 (2020) 460–485.
- M. Kudo, Y. Yokota, Y. Takahashi, M. Yasuda, *Acceleration of index calculus
  for solving ECDLP over prime fields and its limitation*, CANS 2018.
- I. Semaev, *New algorithm for the discrete logarithm problem on elliptic
  curves*, ePrint 2015/310.
- K. Karabina, *Point decomposition problem in binary elliptic curves*, ICISC
  2015.
- M.-D. Huang, M. Kosters, S. L. Yeo, *Last fall degree, HFE, and Weil descent
  attacks on ECDLP*, CRYPTO 2015.
- S. D. Galbraith, S. W. Gebregiyorgis, *Summation polynomial algorithms for
  elliptic curves in characteristic two*, INDOCRYPT 2014.
- M. Bardet, J.-C. Faugère, B. Salvy, *On the complexity of the F5 Gröbner basis
  algorithm*, J. Symbolic Comput. 70 (2015).
- A. Caminata, E. Gorla, *Solving multivariate polynomial systems and an
  invariant from commutative algebra*, WAIFI 2020.
- J. Vélu, *Isogénies entre courbes elliptiques*, C. R. Acad. Sci. Paris 273
  (1971).
