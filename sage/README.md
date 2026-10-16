# SageMath scaffolding

Sage companions to the Rust crypto library. Each file is self-contained and
runnable with `sage <file>.sage`.

## `low_hw_factor_base.sage` — low Hamming-weight factor bases

Enumerates the special factor base `{P : Hamming weight of x(P) <= w}`.

- `ecc2k130` weighs x in the permuted type-II optimal normal basis, which is
  stable under Frobenius. Weight 2 is the fixed-campaign base: 14 signed
  orbits and 3,668 signed points. `--basis polynomial` uses the polynomial
  basis instead; that weight is not Frobenius-stable.
- `p192`, `p224`, and `p256` weigh the integer x in `[0, p)` and keep one
  point of each `{P, -P}` pair.

```
sage low_hw_factor_base.sage ecc2k130 --weight 2 --out ecc2k130-w2.json
sage low_hw_factor_base.sage p192 --weight 2 --out p192-w2.json
sage low_hw_factor_base.sage p224 --weight 2 --out p224-w2.json
sage low_hw_factor_base.sage p256 --weight 2 --out p256-w2.json
```

## `factor_base_dor.sage` — factor base to degree of regularity

Measures the solving degree `D*` (degree of regularity) of the
Weil-descended Semaev `S_3` system as a function of the factor base alone,
so you can search for a base that lowers it. On a non-decomposable target
the system is unsatisfiable and `D*` is the first Macaulay degree whose row
space contains `1`: the Polynomial-Calculus refutation degree, which is
what `RESEARCH_FFD_PROOF_COMPLEXITY.md` identifies as the security number.

The factor base is an `F_2`-subspace `V` of x-coordinates on the Koblitz
curve `y^2 + xy = x^3 + 1` over a toy field `F_{2^n}`, and each summand may
use its own subspace. Structures: `subfield`, `frobenius`
(squaring-stable), `standard` (`span{1, z, ..., z^{n'-1}}`), `random`, and
two normal-basis windows: `nb_window` puts both summands in
`V1 = <θ^(2^j) : j < n'>`, and `nb_break` is the symmetry breaking of
Galbraith–Granger–Merz–Petit (SAC 2020, §3.1): `P1` from `V1`, `P2` from
`V2 = π^{n'}(V1)`, with the relation folded back onto `F1` by Frobenius.
The report includes the fraction of targets that decompose, so the
relation-probability gain and the degree cost appear side by side.

```
sage factor_base_dor.sage --field 8 --dim 4 --selftest 8   # validate first
sage factor_base_dor.sage --field 8 --dim 4 --targets 48 --out dor.json
```

`--summands 3` uses the chained Semaev system (two `S_3` with a full-field
intermediate). Two findings there: symmetry breaking (`nb_break`) lifts the
decomposable fraction toward `3! = 6` (5.4x at n=11, dim 3), but the chained
`D*` is structure-independent (= 5 across subfield/standard/random), because
the full-field intermediate dilutes the factor base.

`--symmetric-s4` is the direct `m = 3` encoding instead: `S_4` written in
the elementary symmetric functions `e1, e2, e3` of the summands (the
formula of `src/cryptanalysis/koblitz_sym_s4.rs`, checked here against the
resultant), Boolean-quadratic after descent. Structures: `geometric`
(`V = θ<1, g, …, g^{l-1}>`, the Galbraith–Gebregiyorgis "special subspace",
`e_k` in a space of dimension `k(l-1)+1`), `standard` (`θ = 1, g = z`),
`subfield` (closed under products, `e_k ∈ V`), `random` (`e_k` in
`span(V^k)`). Satisfiability is labelled exactly with PolyBoRi, and the
report counts *spurious* targets: Boolean-satisfiable with no point
decomposition. Findings (`factor_base_dor_results.json`,
`symmetric_s4_direct`): with `6l - 3 > n` the system is underdetermined and
never unsatisfiable (n = 9, l = 3: 0/16 refutable on every structure), which
is why GG add membership constraints; at the determined cells `D*` for the
geometric window is 4 at n = 9, l = 2 and ≥ 5 at n = 15, l = 3 — no lower
than the chained system — while the subfield base refutes at `D* = 3` with
no spurious solutions (and no coverage at that size). The special subspace
buys unknowns (15 against the chain's 24 at n = 15), not degree.

```
sage factor_base_dor.sage --field 9 --dim 2 --symmetric-s4 --selftest 12
sage factor_base_dor.sage --field 15 --dim 3 --symmetric-s4 --dmax 4 --targets 16 --dstar-sample 6
```

`--delta-sweep D0,D1,...` is the Shantz–Teske delta method: the base
dimension is set by `m·dim = n − δ` and the sweep reports, per `δ`, the
decomposable fraction (beside the predicted `2^{-δ}` factor), the mean `D*`,
the Macaulay shape at `D*`, and a work-per-relation proxy
(`rows × cols / fraction`). At n = 8, 9 (m = 2) `δ = +1, +2` lower `D*` by
about one degree and the proxy 10–13x; `δ = +3` at n = 9 leaves no
decomposable target in 24, so the yield side wins once `2^{-δ}` reaches the
base's coverage. The proxy ignores that a smaller base also needs fewer
relations, so it understates the gain.

```
sage factor_base_dor.sage --field 9 --structures standard --delta-sweep=-1,1,3 --dmax 5 --targets 24
```

`--frobenius-windows` checks the GGMP §3.1 mechanism: the windows
`V_i = <θ^(2^j)>` satisfy `π^dim(V_i) = V_{i+1}` and are disjoint, so a
summand from `F_i` folds into `F_1` as `[λ^((i-1)·dim)] P'`. That yields the
`m!` saving with only `|F_1|` unknowns. It does not yield the `n`-fold
invariance saving: a window is not Frobenius-invariant, and on the prime
subgroup Frobenius acts as `λ`, so a relation's Frobenius images are scalar
multiples of it. The `n`-fold needs each base to be invariant itself (GGMP
Theorems 3.4 and 3.5). The check needs only a normal basis, so it holds at
`n = 131`:

```
sage factor_base_dor.sage --frobenius-windows --field 131 --dim 4 --summands 3
```

`--stacking A,B` asks whether both savings can coexist without an invariant
subspace. It uses the normal-basis weight classes
`W_w = {x : normal-basis weight w}`: Frobenius rotates normal-basis
coordinates, so each class is invariant, and distinct weights are disjoint.
Pairing `W_A` with `W_B` breaks symmetry while keeping the orbit saving. It
reports `trials = (unknowns + 1) / coverage` and an enumeration-work proxy
against the windows. At `n = 11` and `13` (same `x^n - 1` shape as 131) the
weight classes cut unknowns by exactly `n` and trials by about 10–20x, and
cross-class pairing lowers enumeration work a further 1.2–1.9x. Membership is
a cardinality constraint, so this suits the SAT oracle, not a linear descent.

```
sage factor_base_dor.sage --stacking 1,2 --field 11 --dim 4
```

Run `--selftest` before trusting a number. It checks that the descended
system has a Boolean zero exactly when brute force finds a decomposition,
and that a satisfiable target never derives `1`.

What it does not cover:

- **The Frobenius-invariant speedup itself.** At `n = 131` there is no
  Frobenius-stable *subspace* of useful size:
  `x^131 - 1 = (x+1)·(irreducible of degree 130)`, so the `subfield` and
  `frobenius` structures here do not exist there. Frobenius-invariant
  factor *bases* do: the weight-2 ONB base from `low_hw_factor_base.sage`
  (14 orbits) and the Couveignes–Lercier sets GGMP build. GGMP show such a
  base cuts the relations needed by about `n` and the linear algebra by
  about `n^2`, while their invariant vector spaces solve no slower than
  random ones. That saving is in the number of systems, not in `D*`, so
  this engine does not show it. The ONB weight base is not a subspace, so
  it cannot be descended linearly here (membership is a cardinality
  constraint; the repo's SAT oracle handles it).
- **`nb_break` does apply at 131.** It needs only a normal basis, which
  ECC2K-130 has. `--field 11` has the same `x^n - 1` factorization as 131.
- **Symmetric-function rewriting** (Faugère–Gaudry–Huot–Renault). Writing
  the summands in elementary symmetric variables drops binary `S_4` from
  total degree 12 to degree 4. Over a vector-space factor base this does not
  by itself lower the solving degree (Galbraith–Gebregiyorgis 2014; HKY need
  a "special subspace" closed under products, which is the subfield case).
  `--symmetric-s4` measures exactly this and agrees: the geometric special
  subspace does not lower `D*`, the multiplicatively closed base does. What
  the engine still lacks is GG's membership constraints, so the symmetric
  system is only well-posed when `6l − 3 ≤ n`.
- **Prime fields.** There is no subfield to descend into, so the degree of
  regularity is a different object (Diem's degree-`k` extension, or the
  symmetric Semaev system over `F_p`). `--curve p256` prints this note.
- **Total cost.** A lower `D*` is a win only if coverage holds. Shrinking
  `V` lowers `D*` by overdetermining the system and also lowers the chance
  a target decomposes. Compare structures at equal `n` and `n'`.

## Campaign 2026-10-13: factor base → D* yield and correlations

`research/factor_base_dor_campaign_20261013/` runs this engine over the
full structure × dimension grid (n = 9, 11, 13, 15; 3 seeds), the
Shantz–Teske δ sweep at n = 11 and 13, the symmetric-S4 cells at
(n = 9, l = 2) and (n = 15, l = 3), and the n = 131 applicability checks;
its `REPORT.typ`/`REPORT.pdf` carry the measured tables and figures.
Findings (all toy-scale stage diagnostics):

- **Only multiplicative closure lowers `D*`.** The subfield base refutes
  at `D* = 2` exactly at every matched cell (vs 2.35–3.00 standard/random),
  and at `D* = 3` with 9 unknowns in the symmetric-S4 encoding at
  n = 15, l = 3 against `D* = 5` with 15 unknowns for standard *and* the
  geometric special subspace. `nb_window`/`nb_break` never beat random on
  `D*`; the Frobenius-window fold and the weight-class stacking (10–38×
  trial reductions at n = 11/13) are unknowns levers, not degree levers.
- **`D*` tracks the base dimension, not the field degree**
  (Pearson r = 0.87 with `n'`, −0.04 with `n`), rising ≈ one degree per
  +2 of `n'`; the standard base is measurably below random (partial
  closure), e.g. 2.23 vs 2.79 at n = 13, `n' = 4`, all seeds.
- **Work per relation has an interior optimum in `n'`** (≈ `n/3..n/4` at
  toy scale): the largest measured cells are 20–300× worse per relation
  than the optimum (n = 13 standard: 1.6·10⁴ at `n' = 3` against
  4.5·10⁶ at `n' = 6`). The δ sweep improves per-relation work 6.5×
  (n = 11) and 11× (n = 13) from δ = −1 to δ = +3.
- **At n = 131 the closure lever is absent** (`x^131 − 1` has no proper
  factor), so the challenge curve only gets the size (δ) and arity levers;
  the Frobenius-window fold itself is verified at 131
  (`raw/frobwin_n131.log`). The m=83 confidence gate of AGENTS.md §8a
  is not yet run for any of these levers.

**Round 2 (scaling, same study directory):** the grid extended to n = 17,
19, 23 and the δ sweeps to n = 15, 17, 19 (n = 19 shares 131's full-length
cyclotomic order, ord = n − 1). The size laws scale: D* still tracks n′
(r = 0.82 over 133 aggregates), the interior work-per-relation optimum sits
at n′ ≈ n/2.7 with 100–875× spreads at n = 19, and the per-relation δ gain
*grows* with n (6.5× at n = 11 → 19.7× at n = 19, 246× at n = 15 δ = +5),
until the counting floor 2^−δ empties the yield (visible at n ≥ 17, δ ≥ 5).
The new arity+size result: on the symmetric-S4 m = 3 encoding, a heavy δ
drops D* by whole degrees — 5 → 3 from (n = 15, l = 3) to (n = 15, l = 2),
and D* = 4 at (n = 19, l = 3, δ = +10) with yield still 0.0625 — the one
whole-degree lever that exists at n = 131.

**Round 3 (the arity+size law):** seven more symmetric-S4 cells up to
n = 29 (ord₂(29) = 28 = n − 1, 131's full-length order). At heavy δ the law
is `D* = l + 1`, independent of the field degree: D* = 3 at l = 2 for
n = 19, 21, 23; D* = 4 at l = 3 for n = 19, 23, 29; D* = 5 at l = 4 for
n = 29 — the m = 3 oracle's per-system Macaulay cost grows only with the
base dimension, never with n. Measured point-decomposable yield sits far
above the naive counting floor at heavy δ (2^-5 measured against a 2^-20
floor at n = 29, l = 3 — Frobenius phases and permutation orbits amplify
yield ~2^12–2^15). One cell (n = 21, l = 4) was lost to a PolyBoRi
segfault, preserved in the logs. Drastic-speedup recipe at 131: an m ≥ 3
oracle with small l (D* ≤ 5) plus orbit-amplified yield; the m = 83 gate
is the unmet validation step.

**Round 4 (yield is governed by the e-space, not the field):**
high-precision yield cells (256 targets at l = 2 across n = 23..41;
128 at l = 3 across n = 29, 37; n = 37 has 131's full-length order
ord₂(37) = 36). The point-decomposable fraction is *constant in n*:
random 0.00781 and geometric 0.00391 at l = 2 at every n; random 0.0234
and geometric 0.0078 at l = 3 at both n = 29 and 37 — the amplification
A = pfrac·2ⁿ/2^(3l) is exactly 2^(n−13) (random l = 2) / 2^(n−14)
(geometric l = 2), i.e. the counting is against the symmetric e-space
(~2¹³–2¹⁴), not 2ⁿ. At n = 37, l = 3 the geometric base refutes at
D* = 3, *below* the l+1 diagonal, with yield unchanged. Derived
consequence (F15): per-relation oracle cost — systems per relation ×
Macaulay cost per system — is independent of the field degree; the
relation harvest scales with the factor-base size and its Frobenius-orbit
reduction, not with the subgroup size. Remaining obligations: the
orbit-reduced linear algebra and the m = 83 gate.

## `sm2.sage` — SM2 (sm2p256v1, GB/T 32918)

Mirrors `../src/ecc/sm2.rs` and `../src/hash/sm3.rs`. Run:

```
sage sm2.sage
```

What it does:

1. **Curve verification.** Builds the recommended SM2 curve and checks the
   structural / security properties:
   - special prime `p = 2²⁵⁶ − 2²²⁴ − 2⁹⁶ + 2⁶⁴ − 1` (Solinas form),
   - `a = p − 3`, cofactor `h = 1`, prime `n`, `#E = n`,
   - non-anomalous, no low embedding degree (MOV/Frey–Rück safe),
   - large CM discriminant,
   - twist structure: `#E' = 7 · 443 · p₁₀₀ · p₁₄₈` — **not** prime-order.
     Small-subgroup part is only `3101`, but the 100-bit factor means
     invalid-curve attacks recover ~100 bits at ~2⁵⁰ work *if* point
     validation is skipped → SM2 implementations must validate input points.

2. **Faithful signature scheme.** Real SM3 (verified against the canonical
   `SM3("abc")` vector, which pins it to the same function as the Rust impl)
   plus the `ZA = SM3(ENTL ‖ ID ‖ a ‖ b ‖ Gₓ ‖ G_y ‖ PAₓ ‖ PA_y)` identity
   binding and the SM2 sign/verify algebra `s = (1+d)⁻¹(k − rd)`.

3. **SM2-specific HNP attack.** The repo's `hnp_ecdsa.rs` only covers ECDSA.
   SM2's equation rearranges to

   ```
   k = s + (r + s)·d   (mod n)
   ```

   i.e. an HNP instance with constant `tᵢ = sᵢ` and multiplier `Aᵢ = rᵢ + sᵢ`.
   **The message digest `e` cancels out** — unlike ECDSA, recovery needs only
   the `(r, s)` pairs and the nonce bias, not the messages. With the top 16
   bits of each nonce zeroed, 40 signatures + LLL recover the private key
   (Boneh–Venkatesan, same scaled-integer basis as `hnp_ecdsa.rs`).
