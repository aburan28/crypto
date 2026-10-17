# Cost model: Joux–Vitse `(n − 1)`-decomposition on `E(F_{qⁿ})`, `n = 5, 6`

Written 2026-10-15 before the measured runs (`PROTOCOL.md`); the measured
inputs are filled in by `README.md` from `summary.json`.  Statements are
labelled **[paper]** (quoted from a source), **[derived]** (a calculation from
stated inputs), **[measured]** (this study's records) or **[extrapolation]**.

## Sources read

- V. Vitse, *Attaques algébriques du problème du logarithme discret sur
  courbes elliptiques*, PhD thesis, UVSQ, defended 2011-10-20
  (<https://www-fourier.univ-grenoble-alpes.fr/~viva/research/PhD/VitsePhD.pdf>):
  §3.3.1 and Table 3.1 (F4 variants on the `E(F_{p⁵})` systems), §7.2.2
  (Gaudry–Diem), §7.3.1–7.3.4 (the `n − 1` variant, its analysis, the
  `F_{q⁵}` application, oracle-assisted SDHP), chapter 8 (cover and
  decomposition, `F_{p⁶}`).
- Joux–Vitse, *Elliptic curve discrete logarithm problem over small degree
  extension fields*, J. Cryptology 26 (2013), ePrint 2010/157: the `n − 1`
  variant (the thesis's §7.3 is its text).  Joux–Vitse, *Cover and
  decomposition index calculus on elliptic curves made practical*,
  Eurocrypt 2012, ePrint 2011/020: the `F_{p⁶}` result.  The ePrint and
  IACR-archive PDFs were not retrievable from this container (HTTP 403);
  the thesis carries both papers' content and is what is quoted.
- Gaudry 2009 and Diem 2011 through the thesis's §7.2 summary (not read
  directly).

**[paper] Correction to the premise.**  The practical `F_{p⁶}` result
(a 132-bit curve in about 3,700 CPU·hours, thesis ch. 8; reported as 149
bits / 110,000 CPU·hours in the Eurocrypt version) is **cover and
decomposition**: the curve is transferred to the Jacobian of a genus-2 or
genus-3 curve over an intermediate field and decomposed there (Nagao).  It
is not `(n − 1)`-decomposition on `E(F_{p⁶})`.  The practical `F_{q⁵}`
result for the `n − 1` variant is the **oracle-assisted static
Diffie–Hellman** attack on the Oakley curve over `F_{2^{155}} = F_{(2^{31})⁵}`
(§7.3.4: a `4!·2^{31}`-trial search at 22.95 ms a trial), **not a complete
DLP**.  For the DLP, §7.3.3 says: "cette variante reste encore trop lente
pour pouvoir résoudre en temps raisonnable le problème ECDLP sur des corps
de tailles compatibles avec les niveaux de sécurité exigés actuellement".

## The method at `n = 5`

Factor base `F = {P ∈ E(F_{q⁵}) : x(P) ∈ F_q}/±`, `|F| ≈ q/2`.  A trial
takes `R = aP + bQ` and asks whether `R = ±P₁ ± P₂ ± P₃ ± P₄`, `P_i ∈ F`:
`S₅(x₁, …, x₄, x_R) = 0`, symmetrised in `e₁, …, e₄` (total degree 8) and
Weil-restricted: **5 equations of degree 8 in 4 unknowns over `F_q`**
(Vitse (7.9)).

| quantity | value | label |
|:--|:--|:--|
| trials per relation | `(n − 1)! q = 24 q` (Hypothesis 7.3.2) | [paper] |
| counting bound | rate `≤ 16·C(|F| + 3, 4)/#E ≈ (2/3)/q` | [derived] |
| relations needed | `|F| + 1 ≈ q/2` | [paper] |
| relation search | `O((n − 1)! c̃ q²)`; with the above, `≈ 12 q² · c̃` | [paper]/[derived] |
| `c̃`, one trial | `≈ 3·10⁹` `F_p` mults (§7.3.3); F4 `10¹⁰`, F4Remake ~3× fewer (Rem. 7.3.6); 8.4 s/trial at 32-bit `p` (Table 3.1) | [paper] |
| Macaulay bound | `Σ(deg ϕ_i − 1) + 1 = 5·7 + 1 = 36` (Hyp. 7.3.4); 29 F4 steps | [paper] |
| linear algebra | `Õ(n q²)` sparse, `≈ 4 (q/2)²` mults mod `ℓ`, negligible | [paper] |
| rho | `√(π·q⁵/4)` group operations (negation map), `c_add` mults each | [derived] |
| crossover | "plus rapide dès que la cardinalité du corps de base est plus grande que `2^{60}`" (odd `p`), `2^{45}` in characteristic 2 | [paper] |

**[derived] Crossover from the closed forms.**
`IC ≈ (q/2)·κq·c̃` against `rho ≈ √(π/4)·q^{5/2}·c_add`, so
`q* = (κ c̃ / (2·√(π/4)·c_add))²`.  With the paper's `κ = 24`,
`c̃ = 3·10⁹`, and `c_add ≈ 190` (this repository's affine `F_{p⁵}` addition,
measured below): `q* ≈ (24·3·10⁹ / 337)² ≈ 4.6·10¹⁶ ≈ 2^{55.3}`, group
`≈ 2^{276}`.  The paper's `2^{60}` is the same order (it prices rho
differently).  At any `q` this repository can run (`q ≤ 2^{12}`) the
predicted ratio is `IC/rho ≈ √(q*/q) ≳ 10⁶`.

**Why large primes do not help here [derived].**  The relation search
dominates (`≈ 12 q² c̃` against `≈ q²` for the linear algebra), so there is
no linear algebra to rebalance.  Shrinking the factor base by a factor `λ`
multiplies the trials per relation by `λ⁴` (four base points per
relation), and large-prime relations cost a full trial each: the trade goes
the wrong way.

**Comparison with `n = 4` [derived, different groups].**  Gaudry's square
`n = 4` system needs about `4! = 24` trials per relation at any `q`; with
PR #1863's D₄ oracle (`C₄ = 5.5–8.5·10⁷`, open and unmerged) that is
`≈ 1.7·10⁹` per relation.  The `n − 1` variant at `n = 5` pays
`24 q · c̃ ≈ 1.3·10¹¹ q` per relation: about `80 q` times more, on a group
`q` times larger (`q⁵` against `q⁴`), whose rho is `√q` times dearer.

## `n = 6` [derived; not measured]

`S₆` symmetrised in `e₁, …, e₅` has total degree 16 on `C(21, 5) = 20,349`
monomials, times 17 powers of `x_R`.  Interpolating it costs about
`(20,349/495)³ ≈ 6.9·10⁴` times the measured `n = 5` precomputation
(`1.87·10⁹`), so about `1.3·10¹⁴` `F_p` multiplications, on a dense
20,349 × 20,366 system over `F_{p⁶}` (about 20 GB at 6 words an entry;
this container has 15 GB).  The trial system is 6 equations of degree 16
in 5 unknowns, Macaulay bound `6·15 + 1 = 91`, with `C(96, 5) = 6.1·10⁷`
monomials up to that degree (`C(40, 4) = 91,390` at `n = 5`, where F4's
widest matrix is about 3,700 columns).  Trials per relation would be
`5!·q = 120 q`.  Neither the symmetrised `S₆` nor one trial is within this
container's memory, so `n = 6` was not attempted.  This is a resource
limit, not a mathematical obstruction.
