# The cover route for prime-field and binary curves: results (2026-10-08)

Protocol: `PROTOCOL.md`, registered before any run.  Tool: `ghs_magic.rs`
(standalone, `rustc -O`; `./ghs_magic [deployed_binary.txt]`), output
`results_ghs_magic.txt`.  Curve parameters for the deployed binary curves
come from `docs/curves/registry.json` (`deployed_binary.txt` lists name,
degree, modulus, `b`, slug); the RFC 2409 parameters are transcribed from
the RFC's §6.3 and §6.4.

## Part A: deployed curves have no subfield the route can use

**Prime fields.**  `F_p` has no proper subfield.  A cover `C → E` with `C`
over `F_p` enlarges the group, and index calculus on `Jac(C)` of genus `g`
over `F_p` costs `Õ(p^{2 − 2/g}) ≥ p` against rho's `√p`; base-changing to
`F_{p^k}` and descending to `F_p` gives `g ≥ k` and at best `p^{4/3}`.
**A-1 holds** for every deployed prime-field curve, as a statement about
the complexity formulas.  Pairing G2 groups (twists over `F_{p^2}`,
`F_{p^4}`) face Gaudry's `Õ(p)` against `√r ≤ √p`: **A-3 holds**.

**Binary fields of prime degree.**  The only proper subfield is `F_2`.  The
GHS magic number is `m = dim_{F_2} Span{(1, γ^{2^i})}` with `γ = √b`; the
span of the conjugates of `γ` is the cyclic `F_2[σ]`-module it generates,
whose dimension is the degree of its annihilator, a divisor of `x^n − 1`.
For prime `n` the irreducible factors of `x^n − 1` over `F_2` are `x + 1`
and factors of degree `ord_n(2)`, so `m ∈ {1, 2}` or `m ≥ ord_n(2)`.  The
protocol's phrasing "`m ∈ {n, n + 1}` otherwise" was too strong as a
theorem and is corrected here; the measured values below are the maximum
anyway.

| n | ord_n(2) | smallest genus above the degenerate case |
|--:|--:|--:|
| 113 | 28 | 2^27 |
| 131 | 130 | 2^129 |
| 163 | 162 | 2^161 |
| 233 | 29 | 2^28 |
| 239 | 119 | 2^118 |
| 283 | 94 | 2^93 |
| 409 | 204 | 2^203 |
| 571 | 114 | 2^113 |

Measured over `F_2` for every binary curve in FIPS 186-4 / SEC 2:

| curve | b in a proper subfield | m | cover genus |
|:--|:--|--:|--:|
| sect163k1, sect233k1, sect239k1, sect283k1, sect409k1, sect571k1 | yes, `b = 1` | 1 | 1 or 0: the subgroup is not carried |
| sect113r1 | no | 113 | 2^112 |
| sect131r1 | no | 131 | 2^130 |
| sect163r2 | no | 163 | 2^162 |
| sect233r1 | no | 233 | 2^232 |
| sect283r1 | no | 283 | 2^282 |
| sect409r1 | no | 409 | 2^408 |
| sect571r1 | no | 571 | 2^570 |

Every Koblitz curve degenerates (`m = 1`) and every random-`b` curve sits
at the maximum `m = n`, nowhere near the `ord_n(2)` floor that a
specially chosen `b` could reach at `n = 113` or `233`.  The
Galbraith–Hess–Smart isogeny extension does not help: `m` depends on the
annihilator of `√b`, and the only small case is `b ∈ F_2`, which is the
degenerate one.  **A-2 holds**, with the corrected theory.

## Part B: the two composite-degree binary curves ever standardised

RFC 2409's Oakley groups 3 (`F_{2^155}`, `155 = 5 · 31`) and 4
(`F_{2^185}`, `185 = 5 · 37`), both `a = 0`.  Group orders from the RFC:
group 3 has cofactor 12 and a 152-bit prime `r`; group 4 has cofactor 4
and a 183-bit prime `r` (trial division to `2^22`, residual a strong
probable prime to base 2).

| curve | base | n′ | m | cover genus |
|:--|:--|--:|--:|--:|
| Oakley group 3 | F_2 | 155 | 126 | 2^125 |
| Oakley group 3 | F_{2^5} | 31 | 26 | 2^25 |
| Oakley group 3 | F_{2^31} | 5 | 5 | 16 or 15 |
| Oakley group 4 | F_2 | 185 | 185 | 2^184 |
| Oakley group 4 | F_{2^5} | 37 | 37 | 2^36 |
| Oakley group 4 | F_{2^37} | 5 | 5 | 16 or 15 |

- **B-1 holds:** `m = 5`, genus 15 or 16, for both curves over the large
  subfield.  Since `m ≤ n′ = 5` for every curve over these fields, every
  curve over `F_{2^155}` and `F_{2^185}` has such a cover; the group's
  choice of `b` makes no difference here.
- **B-2:** the kill condition (`m ≤ 20`) is not met, but the registered
  value was wrong for group 3: `m = 26`, not 31 or 32, because `√b`'s
  annihilator over `F_{2^5}` omits a degree-5 factor.  Recorded as a
  prediction that was right in its consequence and wrong in its number.
  Group 4 gives `m = 37` as stated.
- **B-3 holds:** `16 · 31 = 496 ≥ 152` and `16 · 37 = 592 ≥ 183`, so the
  prime-order subgroup survives the conorm–norm map to `Jac(C)(F_q)`.
- **B-4 (formulas only).**  Double-large-prime hyperelliptic index
  calculus at `Õ(q^{2 − 2/g})` with `g = 16`:

| curve | q | IC formula, field ops | rho `√(πr/4)`, group ops | formula ratio |
|:--|--:|--:|--:|--:|
| Oakley group 3 | 2^31 | 2^58.1 | 2^75.8 | 2^17.7 |
| Oakley group 4 | 2^37 | 2^69.4 | 2^91.3 | 2^22 |

  Both exceed the registered `2^8`.  The units differ (a group operation
  is roughly ten field multiplications, which favours the cover side
  further), the `Õ` hides constants that grow with `g`, and no
  hyperelliptic index calculus at genus 16 exists in this repository; the
  row is literature formulas, labelled so.  It agrees with the published
  reading of these fields as weak for the GHS route, which is why they
  were never carried into IKEv2.

## What this settles and what it leaves

- The cover route has **no target** among deployed prime-field curves,
  deployed prime-degree binary curves, or pairing G2 groups, by the
  structure of the fields and the measured `m`.  This closes the user's
  question for the curves inside any hybrid scheme.
- The only standardised curves it reaches are RFC 2409's two EC2N groups,
  deprecated with IKEv1, where it beats rho by `2^18` to `2^22` in formula
  terms.  Turning that into a measured row would need a genus-16
  hyperelliptic index calculus over `F_{2^31}`, a self-contained
  engineering project with a known answer; it is not started.
- The one open mathematical thread the earlier census left, which classes
  over `F_{q^3}` have an *explicit* small-degree cover, is untouched by
  anything here and remains the live part of the cover line.

**Class: boundary** throughout.  No speedup is claimed.
