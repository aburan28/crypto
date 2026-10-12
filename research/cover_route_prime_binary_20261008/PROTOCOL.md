# The cover route for prime-field and binary curves: protocol (registered 2026-10-08, before any run)

**Question.**  The Joux–Vitse / GHS cover route transfers an ECDLP on
`E / F_{q^k}` to the Jacobian of a curve over the subfield `F_q`.  Does it
have any target among deployed prime-field and binary curves, and what is
its cost on the only composite-degree binary curves ever standardised?

## Part A: exclusion by structure (theory, no run)

- A prime field `F_p` has no proper subfield.  A cover `C → E` with `C`
  over `F_p` itself only enlarges the group: index calculus on `Jac(C)`
  of genus `g` over `F_p` costs `Õ(p^{2 − 2/g})`, which exceeds rho's
  `√p` for every `g ≥ 2`; base-changing `E` to `F_{p^k}` and covering
  from `F_p` gives `g ≥ k` and `Õ(p^{2 − 2/g}) ≥ p^{4/3}` at `k = 3`.
  **Prediction A-1:** no deployed prime-field curve (the NIST primes,
  secp256k1, Curve25519 and Curve448) has a cover-route cost below rho,
  for any `g`.  This is a statement about the complexity formulas, not a
  measurement.
- A binary field `F_{2^n}` with `n` prime has the single proper subfield
  `F_2`.  GHS over `F_2` has magic number `m ∈ {1, 2}` when `b ∈ F_2` (the
  Koblitz curves) and `m ∈ {n, n + 1}` otherwise, because the `F_2`-degree
  of `√b` divides `n`.  Genus `≤ 2` cannot carry a subgroup of order
  `≈ 2^{n−2}`; genus `2^{n−1}` admits no index calculus.
  **Prediction A-2:** every binary curve in FIPS 186-4 / SEC 2 (`n = 113,
  131, 163, 233, 239, 283, 409, 571`) has prime `n`, so the route is empty
  for all of them; the Galbraith–Hess–Smart isogeny extension cannot
  change `m` because `m` depends only on the subfield degree of `√b` and
  `b ∈ F_2` is the only small case.
- Pairing G2 groups live on twists over `F_{p^2}` (BN254, BLS12-381) or
  `F_{p^4}`; the subgroup has `r ≈ p` or less, and Gaudry's `k = 2` index
  calculus costs `Õ(p)` against rho's `√r ≤ √p`.  **Prediction A-3:** no
  deployed pairing G2 is a target either.

## Part B: the composite-degree binary curves (measured)

RFC 2409 §6.3–6.4 standardise two EC2N groups: Oakley group 3 over
`F_{2^155}`, `u^155 + u^62 + 1`, `a = 0`, `b = 0x7338f`; and Oakley group 4
over `F_{2^185}`, `u^185 + u^69 + 1`, `a = 0`, `b = 0x1ee9`; group orders
as printed there.  `155 = 5 · 31` and `185 = 5 · 37`, so each admits the
GHS descent to `F_{2^31}` or `F_{2^37}` with `n′ = 5`, and the reverse
descent to `F_{2^5}` with `n′ = 31` or `37`.

The GHS magic number for `E: y² + xy = x³ + ax² + b` over `K = F_{q^{n′}}`
is `m = dim_{F_2} Span{(1, σ^i(γ)) : 0 ≤ i < n′}` with `γ = √b` and `σ` the
`q`-Frobenius; the cover has genus `2^{m−1}` or `2^{m−1} − 1`.  The tool
`ghs_magic.rs` computes `m` exactly for the four (curve, subfield) pairs
and reports whether `b` lies in a proper subfield.

**Predictions, registered now.**

- **B-1.** With `q = 2^31` and `q = 2^37` (`n′ = 5`): `m = 5` for both
  curves, genus 15 or 16.  *Falsified if* `m ≤ 4` for either.
- **B-2.** With `q = 2^5` (`n′ = 31` or `37`): `m = n′` or `n′ + 1`,
  genus `≥ 2^30`; the small-base descent is useless.  *Falsified if*
  `m ≤ 20`.
- **B-3 (subgroup carried).**  `g · log₂ q ≥ log₂ r` for the `n′ = 5`
  descents, so the prime-order subgroup survives the conorm–norm map.
- **B-4 (pricing, formulas only).**  With Gaudry–Thériault–Thomé–Diem
  double-large-prime index calculus at `Õ(q^{2 − 2/g})` field operations
  and rho at `√(πr/4)` group operations, the route is at least `2^8`
  below rho for both groups.  *Falsified if* the formula ratio is below
  `2^8`.  This is a literature-formula row and is labelled so; no
  hyperelliptic index calculus is run here.

**Class if all hold:** boundary for Part A, boundary-with-a-formula for
Part B.  Nothing in this protocol touches a curve in a deployed hybrid
scheme; the two Oakley groups were deprecated with IKEv1.
