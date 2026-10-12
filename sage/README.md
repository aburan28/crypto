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
  total degree 12 to degree 4. The engine descends the raw summand
  coordinates, so it does not yet measure this lever.
- **Prime fields.** There is no subfield to descend into, so the degree of
  regularity is a different object (Diem's degree-`k` extension, or the
  symmetric Semaev system over `F_p`). `--curve p256` prints this note.
- **Total cost.** A lower `D*` is a win only if coverage holds. Shrinking
  `V` lowers `D*` by overdetermining the system and also lowers the chance
  a target decomposes. Compare structures at equal `n` and `n'`.

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
