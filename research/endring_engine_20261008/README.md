# endring-engine: curve-side layers at toy parameters

Status: **work in progress, curve side only.**  Standalone research crate
(own `Cargo.toml`, no dependency on the root crate).  Not a cryptosystem,
not constant-time, no security claim.

What is here, each with passing tests:

- `src/fp.rs` — `F_p` and `F_{p²} = F_p(i)` for `p ≡ 3 (mod 4)`, `p < 2⁶³`,
  with square roots and Frobenius.
- `src/curve.rs` — short Weierstrass curves over `F_{p²}`, full `(x, y)`
  arithmetic, torsion bases of `E[N]` for `N | p + 1`, exact
  two-dimensional discrete logarithms in smooth torsion.
- `src/isogeny.rs` — Vélu's formulas with full `(x, y)` images (exact
  homomorphisms, exact kernels) and smooth-degree chains.

Toy parameter used by the tests: `p = 2⁵⁴·3²·5 − 1`, where
`E₀: y² = x³ + x` has `E₀(F_{p²}) ≅ (Z/(p+1))²`.

Not here: the quaternion order, its left ideals, the Deuring
correspondence and KLPT.

```sh
cargo test --manifest-path research/endring_engine_20261008/Cargo.toml
```
