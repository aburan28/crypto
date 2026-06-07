# SageMath scaffolding

Sage companions to the Rust crypto library. Each file is self-contained and
runnable with `sage <file>.sage`.

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
